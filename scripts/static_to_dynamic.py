#!/usr/bin/env python3
"""
Convert a static edge list into a dynamic insertion/deletion stream.

Input accepts whitespace-separated lines:
    u v
    u v weight
    u v timestamp

Output:
    + u v t
    - u v t
"""

from __future__ import annotations

import argparse
import random
from collections import defaultdict, deque
from pathlib import Path


Edge = tuple[int, int]
Event = tuple[str, int, int, int]


def edge_key(u: int, v: int) -> Edge:
    return (u, v) if u < v else (v, u)


def read_edges(path: Path) -> list[Edge]:
    edges: list[Edge] = []
    seen: set[Edge] = set()
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line[0] in "#%":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                u, v = int(parts[0]), int(parts[1])
            except ValueError:
                continue
            if u == v:
                continue
            e = edge_key(u, v)
            if e not in seen:
                seen.add(e)
                edges.append(e)
    return edges


def largest_component(edges: list[Edge]) -> list[Edge]:
    adj: defaultdict[int, list[int]] = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    seen: set[int] = set()
    best: set[int] = set()
    for start in adj:
        if start in seen:
            continue
        comp: set[int] = set()
        q = deque([start])
        seen.add(start)
        while q:
            u = q.popleft()
            comp.add(u)
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    q.append(v)
        if len(comp) > len(best):
            best = comp
    return [e for e in edges if e[0] in best and e[1] in best]


def degree_scores(edges: set[Edge]) -> dict[Edge, int]:
    deg: defaultdict[int, int] = defaultdict(int)
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    return {e: deg[e[0]] + deg[e[1]] for e in edges}


def bfs_tree_edges(edges: set[Edge]) -> set[Edge]:
    adj: defaultdict[int, list[int]] = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    if not adj:
        return set()
    root = next(iter(adj))
    seen = {root}
    q = deque([root])
    tree: set[Edge] = set()
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v in seen:
                continue
            seen.add(v)
            q.append(v)
            tree.add(edge_key(u, v))
    return tree


def choose_delete_edge(rng: random.Random, active: set[Edge], mode: str) -> Edge | None:
    if not active:
        return None
    if mode == "random":
        return rng.choice(tuple(active))
    if mode == "hub":
        scores = degree_scores(active)
        sample = rng.sample(tuple(active), min(len(active), 512))
        return max(sample, key=lambda e: scores[e])
    if mode == "tree":
        tree = bfs_tree_edges(active)
        if tree:
            return rng.choice(tuple(tree))
        return rng.choice(tuple(active))
    if mode == "mixed":
        return choose_delete_edge(rng, active, rng.choice(["random", "hub", "tree"]))
    raise ValueError(f"unknown delete mode: {mode}")


def random_missing_edge(rng: random.Random, vertices: list[int], active: set[Edge]) -> Edge | None:
    n = len(vertices)
    for _ in range(5000):
        u = vertices[rng.randrange(n)]
        v = vertices[rng.randrange(n)]
        if u == v:
            continue
        e = edge_key(u, v)
        if e not in active:
            return e
    return None


def make_stream(
    edges: list[Edge],
    delete_mode: str,
    delete_rate: float,
    update_factor: float,
    seed: int,
) -> list[Event]:
    rng = random.Random(seed)
    edges = largest_component(edges)
    rng.shuffle(edges)
    init_count = max(4, len(edges) // 2)
    initial = edges[:init_count]
    reserve = edges[init_count:]
    active = set(initial)
    vertices = sorted({x for e in edges for x in e})
    events: list[Event] = []

    def add(op: str, e: Edge) -> None:
        events.append((op, e[0], e[1], len(events)))

    for e in initial:
        add("+", e)

    updates = max(1, int(len(edges) * update_factor))
    reserve_idx = 0
    for _ in range(updates):
        if rng.random() < delete_rate and len(active) > 8:
            e = choose_delete_edge(rng, active, delete_mode)
            if e is not None:
                active.remove(e)
                add("-", e)
                continue

        if reserve_idx < len(reserve):
            e = reserve[reserve_idx]
            reserve_idx += 1
        else:
            e = random_missing_edge(rng, vertices, active)
            if e is None:
                continue
        if e not in active:
            active.add(e)
            add("+", e)

    return events


def write_stream(path: Path, events: list[Event]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("# static graph converted to fully dynamic stream\n")
        f.write("# format: op u v timestamp\n")
        for op, u, v, ts in events:
            f.write(f"{op} {u} {v} {ts}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--out", required=True)
    parser.add_argument("--delete-mode", choices=["random", "hub", "tree", "mixed"], default="random")
    parser.add_argument("--delete-rate", type=float, default=0.35)
    parser.add_argument("--update-factor", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    edges = read_edges(Path(args.input))
    events = make_stream(edges, args.delete_mode, args.delete_rate, args.update_factor, args.seed)
    write_stream(Path(args.out), events)
    ins = sum(1 for e in events if e[0] == "+")
    dele = len(events) - ins
    print(f"{args.out}: {len(events)} events ({ins} inserts, {dele} deletes)")


if __name__ == "__main__":
    main()
