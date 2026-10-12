#!/usr/bin/env python3
"""
Generate fully dynamic synthetic graph streams for D-ESCENT.

Output format:
    + u v timestamp
    - u v timestamp

Supported graph models:
    er       Erdos-Renyi-like sparse graph
    ba       Barabasi-Albert preferential attachment
    grid     2D torus/grid-with-diagonals
    regular  random d-regular-like graph

Supported deletion modes:
    random   delete random non-protected edges
    hub      delete edges incident to high-degree vertices
    tree     delete edges from a maintained seed spanning tree when safe
    mixed    rotate random, hub, and tree deletions
"""

from __future__ import annotations

import argparse
import math
import random
from collections import defaultdict
from pathlib import Path


Edge = tuple[int, int]
Event = tuple[str, int, int, int]


def edge_key(u: int, v: int) -> Edge:
    return (u, v) if u < v else (v, u)


def parse_sizes(raw: str) -> list[int]:
    return [int(x.strip()) for x in raw.split(",") if x.strip()]


class DynamicBuilder:
    def __init__(self, seed: int) -> None:
        self.rng = random.Random(seed)
        self.events: list[Event] = []
        self.active: set[Edge] = set()
        self.protected: set[Edge] = set()
        self.seed_tree: set[Edge] = set()
        self.deg: defaultdict[int, int] = defaultdict(int)

    def add_event(self, op: str, e: Edge) -> None:
        self.events.append((op, e[0], e[1], len(self.events)))

    def insert(self, u: int, v: int, protect: bool = False, tree: bool = False) -> bool:
        if u == v:
            return False
        e = edge_key(u, v)
        if e in self.active:
            return False
        self.active.add(e)
        self.deg[e[0]] += 1
        self.deg[e[1]] += 1
        if protect:
            self.protected.add(e)
        if tree:
            self.seed_tree.add(e)
        self.add_event("+", e)
        return True

    def delete(self, e: Edge) -> bool:
        if e not in self.active or e in self.protected:
            return False
        self.active.remove(e)
        self.deg[e[0]] -= 1
        self.deg[e[1]] -= 1
        self.add_event("-", e)
        return True

    def random_missing_edge(self, n: int) -> Edge:
        for _ in range(5000):
            u = self.rng.randrange(n)
            v = self.rng.randrange(n)
            if u == v:
                continue
            e = edge_key(u, v)
            if e not in self.active:
                return e
        raise RuntimeError("could not find a missing edge; graph may be too dense")

    def random_existing_edge(self, mode: str) -> Edge | None:
        candidates = [e for e in self.active if e not in self.protected]
        if not candidates:
            return None

        if mode == "random":
            return self.rng.choice(candidates)

        if mode == "hub":
            sample = self.rng.sample(candidates, min(len(candidates), 256))
            return max(sample, key=lambda e: self.deg[e[0]] + self.deg[e[1]])

        if mode == "tree":
            tree_candidates = [e for e in self.seed_tree if e in self.active and e not in self.protected]
            if tree_candidates:
                return self.rng.choice(tree_candidates)
            return self.rng.choice(candidates)

        raise ValueError(f"unknown deletion mode: {mode}")


def seed_cycle_with_chords(b: DynamicBuilder, n: int) -> None:
    for i in range(n):
        b.insert(i, (i + 1) % n, protect=True, tree=True)
    for i in range(n):
        b.insert(i, (i + 2) % n, protect=True)


def build_er_seed(b: DynamicBuilder, n: int, avg_degree: int) -> None:
    seed_cycle_with_chords(b, n)
    target = max(n * avg_degree // 2, 2 * n)
    while len(b.active) < target:
        b.insert(*b.random_missing_edge(n))


def build_regular_seed(b: DynamicBuilder, n: int, avg_degree: int) -> None:
    seed_cycle_with_chords(b, n)
    d = max(4, avg_degree)
    for offset in range(3, d // 2 + 2):
        for u in range(n):
            b.insert(u, (u + offset) % n)
    target = max(n * avg_degree // 2, len(b.active))
    while len(b.active) < target:
        b.insert(*b.random_missing_edge(n))


def build_grid_seed(b: DynamicBuilder, n: int, avg_degree: int) -> int:
    side = int(math.sqrt(n))
    n_grid = side * side
    for r in range(side):
        for c in range(side):
            u = r * side + c
            right = r * side + ((c + 1) % side)
            down = ((r + 1) % side) * side + c
            diag = ((r + 1) % side) * side + ((c + 1) % side)
            b.insert(u, right, protect=True, tree=True)
            b.insert(u, down, protect=True)
            b.insert(u, diag)
    target = max(n_grid * avg_degree // 2, len(b.active))
    while len(b.active) < target:
        b.insert(*b.random_missing_edge(n_grid))
    return n_grid


def build_ba_seed(b: DynamicBuilder, n: int, avg_degree: int) -> None:
    m = max(2, avg_degree // 2)
    init = max(8, m + 2)
    seed_cycle_with_chords(b, init)
    repeated: list[int] = []
    for e in b.active:
        repeated.extend(e)

    for v in range(init, n):
        targets: set[int] = set()
        while len(targets) < m:
            targets.add(b.rng.choice(repeated) if repeated else b.rng.randrange(v))
        first = True
        for u in targets:
            if b.insert(v, u, protect=first, tree=first):
                repeated.extend([v, u])
                first = False


def build_seed(model: str, b: DynamicBuilder, n: int, avg_degree: int) -> int:
    if model == "er":
        build_er_seed(b, n, avg_degree)
        return n
    if model == "ba":
        build_ba_seed(b, n, avg_degree)
        return n
    if model == "grid":
        return build_grid_seed(b, n, avg_degree)
    if model == "regular":
        build_regular_seed(b, n, avg_degree)
        return n
    raise ValueError(f"unknown model: {model}")


def generate_stream(
    model: str,
    final_vertices: int,
    avg_degree: int,
    updates_per_vertex: int,
    delete_rate: float,
    delete_mode: str,
    seed: int,
) -> list[Event]:
    b = DynamicBuilder(seed)
    n = build_seed(model, b, final_vertices, avg_degree)
    init_event_count = len(b.events)
    total_updates = n * updates_per_vertex
    target_edges = max(n * avg_degree // 2, len(b.active))
    modes = ["random", "hub", "tree"] if delete_mode == "mixed" else [delete_mode]

    for step in range(total_updates):
        mode = modes[step % len(modes)]
        too_dense = len(b.active) > target_edges
        do_delete = too_dense or (b.rng.random() < delete_rate and len(b.active) > len(b.protected) + 16)
        if do_delete:
            e = b.random_existing_edge(mode)
            if e is not None and b.delete(e):
                continue
        try:
            b.insert(*b.random_missing_edge(n))
        except RuntimeError:
            e = b.random_existing_edge("random")
            if e is not None:
                b.delete(e)

    # Keep the initial construction in the first half of the file so the C++
    # benchmark starts from a meaningful active graph.
    if init_event_count > len(b.events) // 2:
        padding = init_event_count - len(b.events) // 2
        for _ in range(padding):
            try:
                b.insert(*b.random_missing_edge(n))
            except RuntimeError:
                break
    return b.events


def write_stream(path: Path, events: list[Event], model: str, delete_mode: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("# fully dynamic synthetic graph stream\n")
        f.write(f"# model={model} delete_mode={delete_mode}\n")
        f.write("# format: op u v timestamp, where op is + insertion or - deletion\n")
        for op, u, v, ts in events:
            f.write(f"{op} {u} {v} {ts}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default="datasets/synthetic")
    parser.add_argument("--models", default="er,ba,grid,regular")
    parser.add_argument("--delete-modes", default="random,hub,tree")
    parser.add_argument("--sizes", default="1000,2000,3000")
    parser.add_argument("--avg-degree", type=int, default=8)
    parser.add_argument("--updates-per-vertex", type=int, default=5)
    parser.add_argument("--delete-rate", type=float, default=0.35)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    models = [x.strip() for x in args.models.split(",") if x.strip()]
    delete_modes = [x.strip() for x in args.delete_modes.split(",") if x.strip()]
    sizes = parse_sizes(args.sizes)
    out_dir = Path(args.out_dir)

    job_id = 0
    for model in models:
        for delete_mode in delete_modes:
            for n in sizes:
                job_id += 1
                events = generate_stream(
                    model=model,
                    final_vertices=n,
                    avg_degree=args.avg_degree,
                    updates_per_vertex=args.updates_per_vertex,
                    delete_rate=args.delete_rate,
                    delete_mode=delete_mode,
                    seed=args.seed + job_id,
                )
                path = out_dir / f"{model}_{delete_mode}_n{n}.txt"
                write_stream(path, events, model, delete_mode)
                ins = sum(1 for e in events if e[0] == "+")
                dele = len(events) - ins
                print(f"{path}: {len(events)} events ({ins} inserts, {dele} deletes)")


if __name__ == "__main__":
    main()
