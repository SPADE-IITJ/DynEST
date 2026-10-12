#!/usr/bin/env python3
"""
Infer insertions/deletions from a timestamped temporal edge list.

For consecutive windows:
    G_i     = edges with timestamp in [start_i, start_i + window)
    G_{i+1} = edges with timestamp in [start_i + stride, start_i + stride + window)

Inserted edges = G_{i+1} - G_i
Deleted edges  = G_i - G_{i+1}
"""

from __future__ import annotations

import argparse
from pathlib import Path


Edge = tuple[int, int]
RawEdge = tuple[int, int, int]


def edge_key(u: int, v: int) -> Edge:
    return (u, v) if u < v else (v, u)


def read_temporal(path: Path) -> list[RawEdge]:
    rows: list[RawEdge] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line[0] in "#%":
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                u, v, ts = int(parts[0]), int(parts[1]), int(float(parts[2]))
            except ValueError:
                continue
            if u != v:
                rows.append((u, v, ts))
    rows.sort(key=lambda x: x[2])
    return rows


def edges_in_window(rows: list[RawEdge], start: int, end: int) -> set[Edge]:
    return {edge_key(u, v) for u, v, ts in rows if start <= ts < end}


def infer_updates(rows: list[RawEdge], window: int, stride: int, max_events: int) -> list[tuple[str, int, int, int]]:
    if not rows:
        return []
    t_min = rows[0][2]
    t_max = rows[-1][2]
    current = edges_in_window(rows, t_min, t_min + window)
    events: list[tuple[str, int, int, int]] = []

    for e in sorted(current):
        events.append(("+", e[0], e[1], len(events)))

    start = t_min
    while start + stride + window <= t_max + 1:
        nxt_start = start + stride
        nxt = edges_in_window(rows, nxt_start, nxt_start + window)
        deleted = current - nxt
        inserted = nxt - current
        for e in sorted(deleted):
            events.append(("-", e[0], e[1], len(events)))
            if max_events and len(events) >= max_events:
                return events
        for e in sorted(inserted):
            events.append(("+", e[0], e[1], len(events)))
            if max_events and len(events) >= max_events:
                return events
        current = nxt
        start = nxt_start
    return events


def write_stream(path: Path, events: list[tuple[str, int, int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("# temporal snapshots converted to fully dynamic stream\n")
        f.write("# format: op u v timestamp\n")
        for op, u, v, ts in events:
            f.write(f"{op} {u} {v} {ts}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--out", required=True)
    parser.add_argument("--window", type=int, required=True)
    parser.add_argument("--stride", type=int)
    parser.add_argument("--max-events", type=int, default=0)
    args = parser.parse_args()

    stride = args.stride if args.stride is not None else args.window
    rows = read_temporal(Path(args.input))
    events = infer_updates(rows, args.window, stride, args.max_events)
    write_stream(Path(args.out), events)
    ins = sum(1 for e in events if e[0] == "+")
    dele = len(events) - ins
    print(f"{args.out}: {len(events)} events ({ins} inserts, {dele} deletes)")


if __name__ == "__main__":
    main()
