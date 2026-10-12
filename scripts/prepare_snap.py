#!/usr/bin/env python3
"""Prepare the ten SNAP inputs listed in Table 1 of the DynEST manuscript.

Nine datasets are gzip-compressed three-column temporal edge lists. Reddit is
provided as two TSVs with string vertex IDs and additional attributes.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
from datetime import datetime, timezone
from pathlib import Path


SNAP_NAMES = (
    "sx-stackoverflow",
    "sx-stackoverflow-c2a",
    "sx-stackoverflow-a2q",
    "wiki-talk-temporal",
    "sx-superuser",
    "sx-askubuntu",
    "soc-RedditHyperlinks",
    "sx-mathoverflow",
    "email-Eu-core-temporal",
    "CollegeMsg",
)
REDDIT_FILES = (
    "soc-redditHyperlinks-body.tsv",
    "soc-redditHyperlinks-title.tsv",
)


def input_names(name: str) -> tuple[str, ...]:
    if name == "soc-RedditHyperlinks":
        return REDDIT_FILES
    return (f"{name}.txt.gz",)


def prepare_gzip(source: Path, destination: Path) -> None:
    with gzip.open(source, "rb") as compressed, destination.open("wb") as output:
        shutil.copyfileobj(compressed, output, length=1024 * 1024)


def prepare_reddit(raw_dir: Path, destination: Path) -> None:
    vertex_ids: dict[str, int] = {}
    events = 0
    with destination.open("w", encoding="utf-8", newline="\n") as output:
        output.write("# Reddit body and title hyperlinks; SRC DST UNIXTS (UTC)\n")
        for name in REDDIT_FILES:
            with (raw_dir / name).open("r", encoding="utf-8-sig", newline="") as source:
                for line_number, row in enumerate(csv.reader(source, delimiter="\t"), 1):
                    if not row or row[0] == "SOURCE_SUBREDDIT":
                        continue
                    if len(row) < 4:
                        raise ValueError(f"{name}:{line_number}: expected at least four TSV columns")
                    src, dst, _, raw_time = row[:4]
                    if src == dst:
                        continue
                    try:
                        timestamp = int(
                            datetime.strptime(raw_time, "%Y-%m-%d %H:%M:%S")
                            .replace(tzinfo=timezone.utc)
                            .timestamp()
                        )
                    except ValueError as exc:
                        raise ValueError(f"{name}:{line_number}: invalid timestamp {raw_time!r}") from exc
                    u = vertex_ids.setdefault(src, len(vertex_ids))
                    v = vertex_ids.setdefault(dst, len(vertex_ids))
                    output.write(f"{u} {v} {timestamp}\n")
                    events += 1
    print(f"  {events} Reddit temporal rows, {len(vertex_ids)} vertices")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", nargs="?", choices=SNAP_NAMES, help="SNAP dataset name")
    parser.add_argument("--all", action="store_true", help="prepare all ten paper datasets")
    parser.add_argument("--raw-dir", type=Path, default=Path("datasets/raw"))
    parser.add_argument("--out-dir", type=Path, default=Path("datasets/temporal"))
    args = parser.parse_args()
    if args.all == bool(args.dataset):
        parser.error("specify one dataset name or --all")

    names = SNAP_NAMES if args.all else (args.dataset,)
    missing = [args.raw_dir / filename for name in names for filename in input_names(name)
               if not (args.raw_dir / filename).is_file()]
    if missing:
        parser.error("missing input file(s): " + ", ".join(map(str, missing)))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name in names:
        destination = args.out_dir / f"{name}.txt"
        if name == "soc-RedditHyperlinks":
            prepare_reddit(args.raw_dir, destination)
        else:
            prepare_gzip(args.raw_dir / input_names(name)[0], destination)
        print(f"Prepared {destination}")


if __name__ == "__main__":
    main()
