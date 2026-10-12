# DynEST: Dynamic and Scalable Construction of Edge-Independent Spanning Trees

CPU research artifact accompanying the manuscript **“DynEST: Dynamic and Scalable Construction of Edge-Independent Spanning Trees”** by Akanksha Dwivedi, Subhanshu Gupta, and Dip Sankar Banerjee (ICDCN 2027 submission). DynEST maintains two edge-independent spanning trees as edges are inserted and deleted, using local maintenance and selective rebuilding on the largest biconnected component.

This repository contains the C++17/OpenMP programs, dataset preparation utilities, plotting script, and archived CSVs. The archived CSVs contain synthetic experiments and a converted CollegeMsg stream; they are **not** the complete ten-dataset evaluation in Table 1 of the manuscript. Raw datasets and the manuscript are not redistributed here.

## Layout

| Path | Purpose |
| --- | --- |
| `src/d_escent.cpp` | Standalone dynamic program (`dynest`) |
| `src/benchmark.cpp` | Dynamic/static comparison and five CSV experiment series (`dynest_benchmark`) |
| `src/static_eist_construction.cpp` | Static root-aware construction program (`static_eist`) |
| `src/test_st.cpp` | Small st-numbering diagnostic (`test_st`) |
| `scripts/prepare_snap.py` | Decompress SNAP inputs and convert Reddit TSVs |
| `scripts/generate_dynamic_synthetic.py` | Generate synthetic insertion/deletion streams |
| `scripts/static_to_dynamic.py` | Convert a static edge list to a simulated update stream |
| `scripts/temporal_to_updates.py` | Convert temporal edges to window-to-window updates |
| `scripts/plot_results.py` | Plot archived benchmark CSVs |
| `results/`, `results_parallel/` | Archived sequential and OpenMP CSV measurements |

## Build

Requires CMake 3.16+, a C++17 compiler with OpenMP support, and Python 3.9+ for the scripts. On Linux (or with a working OpenMP C++ toolchain on another platform):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Executables are placed in `build/` (`.exe` on Windows). `-DDYNEST_NATIVE=ON` enables `-march=native` on supported non-MSVC compilers and produces a host-specific binary. The Python data utilities use only the standard library; plotting needs `python -m pip install -r requirements.txt`. On Windows, `py -3` can be used in place of `python` if the default interpreter does not provide `pip`.

## Get the datasets

Table 1 of the manuscript uses these **ten** SNAP datasets. Download the named file(s) from the linked SNAP page into `datasets/raw/` with the filename shown. For example, put `CollegeMsg.txt.gz` at `datasets/raw/CollegeMsg.txt.gz`. These paths are local inputs and are ignored by Git. SNAP publishes G1–G6 and G8–G10 as compressed `SRC DST UNIXTS` text. G7 is different: download **both** tab-separated Reddit files, which have subreddit names and extra columns.

| Paper ID | Dataset and SNAP source | Download to `datasets/raw/` |
| --- | --- | --- |
| G1 | [sx-stackoverflow](https://snap.stanford.edu/data/sx-stackoverflow.html) | `sx-stackoverflow.txt.gz` |
| G2 | [sx-stackoverflow-c2a](https://snap.stanford.edu/data/sx-stackoverflow.html) | `sx-stackoverflow-c2a.txt.gz` |
| G3 | [sx-stackoverflow-a2q](https://snap.stanford.edu/data/sx-stackoverflow.html) | `sx-stackoverflow-a2q.txt.gz` |
| G4 | [wiki-talk-temporal](https://snap.stanford.edu/data/wiki-talk-temporal.html) | `wiki-talk-temporal.txt.gz` |
| G5 | [sx-superuser](https://snap.stanford.edu/data/sx-superuser.html) | `sx-superuser.txt.gz` |
| G6 | [sx-askubuntu](https://snap.stanford.edu/data/sx-askubuntu.html) | `sx-askubuntu.txt.gz` |
| G7 | [soc-RedditHyperlinks](https://snap.stanford.edu/data/soc-RedditHyperlinks.html) | `soc-redditHyperlinks-body.tsv` **and** `soc-redditHyperlinks-title.tsv` |
| G8 | [sx-mathoverflow](https://snap.stanford.edu/data/sx-mathoverflow.html) | `sx-mathoverflow.txt.gz` |
| G9 | [email-Eu-core-temporal](https://snap.stanford.edu/data/email-Eu-core-temporal.html) | `email-Eu-core-temporal.txt.gz` |
| G10 | [CollegeMsg](https://snap.stanford.edu/data/CollegeMsg.html) | `CollegeMsg.txt.gz` |

Prepare one dataset (recommended to start with the small G10) or all datasets after downloading them:

```sh
python scripts/prepare_snap.py CollegeMsg
python scripts/prepare_snap.py --all
```

The utility writes uncompressed, loader-compatible files to `datasets/temporal/<name>.txt`. `--all` requires every file in the table. For Reddit it combines the body and title files, assigns integer IDs to subreddit names, parses timestamps as UTC, and discards attributes. The original TSVs remain in `datasets/raw/`. Preparing G1 requires substantial disk space because it has over 63 million temporal rows. To change the directories, see `python scripts/prepare_snap.py --help`.

The paper describes timestamp sorting, removing self-loops, retaining the earliest appearance of duplicate edges, and remapping vertices before forming its initial graph and workloads. The C++ programs sort and remap input and, for insertion-only temporal input, retain the first occurrence of each *undirected* pair. They use the remaining events for their own update workload. The raw three-column SNAP files represent interactions, **not explicit deletions**. The archived `CollegeMsg_window604800` CSVs are labeled as a seven-day window conversion; to construct such a stream:

```sh
python scripts/temporal_to_updates.py datasets/temporal/CollegeMsg.txt \
  --out datasets/temporal/CollegeMsg_window604800.txt --window 604800
```

The window converter takes edge-set differences between consecutive windows (`--stride` defaults to `--window`). This is a separate preprocessing choice from the paper's stated earliest-occurrence protocol. The exact settings used to produce the archived CollegeMsg CSVs are not recorded. The manuscript does not specify enough preprocessing detail for the Reddit body/title combination to guarantee bit-for-bit reconstruction of its G7 experiment.

## Run

The dynamic and benchmark programs accept whitespace-separated input with either `u v timestamp` (insertion-only temporal data) or explicit `+ u v timestamp` / `- u v timestamp` events. `I` and `D` are also accepted. Lines beginning with `#` or `%` are comments; self-loops are skipped and vertex IDs are remapped. Files must be **decompressed** before passing them to the C++ programs. The static program instead expects a plain two-column `u v` edge list; it writes `results/results_sequential.csv` and is not the entry point for temporal input.

```sh
./build/dynest datasets/temporal/CollegeMsg.txt --validate
./build/dynest_benchmark datasets/temporal/CollegeMsg.txt \
  --out-dir reproduced_results/CollegeMsg
```

On Windows, use `./build/dynest.exe` and `./build/dynest_benchmark.exe`. Run commands from the repository root. `dynest_benchmark` creates the requested output directory and writes `exp1_speedup_*` through `exp5_triggers_*` CSVs. It currently initializes from the first 50% of events, caps its update pool at 50,000 events, uses a main batch size of 100, and samples at most 300 static checkpoints. Thus a default run on a large dataset is an artifact demonstration, **not a full-stream reproduction of all paper figures**. These settings are constants in `src/benchmark.cpp`.

Synthetic streams corresponding to the archived `er`, `ba`, `grid`, and `regular` names can be generated with:

```sh
python scripts/generate_dynamic_synthetic.py --out-dir datasets/synthetic \
  --sizes 1000,2000,3000 --seed 42
./build/dynest_benchmark datasets/synthetic/er_random_n1000.txt \
  --out-dir reproduced_results/er_random_n1000
```

The synthetic generator exposes graph model, deletion mode, degree, update count, deletion rate, and seed as options. `static_to_dynamic.py` can also simulate updates from a plain static two-column edge list. Neither simulated stream is a substitute for the paper's real temporal datasets.

For timing comparisons, use a Release build and record compiler/version, CPU, `OMP_NUM_THREADS`, thread affinity, and input preparation. For example on Linux:

```sh
export OMP_NUM_THREADS=16
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

## Plot results

```sh
python -m pip install -r requirements.txt
python scripts/plot_results.py
```

The plotting script reads `results/` by default and writes PNG/PDF plots and `summary_table.csv` there. Use `--results-dir results_parallel` to plot the parallel archive or `--results-dir reproduced_results/CollegeMsg` for a new run. Archived timings are measurements, not expected test outputs, and may differ across hardware. See [results/README.md](results/README.md) for the CSV series.

## Citation

If you use this artifact, cite the DynEST manuscript by Akanksha Dwivedi, Subhanshu Gupta, and Dip Sankar Banerjee. Formal publication metadata can be added when available. Cite [SNAP](https://snap.stanford.edu/data/) and the source dataset(s) separately as requested by their maintainers.
