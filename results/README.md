# Archived benchmark results

This directory contains CSV output produced by `src/benchmark.cpp` for the
sequential/reference experiment runs. The dataset suffixes here are synthetic
`er`, `ba`, `grid`, and `regular` streams plus `CollegeMsg_window604800`.
They do not constitute the ten real-world datasets in Table 1 of the paper.
Files are grouped by experiment prefix:

- `exp1_speedup_*`: batch-size timing and speedup
- `exp2_fastpath_*`: fast-path and slow-path update counts
- `exp3_quality_*`: dynamic versus static tree-cost comparison
- `exp4_scaling_*`: runtime and largest-BiCC size over update rounds (not a thread-scaling sweep)
- `exp5_triggers_*`: cumulative reconstruction-trigger counts

These files are archived measurements rather than golden test outputs. Hardware,
compiler, OpenMP configuration, and dataset preprocessing must be reported when
using them in a publication. New runs should use `--out-dir` and write to a
separate directory such as `reproduced_results/`. The sibling
`results_parallel/` directory holds a second archive of OpenMP measurements;
its files use the same CSV schema.
