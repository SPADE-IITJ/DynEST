"""
D-ESCENT Benchmark Visualization

Reads benchmark CSVs and produces plots in the same directory.
Handles missing/partial data gracefully (timeout safeguard).

Usage: python scripts/plot_results.py [--results-dir DIR]
"""

import argparse
import os
import glob
import warnings

warnings.filterwarnings('ignore')

try:
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
except ImportError as e:
    print(f"Missing dependency: {e}")
    print("Install with: pip install pandas matplotlib")
    exit(1)

# ─── Style ───────────────────────────────────────────────────────
plt.rcParams.update({
    'figure.figsize': (8, 5),
    'figure.dpi': 150,
    'font.size': 11,
    'font.family': 'serif',
    'axes.grid': True,
    'grid.alpha': 0.3,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'legend.framealpha': 0.9,
    'legend.edgecolor': '0.8',
})

RESULTS_DIR = 'results'
DATASET_ORDER = [
    'CollegeMsg',
    'email-Eu-core-temporal',
    'sx-mathoverflow',
    'sx-askubuntu',
    'sx-superuser',
    'wiki-talk-temporal',
    'sx-stackoverflow',
]
SHORT_NAMES = {
    'CollegeMsg': 'College',
    'email-Eu-core-temporal': 'Email-EU',
    'sx-mathoverflow': 'MathOF',
    'sx-askubuntu': 'AskUbuntu',
    'sx-superuser': 'SuperUser',
    'wiki-talk-temporal': 'WikiTalk',
    'sx-stackoverflow': 'StackOF',
}

COLORS = plt.cm.tab10.colors


def safe_read(path):
    """Read CSV, return None if missing or empty."""
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path)
        if df.empty:
            return None
        return df
    except Exception:
        return None


def discover_datasets():
    """Find which datasets have exp1 results."""
    found = []
    for name in DATASET_ORDER:
        p = os.path.join(RESULTS_DIR, f'exp1_speedup_{name}.csv')
        if os.path.exists(p):
            found.append(name)
    # Also pick up any not in the predefined order
    for f in glob.glob(os.path.join(RESULTS_DIR, 'exp1_speedup_*.csv')):
        stem = os.path.basename(f).replace('exp1_speedup_', '').replace('.csv', '')
        if stem not in found:
            found.append(stem)
    return found


# ─── Exp-I: Speedup vs Batch Size ────────────────────────────────
def plot_exp1(datasets):
    fig, ax = plt.subplots()
    has_data = False

    for i, ds in enumerate(datasets):
        df = safe_read(os.path.join(RESULTS_DIR, f'exp1_speedup_{ds}.csv'))
        if df is None:
            continue

        # Mean speedup per batch size
        batch_col = 'batch_size' if 'batch_size' in df.columns else 'B'
        grp = df.groupby(batch_col)['speedup'].mean().reset_index()
        if grp.empty:
            continue

        label = SHORT_NAMES.get(ds, ds)
        ax.plot(grp[batch_col], grp['speedup'], 'o-',
                color=COLORS[i % len(COLORS)], label=label,
                markersize=6, linewidth=2)
        has_data = True

    if not has_data:
        plt.close(); return

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Batch Size (B)')
    ax.set_ylabel('Mean Speedup (Static / Dynamic)')
    ax.set_title('Exp-I: Maintenance Speedup vs Batch Size')
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.legend(fontsize=9)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, 'exp1_speedup.png'), dpi=200)
    fig.savefig(os.path.join(RESULTS_DIR, 'exp1_speedup.pdf'))
    plt.close()
    print("  [Exp-I] exp1_speedup.png / .pdf")


# ─── Exp-I supplemental: per-round speedup timeseries ────────────
def plot_exp1_timeseries(datasets):
    for ds in datasets:
        df = safe_read(os.path.join(RESULTS_DIR, f'exp1_speedup_{ds}.csv'))
        if df is None:
            continue

        if 'edges_processed' not in df.columns:
            continue
        batch_col = 'batch_size' if 'batch_size' in df.columns else 'B'
        batch_sizes = sorted(df[batch_col].unique())
        fig, ax = plt.subplots()
        for j, B in enumerate(batch_sizes):
            sub = df[df[batch_col] == B]
            ax.plot(sub['edges_processed'], sub['speedup'],
                    color=COLORS[j % len(COLORS)], label=f'B={B}',
                    linewidth=1.5, alpha=0.8)

        ax.set_xlabel('Edges Processed')
        ax.set_ylabel('Speedup (Static / Dynamic)')
        label = SHORT_NAMES.get(ds, ds)
        ax.set_title(f'Exp-I: Speedup Over Time — {label}')
        ax.axhline(y=1, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS_DIR, f'exp1_timeseries_{ds}.png'), dpi=200)
        plt.close()
    print("  [Exp-I] per-dataset timeseries saved")


# ─── Exp-II: Fast Path Rate Breakdown ────────────────────────────
def plot_exp2(datasets):
    for ds in datasets:
        df = safe_read(os.path.join(RESULTS_DIR, f'exp2_fastpath_{ds}.csv'))
        if df is None:
            continue

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
        label = SHORT_NAMES.get(ds, ds)

        x = df['edges_processed'] if 'edges_processed' in df.columns else df['n_edges']
        fast = df['fast_ins'] if 'fast_ins' in df.columns else df['fast']
        slow = df['slow_ins'] if 'slow_ins' in df.columns else df['slow']
        fast_pct = df['fast_rate_pct'] if 'fast_rate_pct' in df.columns else df['fast_pct']

        # Top: stacked bar of fast vs slow
        w = (x.max() - x.min()) / len(x) * 0.8 if len(x) > 1 else 10
        ax1.bar(x, fast, width=w, label='Fast', color=COLORS[0], alpha=0.8)
        ax1.bar(x, slow, width=w, bottom=fast,
                label='Slow', color=COLORS[3], alpha=0.8)
        ax1.set_ylabel('Operations per Round')
        ax1.set_title(f'Exp-II: Fast vs Slow Path — {label} (B=100)')
        ax1.legend(fontsize=9)

        # Bottom: fast rate %
        ax2.plot(x, fast_pct, color=COLORS[2], linewidth=2)
        ax2.set_xlabel('Edges Processed')
        ax2.set_ylabel('Fast Path Rate (%)')
        ax2.set_ylim(0, 105)

        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS_DIR, f'exp2_fastpath_{ds}.png'), dpi=200)
        plt.close()
    print("  [Exp-II] fast path plots saved")


# ─── Exp-III: Tree Quality Over Time ─────────────────────────────
def plot_exp3(datasets):
    fig, ax = plt.subplots()
    has_data = False

    for i, ds in enumerate(datasets):
        df = safe_read(os.path.join(RESULTS_DIR, f'exp3_quality_{ds}.csv'))
        if df is None:
            continue
        label = SHORT_NAMES.get(ds, ds)
        x = df['edges_processed'] if 'edges_processed' in df.columns else df['n_edges']
        ratio = df['ratio_dyn_to_race'] if 'ratio_dyn_to_race' in df.columns else df['ratio']
        ax.plot(x, ratio,
                color=COLORS[i % len(COLORS)], label=label,
                linewidth=1.5, alpha=0.8)
        has_data = True

    if not has_data:
        plt.close(); return

    ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.set_xlabel('Edges Processed')
    ax.set_ylabel('Cost Ratio (Dynamic / Static RACE)')
    ax.set_title('Exp-III: Dynamic Tree Quality vs Fresh Static (K=32)')
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, 'exp3_quality.png'), dpi=200)
    fig.savefig(os.path.join(RESULTS_DIR, 'exp3_quality.pdf'))
    plt.close()

    # Also per-dataset: dynamic vs static vs arbitrary
    for ds in datasets:
        df = safe_read(os.path.join(RESULTS_DIR, f'exp3_quality_{ds}.csv'))
        if df is None:
            continue
        fig, ax = plt.subplots()
        label = SHORT_NAMES.get(ds, ds)
        x = df['edges_processed'] if 'edges_processed' in df.columns else df['n_edges']
        dyn_col = 'cost_dynamic' if 'cost_dynamic' in df.columns else 'dyn_cost'
        sta_col = 'cost_static_race' if 'cost_static_race' in df.columns else 'sta_cost'
        ax.plot(x, df[dyn_col], label='Dynamic', linewidth=2)
        ax.plot(x, df[sta_col], label='Static RACE (K=32)',
                linewidth=2, linestyle='--')
        if 'cost_static_arb' in df.columns:
            ax.plot(x, df['cost_static_arb'], label='Static Arb (K=1)',
                    linewidth=2, linestyle=':')
        ax.set_xlabel('Edges Processed')
        ax.set_ylabel('EIST Cost')
        ax.set_title(f'Exp-III: Tree Cost Comparison — {label}')
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS_DIR, f'exp3_cost_{ds}.png'), dpi=200)
        plt.close()
    print("  [Exp-III] quality plots saved")


# ─── Exp-IV: Graph-size or thread scaling ────────────────────────
def plot_exp4(datasets):
    for ds in datasets:
        df = safe_read(os.path.join(RESULTS_DIR, f'exp4_scaling_{ds}.csv'))
        if df is None:
            continue

        fig, ax = plt.subplots()
        label = SHORT_NAMES.get(ds, ds)

        if 'threads' not in df.columns:
            x = df['n_edges'] if 'n_edges' in df.columns else df['edges_processed']
            ax.plot(x, df['dyn_ms'], label='Dynamic batch', linewidth=2)
            ax.plot(x, df['sta_ms'], label='Static rebuild', linewidth=2, linestyle='--')
            ax.set_xlabel('Edges Processed')
            ax.set_ylabel('Execution Time (ms)')
            ax.set_title(f'Exp-IV: Scaling with Graph Size - {label}')
            ax.legend(fontsize=9)
            fig.tight_layout()
            fig.savefig(os.path.join(RESULTS_DIR, f'exp4_scaling_{ds}.png'), dpi=200)
            plt.close()
            continue

        for j, K in enumerate(sorted(df['K'].unique())):
            sub = df[df['K'] == K].groupby('threads')['speedup_vs_1thread'].mean().reset_index()
            ax.plot(sub['threads'], sub['speedup_vs_1thread'], 'o-',
                    color=COLORS[j % len(COLORS)], label=f'K={K}',
                    markersize=6, linewidth=2)

        # Ideal line
        thr_max = df['threads'].max()
        ax.plot([1, thr_max], [1, thr_max], 'k--', alpha=0.3, label='Ideal')

        ax.set_xlabel('Threads')
        ax.set_ylabel('Speedup vs 1 Thread')
        ax.set_title(f'Exp-IV: Parallel Scaling — {label}')
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS_DIR, f'exp4_scaling_{ds}.png'), dpi=200)
        plt.close()
    print("  [Exp-IV] scaling plots saved")


# ─── Exp-V: Rebuild Triggers ─────────────────────────────────────
def plot_exp5(datasets):
    for ds in datasets:
        df = safe_read(os.path.join(RESULTS_DIR, f'exp5_triggers_{ds}.csv'))
        if df is None:
            continue

        fig, ax = plt.subplots()
        label = SHORT_NAMES.get(ds, ds)
        x = df['edges_processed'] if 'edges_processed' in df.columns else df['n_edges']
        w = (x.max() - x.min()) / len(x) * 0.8 if len(x) > 1 else 10

        if 'cum_trig_merge' in df.columns:
            merge = df['cum_trig_merge'].diff().fillna(df['cum_trig_merge'])
            periodic = df['cum_trig_periodic'].diff().fillna(df['cum_trig_periodic'])
            pending = df['cum_trig_pending'].diff().fillna(df['cum_trig_pending'])
            ax.bar(x, merge, width=w, label='Merge (Eager)', color=COLORS[1], alpha=0.8)
            ax.bar(x, periodic, width=w, bottom=merge,
                   label='Periodic (Cost)', color=COLORS[2], alpha=0.8)
            ax.bar(x, pending, width=w, bottom=merge + periodic,
                   label='Pending Merge', color=COLORS[4], alpha=0.8)
        else:
            ax.bar(x, df['no_replacement'], width=w, label='No Replacement',
                   color=COLORS[3], alpha=0.8)
            ax.bar(x, df['merge_eager'], width=w, bottom=df['no_replacement'],
                   label='Merge (Eager)', color=COLORS[1], alpha=0.8)
            ax.bar(x, df['periodic_cost'], width=w,
                   bottom=df['no_replacement'] + df['merge_eager'],
                   label='Periodic (Cost)', color=COLORS[2], alpha=0.8)

        ax.set_xlabel('Edges Processed')
        ax.set_ylabel('Rebuilds per Round')
        ax.set_title(f'Exp-V: Rebuild Triggers — {label} (B=100)')
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS_DIR, f'exp5_triggers_{ds}.png'), dpi=200)
        plt.close()
    print("  [Exp-V] trigger plots saved")


# ─── Summary table ───────────────────────────────────────────────
def make_summary(datasets):
    rows = []
    for ds in datasets:
        df = safe_read(os.path.join(RESULTS_DIR, f'exp1_speedup_{ds}.csv'))
        if df is None:
            continue
        label = SHORT_NAMES.get(ds, ds)
        batch_col = 'batch_size' if 'batch_size' in df.columns else 'B'
        for B in sorted(df[batch_col].unique()):
            sub = df[df[batch_col] == B]
            rows.append({
                'Dataset': label,
                'B': int(B),
                'Samples': int(sub['S_samples'].iloc[0]) if 'S_samples' in sub.columns else len(sub),
                'Mean Speedup': f"{sub['speedup'].mean():.2f}",
                'Median Speedup': f"{sub['speedup'].median():.2f}",
                'Dynamic total (ms)': f"{sub['dyn_total_ms'].mean():.2f}" if 'dyn_total_ms' in sub.columns else f"{sub['t_dynamic_ms'].mean():.2f}",
                'Static total (ms)': f"{sub['sta_total_ms'].mean():.2f}" if 'sta_total_ms' in sub.columns else f"{sub['t_static_ms'].mean():.2f}",
            })

    if rows:
        summary = pd.DataFrame(rows)
        summary.to_csv(os.path.join(RESULTS_DIR, 'summary_table.csv'), index=False)
        print(f"\n  Summary table -> {RESULTS_DIR}/summary_table.csv")
        print(summary.to_string(index=False))


# ─── Main ────────────────────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot DynEST benchmark CSVs')
    parser.add_argument('--results-dir', default='results',
                        help='CSV input and figure output directory (default: results)')
    args = parser.parse_args()
    RESULTS_DIR = args.results_dir
    print("D-ESCENT Benchmark Visualization")
    print("=" * 40)

    datasets = discover_datasets()
    if not datasets:
        print(f"No result CSVs found in {RESULTS_DIR}/. Run benchmark first.")
        exit(0)

    print(f"Found data for: {', '.join(datasets)}")

    plot_exp1(datasets)
    plot_exp1_timeseries(datasets)
    plot_exp2(datasets)
    plot_exp3(datasets)
    plot_exp4(datasets)
    plot_exp5(datasets)
    make_summary(datasets)

    print(f"\nAll plots saved to {RESULTS_DIR}/")
