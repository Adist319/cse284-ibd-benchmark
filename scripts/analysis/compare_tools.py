import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


def load_plink_genome(path):
    # df = pd.read_csv(path, delim_whitespace=True)  # deprecated in pandas 3
    df = pd.read_csv(path, sep=r"\s+")
    df['pair'] = df.apply(
        lambda r: tuple(sorted([str(r['IID1']), str(r['IID2'])])), axis=1
    )
    return df


def load_germline_matches(path):
    cols = [
        'FID1', 'IID1', 'FID2', 'IID2', 'CHR',
        'START_BP', 'END_BP', 'START_SNP', 'END_SNP',
        'NUM_SNPS', 'GEN_LENGTH', 'UNITS', 'MISMATCHES',
        'HOM1', 'HOM2',
    ]
    df = pd.read_csv(path, sep=r"\s+", header=None, names=cols)

    # germline -haploid appends .0/.1 to IIDs, need to strip
    df['IID1'] = df['IID1'].astype(str).str.replace(r'\.\d+$', '', regex=True)
    df['IID2'] = df['IID2'].astype(str).str.replace(r'\.\d+$', '', regex=True)

    # canonicalize pair order so (A,B) == (B,A)
    mask = df['IID1'] > df['IID2']
    df.loc[mask, ['IID1', 'IID2']] = df.loc[mask, ['IID2', 'IID1']].values

    agg = (
        df.groupby(['IID1', 'IID2'])
        .agg(
            total_ibd_cm=('GEN_LENGTH', 'sum'),
            num_segments=('GEN_LENGTH', 'count'),
            max_segment_cm=('GEN_LENGTH', 'max'),
            total_snps=('NUM_SNPS', 'sum'),
        )
        .reset_index()
    )
    agg['pair'] = agg.apply(
        lambda r: tuple(sorted([r['IID1'], r['IID2']])), axis=1
    )
    return agg


def load_known_rels(path):
    df = pd.read_csv(path, sep='\t')
    df.columns = [c.lower().strip() for c in df.columns]

    # column names vary between files, just try both
    id1 = 'iid1' if 'iid1' in df.columns else 'sample1'
    id2 = 'iid2' if 'iid2' in df.columns else 'sample2'
    rel = 'relationship' if 'relationship' in df.columns else 'relationship_type'
    df = df.rename(columns={id1: 'sample1', id2: 'sample2', rel: 'relationship_type'})

    df["pair"] = df.apply(
        lambda r: tuple(sorted([str(r["sample1"]), str(r["sample2"])])), axis=1
    )
    return df


def classify_plink_rel(row):
    pi_hat = row['PI_HAT']
    z0 = row.get('Z0', 0)
    z1 = row.get('Z1', 0)
    z2 = row.get('Z2', 0)

    if pi_hat > 0.9:
        return 'duplicate'
    elif pi_hat > 0.4:
        # parent-child has z0 ~ 0, siblings have z0 ~ 0.25
        if z0 < 0.1:
            return 'parent-child'
        else:
            return 'full-sibling'
    elif pi_hat > 0.17:
        return "second-degree"
    elif pi_hat > 0.08:
        return "third-degree"
    else:
        return "unrelated"


# TODO add confusion matrix output
def classify_germline_rel(row, chr22_cm=55.0):
    total = row['total_ibd_cm']
    n_seg = row['num_segments']
    max_seg = row['max_segment_cm']
    prop = total / chr22_cm  # 55 cM for chr22

    if prop > 0.9:
        return 'duplicate'
    elif prop > 0.40:
        # parent-child = few long segments, sibs = many shorter
        if n_seg <= 3 and max_seg > 15:
            return 'parent-child'
        else:
            return 'full-sibling'
    elif prop > 0.15:
        return 'second-degree'
    elif prop > 0.07:
        return 'third-degree'
    else:
        return 'unrelated'


def plot_pihat_vs_germline(merged, output_path):
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    scatter = ax.scatter(
        merged['PI_HAT'],
        merged['total_ibd_cm'],
        c=merged['known_rel_code'],
        cmap='Set1',
        alpha=0.6, s=20,
        edgecolors='none',
    )
    ax.set_xlabel('PLINK PI_HAT', fontsize=12)
    ax.set_ylabel('GERMLINE Total IBD (cM)', fontsize=12)
    ax.set_title('PLINK PI_HAT vs GERMLINE IBD Sharing (chr22)', fontsize=13)

    rel_types = merged['known_relationship'].unique()
    colors = plt.cm.Set1(np.linspace(0, 1, len(rel_types)))
    handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=8, label=r)
        for r, c in zip(rel_types, colors)
    ]
    ax.legend(handles=handles, title='Known Relationship', loc='upper left')
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: %s" % output_path)

def plot_z0_z1(plink_df, output_path):
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    rels = plink_df['plink_class'].unique()
    palette = sns.color_palette('Set2', len(rels))
    for i, rel in enumerate(sorted(rels)):
        sub = plink_df[plink_df['plink_class'] == rel]
        ax.scatter(sub['Z0'], sub['Z1'], label=rel, alpha=0.6, s=20, color=palette[i])

    ax.set_xlabel('Z0 (P(IBD=0))', fontsize=12)
    ax.set_ylabel('Z1 (P(IBD=1))', fontsize=12)
    ax.set_title("PLINK Z0 vs Z1 - Relationship Classification", fontsize=13)
    ax.legend(title='Classified As', fontsize=9)
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_param_sensitivity(bench_path, output_path, tool_name):
    df = pd.read_csv(bench_path, sep='\t')
    if df.empty:
        print("No benchmark data for " + tool_name)
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for cohort in df['cohort'].unique():
        sub = df[df['cohort'] == cohort]
        axes[0].plot(sub.iloc[:, 2], sub['wall_time_sec'], 'o-', label=cohort)
        axes[1].plot(sub.iloc[:, 2], sub['peak_memory_kb'] / 1024, 'o-', label=cohort)

    axes[0].set_xlabel('Parameter Value')
    axes[0].set_ylabel('Wall Time (seconds)')
    axes[0].set_title('{} - Runtime vs Parameters'.format(tool_name))
    axes[0].legend()
    axes[1].set_xlabel('Parameter Value')
    axes[1].set_ylabel('Peak Memory (MB)')
    axes[1].set_title(f'{tool_name} - Memory vs Parameters')
    axes[1].legend()

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_network(classified_df, output_path, title, max_pairs=200):
    try:
        import networkx as nx
    except ImportError:
        print("networkx not available, skipping network plot")
        return

    related = classified_df[classified_df['relationship'] != 'unrelated']
    if len(related) > max_pairs:
        related = related.head(max_pairs)
    if related.empty:
        print(f"No related pairs to plot for {title}")
        return

    G = nx.Graph()
    cmap = {
        'parent-child': '#e41a1c',
        'full-sibling': '#377eb8',
        'second-degree': '#4daf4a',
        'third-degree': '#984ea3',
        'duplicate': '#ff7f00',
    }

    for _, row in related.iterrows():
        rel = row['relationship']
        G.add_edge(row['IID1'], row['IID2'], relationship=rel, color=cmap.get(rel, 'gray'))

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    pos = nx.spring_layout(G, seed=42, k=2)
    edge_colors = [G[u][v]['color'] for u, v in G.edges()]

    nx.draw_networkx(
        G, pos, ax=ax,
        node_size=100, node_color='lightblue',
        edge_color=edge_colors, width=1.5,
        font_size=6, alpha=0.8,
    )

    handles = [
        Line2D([0], [0], color=c, linewidth=2, label=r)
        for r, c in cmap.items()
        if r in related['relationship'].values
    ]
    ax.legend(handles=handles, title='Relationship', loc='upper left')
    ax.set_title(title, fontsize=14)
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def compute_accuracy(classified_df, known_df, label):
    # merge on pair key to match predicted vs known
    merged = classified_df.merge(known_df[['pair', 'relationship_type']], on='pair', how='left')
    merged['known'] = merged['relationship_type'].fillna('unrelated')
    # print(merged.head())

    res = []
    for rel_type in merged['relationship'].unique():
        pred_pos = merged['relationship'] == rel_type
        actual_pos = merged['known'] == rel_type

        tp = (pred_pos & actual_pos).sum()
        fp = (pred_pos & ~actual_pos).sum()
        fn = (~pred_pos & actual_pos).sum()

        prec = tp / (tp + fp) if (tp + fp) > 0 else 0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0

        res.append({
            'tool': label,
            'relationship': rel_type,
            'true_positives': tp,
            'false_positives': fp,
            'false_negatives': fn,
            'precision': round(prec, 3),
            'recall': round(rec, 3),
            'f1': round(f1, 3),
        })

    return (pd.DataFrame(res))


def main():
    parser = argparse.ArgumentParser(description='Compare Plink and GERMLINE IBD results')
    parser.add_argument('--cohort', required=True, help='cohort name')
    parser.add_argument('--plink-genome', required=True, help='PLINK .genome file')
    parser.add_argument('--germline-match', required=True, help='GERMLINE .match file')
    parser.add_argument('--known-rels', required=True, help='known_relationships.tsv')
    parser.add_argument('--output-dir', default=str(PROJECT_ROOT / 'results' / 'comparison'))
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir = PROJECT_ROOT / 'results' / 'figures'
    figdir.mkdir(parents=True, exist_ok=True)

    cohort = args.cohort
    print(f"\nComparing PLINK vs GERMLINE for cohort: {cohort}\n")

    print('Loading PLINK results...')
    plink = load_plink_genome(args.plink_genome)
    plink['plink_class'] = plink.apply(classify_plink_rel, axis=1)
    print(f"{len(plink)} pairs loaded")

    print("Loading GERMLINE results...")
    germline = load_germline_matches(args.germline_match)
    germline['germline_class'] = germline.apply(classify_germline_rel, axis=1)
    print("%d pairs with shared segments" % len(germline))

    print('Loading known relationships...')
    known = load_known_rels(args.known_rels)
    print(f"{len(known)} known relationships")

    # merge everything together
    # merged = plink.merge(germline, on=['IID1','IID2'])  # doesn't work, need pair keys
    merged = plink.merge(germline, on='pair', how='outer', suffixes=('_plink', '_germline'))
    merged = merged.merge(known[['pair', 'relationship_type']], on='pair', how='left')
    merged['known_relationship'] = merged['relationship_type'].fillna('unrelated')
    rel_codes = {r: i for i, r in enumerate(merged['known_relationship'].unique())}
    merged['known_rel_code'] = merged['known_relationship'].map(rel_codes)
    merged['PI_HAT'] = merged['PI_HAT'].fillna(0)
    merged['total_ibd_cm'] = merged['total_ibd_cm'].fillna(0)

    print(f"\nConcordance summary:")
    both = merged[(merged['PI_HAT'] > 0.05) | (merged['total_ibd_cm'] > 3)]
    if len(both) > 0:
        corr = both[['PI_HAT', 'total_ibd_cm']].corr().iloc[0, 1]
        print(f"Pearson correlation (PI_HAT vs total IBD cM): {corr:.3f}")
    print(f"Pairs found by PLINK (PI_HAT > 0.05): {(plink['PI_HAT'] > 0.05).sum()}")
    print(f"Pairs found by GERMLINE: {len(germline)}")

    agree = merged.dropna(subset=['plink_class', 'germline_class'])
    agreement = (agree['plink_class'] == agree['germline_class']).mean()
    print(f"Classification agreement: {agreement:.1%}")

    out_path = outdir / f"{cohort}_merged_comparison.tsv"
    merged.to_csv(out_path, sep='\t', index=False)
    print(f"\nSaved merged comparison: {out_path}")

    print('\nGenerating plots...')
    plot_pihat_vs_germline(
        merged[merged['PI_HAT'] > 0.01],
        figdir / f'{cohort}_pihat_vs_germline_ibd.png',
    )
    plot_z0_z1(plink, figdir / f'{cohort}_z0_z1_scatter.png')

    # network graphs
    plink_classified = plink[['IID1', 'IID2', 'pair', 'plink_class']].rename(
        columns={'plink_class': 'relationship'}
    )
    germline_classified = germline[['IID1', 'IID2', 'pair', 'germline_class']].rename(
        columns={'germline_class': 'relationship'}
    )
    plot_network(
        plink_classified, figdir / f'{cohort}_network_plink.png',
        f'PLINK Family Network - {cohort}',
    )
    plot_network(
        germline_classified, figdir / f'{cohort}_network_germline.png',
        f'GERMLINE Family Network - {cohort}',
    )

    # TODO should probably do per-cohort accuracy separately instead of lumping
    print('\nComputing accuracy...')
    plink_acc = compute_accuracy(plink_classified, known, 'PLINK')
    germline_acc = compute_accuracy(germline_classified, known, 'GERMLINE')
    accuracy = pd.concat([plink_acc, germline_acc], ignore_index=True)
    acc_path = outdir / f'{cohort}_accuracy.tsv'
    accuracy.to_csv(acc_path, sep='\t', index=False)
    print(f"Saved: {acc_path}")
    print(accuracy.to_string(index=False))

    plink_bench = PROJECT_ROOT / 'results' / 'plink' / 'benchmarks.tsv'
    germline_bench = PROJECT_ROOT / 'results' / 'germline' / 'benchmarks.tsv'
    if plink_bench.exists():
        plot_param_sensitivity(plink_bench, figdir / 'plink_benchmarks.png', 'PLINK')
    if germline_bench.exists():
        plot_param_sensitivity(germline_bench, figdir / 'germline_benchmarks.png', 'GERMLINE')

    print(f"\nComparison complete for {cohort}")


if __name__ == '__main__':
    main()
