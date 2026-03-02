import argparse
import sys
from pathlib import Path

import pandas as pd

MATCH_COLS = [
    'fid1', 'iid1', 'fid2', 'iid2',
    'chromosome', 'seg_start_bp', 'seg_end_bp',
    'seg_start_snp', 'seg_end_snp', 'total_snps',
    'genetic_length', 'units', 'mismatches',
    'is_hom1', 'is_hom2',
]

CHR22_CM = 55.0

# thresholds calibrated for chr22 only, won't generalize
# originally had full_sibling min at 25 but that missed too many
THRESHOLDS = {
    'parent_child': {'min_ibd_cm': 33.0, 'max_segments': 4},
    'full_sibling': {'min_ibd_cm': 20.0, 'max_segments': 10},
    'second_degree': {'min_ibd_cm': 10.0},
    'third_degree': {'min_ibd_cm': 5.0},
}


def read_match_file(path):
    df = pd.read_csv(
        path, sep=r"\s+", header=None, names=MATCH_COLS,
        dtype={'fid1': str, 'iid1': str, 'fid2': str, 'iid2': str, 'chromosome': str},
    )
    return df

def aggregate_pairs(df):
    if df.empty:
        return pd.DataFrame(columns=[
            'fid1', 'iid1', 'fid2', 'iid2',
            'total_ibd_cm', 'num_segments', 'max_segment_cm',
            'proportion_shared',
        ])

    # strip .0/.1 haplotype suffixes from germline -haploid mode
    df = df.copy()
    df['iid1'] = df['iid1'].str.replace(r'\.\d+$', '', regex=True)
    df['iid2'] = df['iid2'].str.replace(r'\.\d+$', '', regex=True)
    df['fid1'] = df['fid1'].str.replace(r'\.\d+$', '', regex=True)
    df['fid2'] = df['fid2'].str.replace(r'\.\d+$', '', regex=True)

    # canonicalize so pair order is consistent
    mask = df['iid1'] > df['iid2']
    df.loc[mask, ['fid1', 'iid1', 'fid2', 'iid2']] = df.loc[mask, ['fid2', 'iid2', 'fid1', 'iid1']].values

    grouped = df.groupby(['fid1', 'iid1', 'fid2', 'iid2'], as_index=False)
    pairs = grouped.agg(
        total_ibd_cm=('genetic_length', 'sum'),
        num_segments=('genetic_length', 'count'),
        max_segment_cm=('genetic_length', 'max'),
    )
    pairs['proportion_shared'] = pairs['total_ibd_cm'] / 55.0  # chr22
    return pairs


def classify_rel(row):
    total = row['total_ibd_cm']
    n_seg = row['num_segments']

    if total >= THRESHOLDS['parent_child']['min_ibd_cm']:
        # parent-child = few long segments vs siblings = more from recombination
        if n_seg <= THRESHOLDS['parent_child']['max_segments']:
            return 'parent-child'
        return 'full-sibling'

    if total >= THRESHOLDS['full_sibling']['min_ibd_cm']:
        return 'full-sibling'
    if total >= THRESHOLDS['second_degree']['min_ibd_cm']:
        return 'second-degree'
    if total >= THRESHOLDS['third_degree']['min_ibd_cm']:
        return 'third-degree'
    return 'unrelated'


def load_known_rels(path):
    df = pd.read_csv(path, sep='\t', dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]
    # handle different column naming conventions
    renames = {}
    if 'sample1' in df.columns:
        renames['sample1'] = 'iid1'
    if 'sample2' in df.columns:
        renames['sample2'] = 'iid2'
    if 'relationship_type' in df.columns:
        renames['relationship_type'] = 'relationship'
    if renames:
        df = df.rename(columns=renames)
    return df


def evaluate_accuracy(classified, known):
    def make_pair_key(df, id1, id2):
        left = df[id1].astype(str)
        right = df[id2].astype(str)
        a = left.where(left <= right, right)
        b = right.where(left <= right, left)
        return a + '_' + b

    classified = classified.copy()
    classified['pair_key'] = make_pair_key(classified, 'iid1', 'iid2')
    known = known.copy()
    known['pair_key'] = make_pair_key(known, 'iid1', 'iid2')

    merged = classified.merge(
        known[['pair_key', 'relationship']],
        on='pair_key', how='inner',
        suffixes=('', '_known'),
    )
    merged = merged.rename(columns={'relationship': 'true_relationship'})
    merged['correct'] = merged['predicted_relationship'] == merged['true_relationship']
    # print(f"DEBUG evaluate_accuracy: {len(merged)} merged pairs")
    return merged


def process_cohort(match_path, cohort, known_path, output_dir):
    print(f"\nProcessing cohort: {cohort}")
    print(f"Match file: {match_path}")

    df = read_match_file(match_path)
    print("Total IBD segments: %d" % len(df))

    pairs = aggregate_pairs(df)
    print(f"Unique pairs with IBD: {len(pairs)}")

    pairs['predicted_relationship'] = pairs.apply(classify_rel, axis=1)

    pairs_out = output_dir / f'{cohort}_pairs_summary.tsv'
    pairs.to_csv(pairs_out, sep='\t', index=False)
    print(f"Pairs summary: {pairs_out}")

    classified_out = output_dir / f'{cohort}_classified.tsv'
    out_cols = [
        'fid1', 'iid1', 'fid2', 'iid2',
        'total_ibd_cm', 'num_segments', 'predicted_relationship',
    ]
    pairs[out_cols].to_csv(classified_out, sep='\t', index=False)
    print(f"Classified output: {classified_out}")

    dist = pairs['predicted_relationship'].value_counts()
    print(f"Classification distribution:")
    for rel, cnt in dist.items():
        print(f"  {rel}: {cnt}")

    res = {
        'cohort': cohort,
        'total_pairs': len(pairs),
        'total_segments': len(df),
    }

    if known_path is not None and known_path.exists():
        known = load_known_rels(known_path)
        merged = evaluate_accuracy(pairs, known)

        if len(merged) > 0:
            acc = merged['correct'].mean()
            n_eval = len(merged)
            print(f"Accuracy vs ground truth: {acc:.3f} ({n_eval} pairs evaluated)")

            res['pairs_evaluated'] = n_eval
            res['correct'] = int(merged['correct'].sum())
            res['accuracy'] = round(acc, 4)

            for rel in merged['true_relationship'].unique():
                sub = merged[merged['true_relationship'] == rel]
                rel_acc = sub['correct'].mean()
                rel_n = len(sub)
                print("  {}: {:.3f} ({} pairs)".format(rel, rel_acc, rel_n))
                res[f'accuracy_{rel}'] = round(rel_acc, 4)
                res[f'n_{rel}'] = rel_n
        else:
            print("No overlapping pairs with ground truth.")
            res['pairs_evaluated'] = 0
            res['accuracy'] = None
    else:
        print("No ground truth file provided, skipping accuracy evaluation.")
        res['pairs_evaluated'] = 0
        res['accuracy'] = None

    return res


def main():
    parser = argparse.ArgumentParser(
        description='Classify relationships from GERMLINE IBD output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="NOTE: thresholds are for chr22 only (~55 cM). don't use for whole-genome.",
    )
    parser.add_argument('--match-files', nargs='+', required=True,
                        help='cohort:path format, e.g. trios:results/trios.match')
    parser.add_argument('--known', type=Path, default=None, help='ground truth tsv')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    acc_results = []

    for match_spec in args.match_files:
        if ':' in match_spec:
            cohort, mpath = match_spec.split(':', 1)
        else:
            mpath = match_spec
            cohort = Path(match_spec).stem.split('_')[0]

        match_path = Path(mpath)
        if not match_path.exists():
            print(f"WARNING: Match file not found: {match_path}, skipping.", file=sys.stderr)
            continue

        res = process_cohort(match_path, cohort, args.known, args.output_dir)
        acc_results.append(res)

    if acc_results:
        acc_df = pd.DataFrame(acc_results)
        acc_out = args.output_dir / 'accuracy_summary.tsv'
        acc_df.to_csv(acc_out, sep='\t', index=False)
        print(f"\nAccuracy summary written to: {acc_out}")

        print('\nOverall summary:')
        for row in acc_results:
            acc_str = f"{row['accuracy']:.4f}" if row.get('accuracy') is not None else 'N/A'
            print(
                f"  {row['cohort']}: "
                f"{row['total_pairs']} pairs, "
                f"{row['total_segments']} segments, "
                f"accuracy={acc_str}"
            )

        # thresholds only work for chr22, not genome-wide
        print(
            "\nNOTE: Thresholds are approximate for chr22 only. "
            "Real analysis should use genome-wide IBD data."
        )


if __name__ == '__main__':
    main()
