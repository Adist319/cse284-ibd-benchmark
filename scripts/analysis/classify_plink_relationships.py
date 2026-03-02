import argparse
import sys
from pathlib import Path

import pandas as pd

REL_TYPES = [
    "identical_twin",
    "parent_child",
    "full_sibling",
    "second_degree",
    "third_degree",
    "unrelated",
]


def classify_rel(row):
    pi_hat = row["PI_HAT"]
    z1 = row["Z1"]
    z2 = row["Z2"]

    # used to threshold at 0.20 but 0.17 works better empirically
    if pi_hat > 0.9:
        return "identical_twin"
    if pi_hat > 0.4 and z1 > 0.8:
        return "parent_child"
    if pi_hat > 0.35 and 0.4 < z1 < 0.7 and z2 > 0.1:
        return "full_sibling"
    if 0.17 < pi_hat < 0.35:
        return "second_degree"
    if 0.08 < pi_hat < 0.17:
        return "third_degree"
    return "unrelated"


def read_genome(path):
    data = pd.read_csv(path, sep=r"\s+")
    data.columns = data.columns.str.strip()
    return data


def read_known_rels(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()

    col_map = {c: c.lower().strip() for c in df.columns}
    df = df.rename(columns=col_map)
    id1 = "iid1" if "iid1" in df.columns else df.columns[0]
    id2 = "iid2" if "iid2" in df.columns else df.columns[1]
    rel = "relationship" if "relationship" in df.columns else df.columns[2]

    # originally had this as a dict comprehension but maybe a little too unreadable
    pairs = pd.DataFrame(
        {
            "IID1": df.apply(lambda r: min(str(r[id1]), str(r[id2])), axis=1),
            "IID2": df.apply(lambda r: max(str(r[id1]), str(r[id2])), axis=1),
            "known_relationship": df[rel],
        }
    )
    return pairs


def classify_cohort(genome_path):
    df = read_genome(genome_path)
    df["predicted_relationship"] = df.apply(classify_rel, axis=1)
    return df


def normalize_rel(rel):
    r = rel.strip().lower()
    mapping = {
        "parent-child": "parent_child",
        "parent_child": "parent_child",
        "po": "parent_child",
        "full_sibling": "full_sibling",
        "full-sibling": "full_sibling",
        "sibling": "full_sibling",
        "fs": "full_sibling",
        "identical_twin": "identical_twin",
        "identical-twin": "identical_twin",
        "duplicate": "identical_twin",
        "mz_twin": "identical_twin",
        "second_degree": "second_degree",
        "second-degree": "second_degree",
        "half_sibling": "second_degree",
        "half-sibling": "second_degree",
        "avuncular": "second_degree",
        "grandparent": "second_degree",
        "hs": "second_degree",
        "third_degree": "third_degree",
        "third-degree": "third_degree",
        "first_cousin": "third_degree",
        "first-cousin": "third_degree",
        "cousin": "third_degree",
        "unrelated": "unrelated",
    }
    return mapping.get(r, r)


def evaluate_accuracy(classified_df, known_df):
    # normalize pair ordering so (A,B) matches (B,A)
    c = classified_df.copy()
    c["IID1_norm"] = c.apply(lambda r: min(str(r["IID1"]), str(r["IID2"])), axis=1)
    c["IID2_norm"] = c.apply(lambda r: max(str(r["IID1"]), str(r["IID2"])), axis=1)

    # merge on IID
    merged = c.merge(
        known_df,
        left_on=["IID1_norm", "IID2_norm"],
        right_on=["IID1", "IID2"],
        how="inner",
        suffixes=("", "_known"),
    )
    # print(merged.head())

    if merged.empty:
        print("WARNING: No overlapping pairs found with ground truth")
        return pd.DataFrame(
            columns=["relationship", "precision", "recall", "f1", "tp", "fp", "fn"]
        )

    res = []
    all_known = known_df["known_relationship"].apply(normalize_rel)
    all_pred = merged["predicted_relationship"]
    all_true = merged["known_relationship"].apply(normalize_rel)

    for rel in REL_TYPES:
        tp = ((all_pred == rel) & (all_true == rel)).sum()
        fp = ((all_pred == rel) & (all_true != rel)).sum()
        fn = ((all_pred != rel) & (all_true == rel)).sum()

        prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0

        res.append(
            {
                "relationship": rel,
                "precision": round(prec, 4),
                "recall": round(rec, 4),
                "f1": round(f1, 4),
                "tp": int(tp),
                "fp": int(fp),
                "fn": int(fn),
            }
        )

    return pd.DataFrame(res)


def main():
    parser = argparse.ArgumentParser(
        description="Classify relationships from PLINK .genome IBD output"
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/plink"),
        help="Directory containing .genome files (default: results/plink)",
    )
    parser.add_argument(
        "--known-relationships",
        type=Path,
        default=Path("data/processed/known_relationships.tsv"),
        help="Path to known relationships TSV (default: data/processed/known_relationships.tsv)",
    )
    parser.add_argument(
        "--cohorts",
        nargs="+",
        default=["admixed", "homogeneous", "trios"],
    )
    parser.add_argument(
        "--suffix",
        default="default",
        help="Genome file suffix, e.g. 'default' for {cohort}_default.genome",
    )
    args = parser.parse_args()

    results_dir = args.results_dir
    if not results_dir.exists():
        print("ERROR: Results directory not found: %s" % results_dir, file=sys.stderr)
        sys.exit(1)

    known_df = None
    if args.known_relationships.exists():
        known_df = read_known_rels(args.known_relationships)
        print(f"Loaded {len(known_df)} known relationship pairs")
    else:
        print(
            f"WARNING: Known relationships file not found: {args.known_relationships}"
        )
        print("Skipping accuracy evaluation")

    all_acc = []

    for cohort in args.cohorts:
        genome_file = results_dir / f"{cohort}_{args.suffix}.genome"
        if not genome_file.exists():
            print(f"\nWARNING: Genome file not found: {genome_file}, skipping")
            continue

        print(f"\nProcessing: {cohort} ({args.suffix})")

        classified = classify_cohort(genome_file)
        related = classified[classified["predicted_relationship"] != "unrelated"]

        print("Total pairs: " + str(len(classified)))
        print(f"Related pairs: {len(related)}")
        print(f"Relationship distribution:")
        dist = classified["predicted_relationship"].value_counts()
        for rel_type, cnt in dist.items():
            print(f"  {rel_type}: {cnt}")

        out_path = results_dir / f"{cohort}_classified.tsv"
        cols = [
            "FID1",
            "IID1",
            "FID2",
            "IID2",
            "PI_HAT",
            "Z0",
            "Z1",
            "Z2",
            "predicted_relationship",
        ]
        avail = [c for c in cols if c in classified.columns]
        classified[avail].to_csv(out_path, sep="\t", index=False)
        print(f"Saved: {out_path}")

        # TODO probably should normalize these labels somewhere else
        if known_df is not None:
            acc = evaluate_accuracy(classified, known_df)
            has_data = not acc.empty
            if has_data:
                acc.insert(0, "cohort", cohort)
                all_acc.append(acc)
                print(f"Accuracy:")
                for _, row in acc.iterrows():
                    if row["tp"] + row["fp"] + row["fn"] > 0:
                        print(
                            f"  {row['relationship']}: "
                            f"P={row['precision']:.3f} R={row['recall']:.3f} "
                            f"F1={row['f1']:.3f} "
                            f"(TP={row['tp']} FP={row['fp']} FN={row['fn']})"
                        )

    if all_acc:
        summary = pd.concat(all_acc, ignore_index=True)
        summary_path = results_dir / "accuracy_summary.tsv"
        summary.to_csv(summary_path, sep="\t", index=False)
        print(f"\nAccuracy summary saved: {summary_path}")

        print("\nOverall accuracy by relationship type:")
        overall = summary.groupby("relationship")[["tp", "fp", "fn"]].sum()
        overall["precision"] = overall["tp"] / (overall["tp"] + overall["fp"]).replace(
            0, float("nan")
        )
        overall["recall"] = overall["tp"] / (overall["tp"] + overall["fn"]).replace(
            0, float("nan")
        )
        overall["f1"] = (
            2
            * overall["precision"]
            * overall["recall"]
            / (overall["precision"] + overall["recall"]).replace(0, float("nan"))
        )
        print(overall.round(4).to_string())
    else:
        print("\nNo accuracy data to summarize.")


if __name__ == "__main__":
    main()
