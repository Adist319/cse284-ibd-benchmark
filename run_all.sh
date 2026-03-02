#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"

echo "IBD benchmarking pipeline"
echo ""

# quick sanity check
"$ROOT/tools/plink19/plink" --version 2>&1 | head -1
"$ROOT/tools/plink2/plink2" --version 2>&1 | head -1
"$ROOT/tools/germline/germline" --help 2>&1 | head -1
bcftools --version | head -1
echo

echo "downloading data..."
bash "$ROOT/scripts/preprocessing/download_data.sh"
echo

echo "preprocessing..."
bash "$ROOT/scripts/preprocessing/preprocess.sh"
echo ""

echo "running PLINK IBD..."
bash "$ROOT/scripts/analysis/run_plink_ibd.sh"
echo ""

# germline trios takes a while
echo "running GERMLINE IBD..."
bash "$ROOT/scripts/analysis/run_germline_ibd.sh"
echo

echo "classifying germline segments..."
source "$ROOT/.venv/bin/activate"

MATCH_ARGS=""
for cohort in trios admixed homogeneous; do
    mf="$ROOT/results/germline/${cohort}_default.match"
    if [ -f "$mf" ]; then
        MATCH_ARGS="$MATCH_ARGS ${cohort}:${mf}"
    fi
done

if [ -n "$MATCH_ARGS" ]; then
    python3 "$ROOT/scripts/analysis/classify_germline_relationships.py" \
        --match-files $MATCH_ARGS \
        --known "$ROOT/data/processed/known_relationships.tsv" \
        --output-dir "$ROOT/results/germline"
fi
echo ""

echo "comparing tools..."
for cohort in trios admixed homogeneous; do
    pf="$ROOT/results/plink/${cohort}_default.genome"
    gf="$ROOT/results/germline/${cohort}_default.match"
    kf="$ROOT/data/processed/known_relationships.tsv"

    if [ -f "$pf" ] && [ -f "$gf" ]; then
        python3 "$ROOT/scripts/analysis/compare_tools.py" \
            --cohort "$cohort" \
            --plink-genome "$pf" \
            --germline-match "$gf" \
            --known-rels "$kf"
    else
        echo "skipping $cohort - missing files"
    fi
done

echo ""
echo "done. results in $ROOT/results/"
