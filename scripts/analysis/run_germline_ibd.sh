#!/bin/bash
set -euo pipefail

PROJ="$(cd "$(dirname "$0")/../.." && pwd)"
GERMLINE="$PROJ/tools/germline/germline"
DATA="$PROJ/data/processed"
OUT="$PROJ/results/germline"
BENCH="$OUT/benchmarks.tsv"

COHORTS=(admixed homogeneous trios)
MIN_M_VALUES=(1 2 3 5 10)
# ERR_VALUES=(0 1 2 3 4)  # 3 didn't add much, dropped it
ERR_VALUES=(0 1 2 4)

mkdir -p "$OUT"

if [[ ! -x "$GERMLINE" ]]; then
    echo "error: germline not found at $GERMLINE" >&2
    exit 1
fi

printf "cohort\tparam_name\tparam_value\twall_time_sec\tpeak_memory_kb\tnum_segments\tnum_pairs\n" > "$BENCH"

# macOS /usr/bin/time -l format
parse_time_output() {
    local tf="$1"
    local wt pm_bytes

    wt=$(grep 'real' "$tf" | awk '{print $1}')
    pm_bytes=$(grep 'peak memory footprint' "${tf}" | awk '{print $1}')

    local pm_kb
    if [[ -n "${pm_bytes}" ]]; then
        pm_kb=$(( pm_bytes / 1024 ))
    else
        pm_kb="NA"
    fi
    echo "$wt" "$pm_kb"
}

count_match_stats() {
    local mf=$1
    local nseg npairs

    if [[ -f "$mf" ]]; then
        nseg=$(wc -l < "$mf" | tr -d ' ')
        npairs=$(awk '{print $1"_"$2"_"$3"_"$4}' "$mf" | sort -u | wc -l | tr -d ' ')
    else
        nseg=0
        npairs=0
    fi
    echo "$nseg" "$npairs"
}

run_germline() {
    local cohort=$1
    local pname=$2
    local pval=$3
    local out_prefix=$4
    shift 4
    local extra=("$@")

    local ped="$DATA/${cohort}_chr22_phased.ped"
    local map="$DATA/${cohort}_chr22_phased.map"

    if [[ ! -f "$ped" ]]; then
        echo "  skipping $cohort: no ped file" >&2
        return 0
    fi
    if [[ ! -f "$map" ]]; then
        echo "  skipping $cohort: no map file" >&2
        return 0
    fi

    local tf
    tf=$(mktemp "$OUT/time_XXXXXX.txt")

    echo "  running: $cohort $pname=$pval"

    /usr/bin/time -l "$GERMLINE" \
        -haploid \
        "${extra[@]}" \
        -input "$ped" "$map" \
        -output "$out_prefix" \
        2> "$tf" || {
            echo "  WARNING: germline failed for $cohort $pname=$pval" >&2
            rm -f "$tf"
            return 0
        }

    local match="${out_prefix}.match"
    read -r wt pm_kb <<< "$(parse_time_output "$tf")"
    read -r nseg npairs <<< "$(count_match_stats "$match")"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$cohort" "$pname" "$pval" \
        "$wt" "$pm_kb" \
        "$nseg" "$npairs" >> "$BENCH"

    rm -f "$tf"
}

echo "GERMLINE IBD analysis"
echo "output: $OUT"
echo ""

echo "running defaults..."
for cohort in "${COHORTS[@]}"; do
    run_germline "$cohort" "default" "default" \
        "$OUT/${cohort}_default"
done

echo ""
echo "running -min_m sweep..."
for cohort in "${COHORTS[@]}"; do
    for min_m in "${MIN_M_VALUES[@]}"; do
        run_germline "$cohort" "min_m" "$min_m" \
            "$OUT/${cohort}_minm${min_m}" \
            -min_m "$min_m"
    done
done

# TODO could skip some of these combos, takes a while for trios
echo ""
echo "running error tolerance sweep..."
for cohort in "${COHORTS[@]}"; do
    for err in "${ERR_VALUES[@]}"; do
        run_germline "$cohort" "err" "$err" \
            "$OUT/${cohort}_err${err}" \
            -err_hom "$err" -err_het "$err"
    done
done

echo ""
echo "done. benchmarks: $BENCH"
