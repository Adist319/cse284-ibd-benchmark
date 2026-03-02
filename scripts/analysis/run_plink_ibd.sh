#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJ="$(cd "$SCRIPT_DIR/../.." && pwd)"

PLINK="$PROJ/tools/plink19/plink"
DATA="$PROJ/data/processed"
OUT="$PROJ/results/plink"
BENCH="$OUT/benchmarks.tsv"

COHORTS=(admixed homogeneous trios)
# THRESHOLDS=(0.05 0.10 0.125 0.15 0.20 0.25)  # 0.125 was too many combos
THRESHOLDS=(0.05 0.10 0.15 0.20 0.25)

mkdir -p "$OUT"

if [ ! -x "$PLINK" ]; then
  echo "ERROR: PLINK not found at $PLINK"
  exit 1
fi

echo -e "cohort\tmin_threshold\twall_time_sec\tpeak_memory_kb\tnum_pairs" > "$BENCH"

# macOS time output is different from linux
parse_time_output() {
  local tf="$1"
  local wt pm

  wt=$(grep "real" "$tf" | head -1 | sed 's/[^0-9.]//g' | head -1)
  if [ -z "$wt" ]; then
    wt=$(head -1 "$tf" | sed 's/[^0-9.]//g')
  fi

  # macOS gives bytes not KB
  pm=$(grep "maximum resident set size" "$tf" | sed 's/[^0-9]//g')
  if [ -n "$pm" ]; then
    pm=$((pm / 1024))
  else
    pm="NA"
  fi

  echo "${wt:-NA}\t${pm}"
}

count_pairs() {
  local gf="$1"
  if [ -f "$gf" ]; then
    tail -n +2 "$gf" | wc -l | tr -d ' '
  else
    echo "0"
  fi
}

run_plink_genome() {
  local cohort="$1"
  local threshold="$2"
  local suffix="$3"
  local input="${DATA}/${cohort}_chr22_pruned"
  local output="$OUT/${cohort}_${suffix}"
  local tlog="$OUT/.time_${cohort}_${suffix}.log"

  if [ ! -f "${input}.bed" ]; then
    echo "  WARNING: ${input}.bed not found, skipping"
    return
  fi

  echo "running: ${cohort} (${suffix})..."
  # echo "DEBUG: input=$input threshold=$threshold"

  local args=(
    --bfile "$input"
    --genome
    --out "$output"
    --allow-extra-chr
  )

  if [ "$threshold" != "none" ]; then
    args+=(--min "$threshold")
  fi

  /usr/bin/time -l "$PLINK" "${args[@]}" \
    2> "$tlog"

  local timing
  timing=$(parse_time_output "$tlog")
  local wt pm
  wt=$(echo -e "$timing" | cut -f1)
  pm=$(echo -e "$timing" | cut -f2)

  local np
  np=$(count_pairs "${output}.genome")

  local thr_label="${threshold}"
  if [ "$threshold" = "none" ]; then
    thr_label="default"
  fi

  echo -e "${cohort}\t${thr_label}\t${wt}\t${pm}\t${np}" >> "$BENCH"
  printf "    Pairs: %s | Time: %ss | Mem: %sKB\n" "$np" "$wt" "$pm"

  rm -f "$tlog"
}

echo "PLINK IBD analysis"
echo "output: $OUT"
echo ""

echo "running defaults..."
for cohort in "${COHORTS[@]}"; do
  run_plink_genome "$cohort" "none" "default"
done
echo ""

# TODO could probably skip some of these thresholds, 0.05 and 0.10 produce a lot of pairs
echo "running --min PI_HAT sweep..."
for cohort in "${COHORTS[@]}"; do
  for threshold in "${THRESHOLDS[@]}"; do
    run_plink_genome "$cohort" "$threshold" "min${threshold}"
  done
done
echo ""

echo "done. benchmarks: $BENCH"
echo ""
column -t -s $'\t' "$BENCH"
