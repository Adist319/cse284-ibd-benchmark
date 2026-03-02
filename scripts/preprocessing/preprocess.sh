#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
RAW_DIR="${PROJECT_DIR}/data/raw"
PROC_DIR="${PROJECT_DIR}/data/processed"
REF_DIR="${PROJECT_DIR}/data/reference"
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"

PLINK="${PROJECT_DIR}/tools/plink19/plink"
PLINK2="${PROJECT_DIR}/tools/plink2/plink2"

VCF="${RAW_DIR}/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
PANEL="${REF_DIR}/integrated_call_samples_v3.20130502.ALL.panel"
PED_FILE="${REF_DIR}/integrated_call_samples.20130502.ALL.ped"
RELATED_FILE="${REF_DIR}/20140625_related_individuals.txt"
GMAP="${REF_DIR}/plink.chr22.GRCh38.map"

mkdir -p "${PROC_DIR}"

echo "checking prereqs..."

for tool in bcftools tabix python3; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: $tool not found in PATH"
        exit 1
    fi
done

if [ ! -x "${PLINK}" ]; then
    echo "ERROR: PLINK 1.9 not found at ${PLINK}"
    exit 1
fi

for f in "${VCF}" "${PANEL}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

if [ ! -f "${GMAP}" ]; then
    echo "WARNING: Genetic map not found, GERMLINE prep will be skipped"
    SKIP_GERMLINE=true
else
    SKIP_GERMLINE=false
fi

# the download sometimes returns HTML 404 instead of actual data
PED_VALID=false
if [ -f "${PED_FILE}" ]; then
    if head -1 "${PED_FILE}" | grep -q "<!DOCTYPE\|<html"; then
        echo "WARNING: Pedigree file looks corrupted (got HTML). Using related file instead."
    else
        PED_VALID=true
    fi
else
    echo "WARNING: No pedigree file, using related individuals file"
fi

VCF_SAMPLES=$(bcftools query -l "${VCF}")
echo "  VCF has $(echo "${VCF_SAMPLES}" | wc -l | tr -d ' ') samples"

echo ""
echo "creating sample lists..."

echo "${VCF_SAMPLES}" | sort > "${PROC_DIR}/.vcf_samples_sorted.tmp"

# admixed: AMR superpop + admixed african
echo "  extracting admixed (PUR, CLM, MXL, PEL, ASW, ACB)..."
awk -F'\t' 'NR > 1 && ($2 == "PUR" || $2 == "CLM" || $2 == "MXL" || $2 == "PEL" || $2 == "ASW" || $2 == "ACB") {print $1}' \
    "${PANEL}" | sort | comm -12 - "${PROC_DIR}/.vcf_samples_sorted.tmp" \
    > "${PROC_DIR}/samples_admixed.txt"

ADMIXED_N=$(wc -l < "${PROC_DIR}/samples_admixed.txt")
echo "  Admixed: ${ADMIXED_N} samples"

echo "  extracting homogeneous (CEU, GBR, TSI)..."
awk -F'\t' 'NR > 1 && ($2 == "CEU" || $2 == "GBR" || $2 == "TSI") {print $1}' \
    "${PANEL}" | sort | comm -12 - "${PROC_DIR}/.vcf_samples_sorted.tmp" \
    > "${PROC_DIR}/samples_homogeneous.txt"

HOMO_N=$(wc -l < "${PROC_DIR}/samples_homogeneous.txt")
echo "  Homogeneous: ${HOMO_N} samples"

echo "  extracting trios from pedigree..."

if [ "${PED_VALID}" = true ]; then
    python3 - "${PED_FILE}" "${PROC_DIR}" "${PROC_DIR}/.vcf_samples_sorted.tmp" << 'PARSE_PED'
import sys
from collections import defaultdict

ped_file = sys.argv[1]
out_dir = sys.argv[2]
vcf_samples_file = sys.argv[3]

with open(vcf_samples_file) as f:
    vcf_samples = set(line.strip() for line in f)

trios = []
all_trio_samples = set()

with open(ped_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#') or line.lower().startswith('family'):
            continue
        parts = line.split('\t')
        if len(parts) < 5:
            parts = line.split()
        if len(parts) < 5:
            continue

        fid, iid, pat, mat = parts[0], parts[1], parts[2], parts[3]
        has_father = pat != '0' and pat != ''
        has_mother = mat != '0' and mat != ''

        if has_father or has_mother:
            members_in_vcf = iid in vcf_samples
            if has_father:
                members_in_vcf = members_in_vcf and pat in vcf_samples
            if has_mother:
                members_in_vcf = members_in_vcf and mat in vcf_samples

            if members_in_vcf:
                trios.append((iid, pat if has_father else None, mat if has_mother else None, fid))
                all_trio_samples.add(iid)
                if has_father:
                    all_trio_samples.add(pat)
                if has_mother:
                    all_trio_samples.add(mat)

with open(f'{out_dir}/samples_trios.txt', 'w') as f:
    for s in sorted(all_trio_samples):
        f.write(s + '\n')

seen_pairs = set()
with open(f'{out_dir}/known_relationships.tsv', 'w') as f:
    f.write('sample1\tsample2\trelationship_type\n')
    for child, father, mother, fid in trios:
        if father:
            pair = tuple(sorted([child, father]))
            if pair not in seen_pairs:
                f.write(f'{child}\t{father}\tparent-child\n')
                seen_pairs.add(pair)
        if mother:
            pair = tuple(sorted([child, mother]))
            if pair not in seen_pairs:
                f.write(f'{child}\t{mother}\tparent-child\n')
                seen_pairs.add(pair)

    parent_to_children = defaultdict(list)
    for child, father, mother, fid in trios:
        if father and mother:
            parent_key = tuple(sorted([father, mother]))
            parent_to_children[parent_key].append(child)

    for parent_key, children in parent_to_children.items():
        if len(children) > 1:
            for i in range(len(children)):
                for j in range(i + 1, len(children)):
                    pair = tuple(sorted([children[i], children[j]]))
                    if pair not in seen_pairs:
                        f.write(f'{children[i]}\t{children[j]}\tfull-sibling\n')
                        seen_pairs.add(pair)

print(f'  Found {len(trios)} children with parents, {len(all_trio_samples)} total samples')
PARSE_PED

else
    # fallback when ped file is missing/corrupted
    echo "  Using related individuals file as fallback..."
    if [ -f "${RELATED_FILE}" ]; then
        python3 - "${RELATED_FILE}" "${PROC_DIR}" "${PROC_DIR}/.vcf_samples_sorted.tmp" << 'PARSE_RELATED'
import sys

related_file = sys.argv[1]
out_dir = sys.argv[2]
vcf_samples_file = sys.argv[3]

with open(vcf_samples_file) as f:
    vcf_samples = set(line.strip() for line in f)

all_samples = set()
relationships = []

with open(related_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('Sample'):
            continue
        parts = line.split('\t')
        if len(parts) < 4:
            continue

        sample = parts[0].strip()
        reason = parts[3].strip()

        if reason.startswith('trio:'):
            parents = reason[5:].split(',')
            if len(parents) == 2:
                father, mother = parents[0].strip(), parents[1].strip()
                if sample in vcf_samples and father in vcf_samples and mother in vcf_samples:
                    all_samples.update([sample, father, mother])
                    relationships.append((sample, father, 'parent-child'))
                    relationships.append((sample, mother, 'parent-child'))

        elif reason.startswith('Sibling:'):
            sibling = reason[8:].strip()
            if sample in vcf_samples and sibling in vcf_samples:
                all_samples.update([sample, sibling])
                relationships.append((sample, sibling, 'full-sibling'))

        elif reason.startswith('Parent:'):
            child = reason[7:].strip()
            if sample in vcf_samples and child in vcf_samples:
                all_samples.update([sample, child])
                relationships.append((sample, child, 'parent-child'))

with open(f'{out_dir}/samples_trios.txt', 'w') as f:
    for s in sorted(all_samples):
        f.write(s + '\n')

seen_pairs = set()
with open(f'{out_dir}/known_relationships.tsv', 'w') as f:
    f.write('sample1\tsample2\trelationship_type\n')
    for s1, s2, rel in relationships:
        pair = tuple(sorted([s1, s2]))
        if pair not in seen_pairs:
            f.write(f'{s1}\t{s2}\t{rel}\n')
            seen_pairs.add(pair)

print(f'  Found {len(all_samples)} related samples, {len(seen_pairs)} unique relationships')
PARSE_RELATED
    else
        echo "  WARNING: No relationship data available. Creating empty trios files."
        echo -e "sample1\tsample2\trelationship_type" > "${PROC_DIR}/known_relationships.tsv"
        touch "${PROC_DIR}/samples_trios.txt"
    fi
fi

rm -f "${PROC_DIR}/.vcf_samples_sorted.tmp"

TRIOS_N=$(wc -l < "${PROC_DIR}/samples_trios.txt")
REL_N=$(( $(wc -l < "${PROC_DIR}/known_relationships.tsv") - 1 ))
echo "  Trios: ${TRIOS_N} samples, ${REL_N} known relationships"

echo ""
echo "filtering VCFs by population..."

for cohort in admixed homogeneous trios; do
    SAMPLES="${PROC_DIR}/samples_${cohort}.txt"
    OUT_VCF="${PROC_DIR}/${cohort}_chr22.vcf.gz"

    N=$(wc -l < "${SAMPLES}")
    if [ "${N}" -eq 0 ]; then
        echo "  Skipping ${cohort}: no samples"
        continue
    fi

    if [ -f "${OUT_VCF}" ] && [ -f "${OUT_VCF}.tbi" ]; then
        echo "  ${cohort}: already exists ($(bcftools query -l "${OUT_VCF}" | wc -l | tr -d ' ') samples)"
        continue
    fi

    echo "  Filtering ${cohort} (${N} samples)..."
    bcftools view \
        --samples-file "${SAMPLES}" \
        --force-samples \
        --min-ac 1 \
        -m2 -M2 -v snps \
        "${VCF}" \
    | bcftools view \
        --min-af 0.05:minor \
        -Oz -o "${OUT_VCF}"

    tabix -p vcf "${OUT_VCF}"

    NVAR=$(bcftools index -n "${OUT_VCF}")
    echo "  ${cohort}: ${NVAR} variants, ${N} samples"
done

echo ""
echo "preparing PLINK files..."

for cohort in admixed homogeneous trios; do
    CVCF="${PROC_DIR}/${cohort}_chr22.vcf.gz"
    CBED="${PROC_DIR}/${cohort}_chr22"
    CPRUNED="${PROC_DIR}/${cohort}_chr22_pruned"

    if [ ! -f "${CVCF}" ]; then
        echo "  Skipping ${cohort}: no filtered VCF"
        continue
    fi

    if [ ! -f "${CBED}.bed" ]; then
        echo "  converting ${cohort} VCF to PLINK..."
        "${PLINK}" \
            --vcf "${CVCF}" \
            --make-bed \
            --out "${CBED}" \
            --allow-extra-chr \
            --const-fid 0 \
            --memory 4096 \
            2>&1 | tail -1
    else
        echo "  ${cohort}: PLINK binary exists"
    fi

    # update .fam with pedigree info
    echo "  updating ${cohort} .fam..."
    if [ "${PED_VALID}" = true ]; then
        python3 - "${PED_FILE}" "${CBED}.fam" << 'UPDATE_FAM_PED'
import sys

ped_file = sys.argv[1]
fam_file = sys.argv[2]

ped_info = {}
with open(ped_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#') or line.lower().startswith('family'):
            continue
        parts = line.split('\t')
        if len(parts) < 5:
            parts = line.split()
        if len(parts) >= 5:
            fid, iid, pat, mat, sex = parts[0], parts[1], parts[2], parts[3], parts[4]
            if sex.lower() in ('male', 'm', '1'):
                sex_code = '1'
            elif sex.lower() in ('female', 'f', '2'):
                sex_code = '2'
            else:
                sex_code = '0'
            ped_info[iid] = (fid, pat, mat, sex_code)

lines = []
updated = 0
with open(fam_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 6:
            iid = parts[1]
            if iid in ped_info:
                fid, pat, mat, sex = ped_info[iid]
                parts[0] = fid
                parts[2] = pat
                parts[3] = mat
                parts[4] = sex
                updated += 1
        lines.append('\t'.join(parts))

with open(fam_file, 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f'  Updated {updated} samples with pedigree info')
UPDATE_FAM_PED
    else
        python3 - "${PANEL}" "${CBED}.fam" << 'UPDATE_FAM_PANEL'
import sys

panel_file = sys.argv[1]
fam_file = sys.argv[2]

panel_info = {}
with open(panel_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 4 and parts[0] != 'sample':
            gender = parts[3].strip()
            sex_code = '1' if gender == 'male' else ('2' if gender == 'female' else '0')
            panel_info[parts[0]] = sex_code

lines = []
with open(fam_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 6:
            iid = parts[1]
            if iid in panel_info:
                parts[4] = panel_info[iid]
        lines.append('\t'.join(parts))

with open(fam_file, 'w') as f:
    f.write('\n'.join(lines) + '\n')
UPDATE_FAM_PANEL
    fi

    # LD pruning
    if [ ! -f "${CPRUNED}.bed" ]; then
        echo "  LD pruning ${cohort}..."
        "${PLINK}" \
            --bfile "${CBED}" \
            --indep-pairwise 50 5 0.2 \
            --out "${CBED}_ldprune" \
            --allow-extra-chr \
            --memory 4096 \
            2>&1 | tail -1

        "${PLINK}" \
            --bfile "${CBED}" \
            --extract "${CBED}_ldprune.prune.in" \
            --make-bed \
            --out "${CPRUNED}" \
            --allow-extra-chr \
            --memory 4096 \
            2>&1 | tail -1

        # plink resets .fam on --make-bed so we have to redo this
        if [ "${PED_VALID}" = true ]; then
            python3 - "${PED_FILE}" "${CPRUNED}.fam" << 'UPDATE_FAM_PED2'
import sys

ped_file = sys.argv[1]
fam_file = sys.argv[2]

ped_info = {}
with open(ped_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#') or line.lower().startswith('family'):
            continue
        parts = line.split('\t')
        if len(parts) < 5:
            parts = line.split()
        if len(parts) >= 5:
            fid, iid, pat, mat, sex = parts[0], parts[1], parts[2], parts[3], parts[4]
            if sex.lower() in ('male', 'm', '1'):
                sex_code = '1'
            elif sex.lower() in ('female', 'f', '2'):
                sex_code = '2'
            else:
                sex_code = '0'
            ped_info[iid] = (fid, pat, mat, sex_code)

lines = []
with open(fam_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 6:
            iid = parts[1]
            if iid in ped_info:
                fid, pat, mat, sex = ped_info[iid]
                parts[0] = fid
                parts[2] = pat
                parts[3] = mat
                parts[4] = sex
        lines.append('\t'.join(parts))

with open(fam_file, 'w') as f:
    f.write('\n'.join(lines) + '\n')
UPDATE_FAM_PED2
        fi

        rm -f "${CBED}_ldprune.prune.in" "${CBED}_ldprune.prune.out" \
              "${CBED}_ldprune.log" "${CBED}_ldprune.nosex"

        PRUNED_N=$(wc -l < "${CPRUNED}.bim")
        ORIG_N=$(wc -l < "${CBED}.bim")
        echo "  ${cohort}: ${ORIG_N} -> ${PRUNED_N} SNPs after pruning"
    else
        echo "  ${cohort}: pruned files exist"
    fi
done

echo ""
echo "preparing GERMLINE PED/MAP files..."

if [ "${SKIP_GERMLINE}" = true ]; then
    echo "  Skipped (no genetic map). Run download_data.sh first."
    echo ""
    echo "done (GERMLINE files skipped)."
    exit 0
fi

CONVERTER="${SCRIPTS_DIR}/vcf_to_germline_ped.py"
cat > "${CONVERTER}" << 'PYTHON_CONVERTER'
#!/usr/bin/env python3
"""Convert phased VCF to GERMLINE PED/MAP format."""
import sys
import gzip
from bisect import bisect_left


def load_genetic_map(map_file, target_chrom):
    positions = []
    cm_values = []
    target_clean = target_chrom.replace('chr', '')

    with open(map_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom_clean = parts[0].replace('chr', '')
            if chrom_clean != target_clean:
                continue
            bp = int(parts[3])
            cm = float(parts[2])
            positions.append(bp)
            cm_values.append(cm)

    return positions, cm_values


def interpolate_cm(pos, positions, cm_values):
    if not positions:
        return 0.0

    idx = bisect_left(positions, pos)

    if idx == 0:
        return cm_values[0]
    if idx >= len(positions):
        return cm_values[-1]

    bp1, bp2 = positions[idx - 1], positions[idx]
    cm1, cm2 = cm_values[idx - 1], cm_values[idx]

    if bp2 == bp1:
        return cm1

    frac = (pos - bp1) / (bp2 - bp1)
    return cm1 + frac * (cm2 - cm1)


def load_panel_sex(panel_file):
    sex_map = {}
    if not panel_file:
        return sex_map
    try:
        with open(panel_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4 and parts[0] != 'sample':
                    gender = parts[3].strip()
                    sex_code = '1' if gender == 'male' else ('2' if gender == 'female' else '0')
                    sex_map[parts[0]] = sex_code
    except FileNotFoundError:
        pass
    return sex_map


def main():
    if len(sys.argv) < 4:
        print("Usage: python3 vcf_to_germline_ped.py <input.vcf.gz> <genetic_map> <output_prefix> [panel_file]",
              file=sys.stderr)
        sys.exit(1)

    vcf_file = sys.argv[1]
    map_file = sys.argv[2]
    out_prefix = sys.argv[3]
    panel_file = sys.argv[4] if len(sys.argv) > 4 else None

    sex_map = load_panel_sex(panel_file)

    opener = gzip.open if vcf_file.endswith('.gz') else open
    samples = []
    variant_info = []
    sample_gts = None

    print(f"  Reading VCF: {vcf_file}", file=sys.stderr)
    n_variants = 0

    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                samples = parts[9:]
                sample_gts = [[] for _ in range(len(samples))]
                print(f"  Found {len(samples)} samples", file=sys.stderr)
                continue

            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            snp_id = parts[2] if parts[2] != '.' else f'{chrom}:{pos}'
            ref = parts[3]
            alt = parts[4]
            alleles = [ref, alt]

            variant_info.append((chrom, pos, snp_id))

            for s_idx, gt_field in enumerate(parts[9:]):
                gt = gt_field.split(':')[0]
                if '|' in gt:
                    idx1, idx2 = gt.split('|')
                elif '/' in gt:
                    idx1, idx2 = gt.split('/')
                else:
                    idx1, idx2 = '.', '.'

                if idx1 == '.' or idx2 == '.':
                    sample_gts[s_idx].append(('0', '0'))
                else:
                    sample_gts[s_idx].append((alleles[int(idx1)], alleles[int(idx2)]))

            n_variants += 1
            if n_variants % 50000 == 0:
                print(f"  Processed {n_variants} variants...", file=sys.stderr)

    print(f"  Total: {n_variants} variants", file=sys.stderr)

    if not variant_info:
        print("  ERROR: No variants found in VCF", file=sys.stderr)
        sys.exit(1)

    target_chrom = variant_info[0][0]
    print(f"  Loading genetic map for {target_chrom}...", file=sys.stderr)
    gmap_pos, gmap_cm = load_genetic_map(map_file, target_chrom)
    print(f"  Genetic map: {len(gmap_pos)} positions loaded", file=sys.stderr)

    print(f"  Writing MAP file...", file=sys.stderr)
    with open(out_prefix + '.map', 'w') as f:
        for chrom, pos, snp_id in variant_info:
            cm = interpolate_cm(pos, gmap_pos, gmap_cm)
            f.write(f'{chrom}\t{snp_id}\t{cm:.6f}\t{pos}\n')

    print(f"  Writing PED file...", file=sys.stderr)
    with open(out_prefix + '.ped', 'w') as f:
        for s_idx, sample in enumerate(samples):
            sex = sex_map.get(sample, '0')
            fields = [sample, sample, '0', '0', sex, '-9']

            for a1, a2 in sample_gts[s_idx]:
                fields.append(a1)
                fields.append(a2)

            f.write(' '.join(fields) + '\n')

            if (s_idx + 1) % 100 == 0:
                print(f"  Written {s_idx + 1}/{len(samples)} samples...", file=sys.stderr)

    print(f"  Done: {len(samples)} samples x {n_variants} variants", file=sys.stderr)


if __name__ == '__main__':
    main()
PYTHON_CONVERTER

chmod +x "${CONVERTER}"

for cohort in admixed homogeneous trios; do
    CVCF="${PROC_DIR}/${cohort}_chr22.vcf.gz"
    OUT="${PROC_DIR}/${cohort}_chr22_phased"

    if [ ! -f "${CVCF}" ]; then
        echo "  Skipping ${cohort}: no filtered VCF"
        continue
    fi

    if [ -f "${OUT}.ped" ] && [ -f "${OUT}.map" ]; then
        echo "  ${cohort}: GERMLINE files exist"
        continue
    fi

    echo "  converting ${cohort} to GERMLINE PED/MAP..."
    # TODO trios cohort takes forever here (~2hrs), maybe subsample?
    python3 "${CONVERTER}" \
        "${CVCF}" \
        "${GMAP}" \
        "${OUT}" \
        "${PANEL}"

    if [ -f "${OUT}.ped" ] && [ -f "${OUT}.map" ]; then
        PED_N=$(wc -l < "${OUT}.ped")
        MAP_N=$(wc -l < "${OUT}.map")
        echo "  ${cohort}: ${PED_N} samples, ${MAP_N} variants"
    else
        echo "  ERROR: Failed to create GERMLINE files for ${cohort}"
        exit 1
    fi
done

echo ""
echo "preprocessing done. output in ${PROC_DIR}/:"
echo ""
echo "  sample lists:"
for cohort in admixed homogeneous trios; do
    f="${PROC_DIR}/samples_${cohort}.txt"
    if [ -f "$f" ]; then
        echo "    samples_${cohort}.txt: $(wc -l < "$f") samples"
    fi
done
echo "    known_relationships.tsv: $(( $(wc -l < "${PROC_DIR}/known_relationships.tsv") - 1 )) relationships"
echo ""
echo "  filtered VCFs (biallelic SNPs, MAF >= 0.05):"
for cohort in admixed homogeneous trios; do
    f="${PROC_DIR}/${cohort}_chr22.vcf.gz"
    if [ -f "$f" ]; then
        n=$(bcftools index -n "$f" 2>/dev/null || echo "?")
        echo "    ${cohort}_chr22.vcf.gz: ${n} variants"
    fi
done
echo ""
echo "  PLINK files:"
for cohort in admixed homogeneous trios; do
    f="${PROC_DIR}/${cohort}_chr22.bed"
    fp="${PROC_DIR}/${cohort}_chr22_pruned.bed"
    if [ -f "$f" ]; then
        echo "    ${cohort}_chr22.{bed,bim,fam}: $(wc -l < "${PROC_DIR}/${cohort}_chr22.bim") SNPs"
    fi
    if [ -f "$fp" ]; then
        echo "    ${cohort}_chr22_pruned.{bed,bim,fam}: $(wc -l < "${PROC_DIR}/${cohort}_chr22_pruned.bim") SNPs"
    fi
done
if [ "${SKIP_GERMLINE}" = false ]; then
    echo ""
    echo "  GERMLINE PED/MAP:"
    for cohort in admixed homogeneous trios; do
        f="${PROC_DIR}/${cohort}_chr22_phased.ped"
        if [ -f "$f" ]; then
            echo "    ${cohort}_chr22_phased.{ped,map}"
        fi
    done
fi
echo ""
echo "done."
