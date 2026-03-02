# converts phased VCF -> GERMLINE PED/MAP format
# usage: python3 vcf_to_germline_ped.py <input.vcf.gz> <genetic_map> <output_prefix> [panel_file]
import sys
import gzip
from bisect import bisect_left


def load_genetic_map(map_file, chrom):
    positions = []
    cm_vals = []
    chrom_clean = chrom.replace('chr', '')
    with open(map_file) as f:
        for line in f:
            p = line.strip().split()
            if len(p) < 4:
                continue
            if p[0].replace('chr', '') != chrom_clean:
                continue
            positions.append(int(p[3]))
            cm_vals.append(float(p[2]))

    ret = (positions, cm_vals)
    return ret


def interpolate_cm(pos, positions, cm_vals):
    if not positions:
        return 0.0

    idx = bisect_left(positions, pos)

    if idx == 0:
        return cm_vals[0]
    if idx >= len(positions):
        return cm_vals[-1]

    bp1, bp2 = positions[idx - 1], positions[idx]
    cm1, cm2 = cm_vals[idx - 1], cm_vals[idx]
    if (bp2 == bp1):
        return cm1

    frac = (pos - bp1) / (bp2 - bp1)
    return cm1 + frac * (cm2 - cm1)


# loads sex info from 1000 genomes panel file
def load_panel_sex(panel_file):
    sex_map = {}
    if not panel_file:
        return sex_map
    # was using pd.read_csv here but it's overkill for this
    try:
        with open(panel_file) as f:
            for line in f:
                p = line.strip().split('\t')
                if len(p) >= 4 and p[0] != 'sample':
                    gender = p[3].strip()
                    sex = '1' if gender == 'male' else ('2' if gender == 'female' else '0')
                    sex_map[p[0]] = sex
    except FileNotFoundError:
        pass
    return (sex_map)


def main():
    if len(sys.argv) < 4:
        print("usage: vcf_to_germline_ped.py <vcf> <map> <out_prefix> [panel]",
              file=sys.stderr)
        sys.exit(1)

    vcf_file = sys.argv[1]
    map_file = sys.argv[2]
    out_prefix = sys.argv[3]
    panel_file = sys.argv[4] if len(sys.argv) > 4 else None

    sex_map = load_panel_sex(panel_file)

    opener = gzip.open if vcf_file.endswith('.gz') else open
    samples = []
    variants = []
    # gt_matrix = np.zeros(...)  # tried numpy first, too slow for this size
    sample_gts = None

    print(f"Reading VCF: {vcf_file}", file=sys.stderr)
    n_vars = 0

    # TODO should validate the VCF header format
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith("#CHROM"):
                parts = line.strip().split('\t')
                samples = parts[9:]
                sample_gts = [[] for _ in range(len(samples))]
                # print(f"DEBUG: {len(samples)} samples found")
                print(f"Found {len(samples)} samples", file=sys.stderr)
                continue

            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            snp_id = parts[2] if parts[2] != '.' else f'{chrom}:{pos}'
            ref = parts[3]
            alt = parts[4]
            alleles = [ref, alt]

            # FIXME handle multi-allelic sites properly
            variants.append((chrom, pos, snp_id))

            for si, gt_field in enumerate(parts[9:]):
                gt = gt_field.split(':')[0]
                if '|' in gt:
                    i1, i2 = gt.split('|')
                elif '/' in gt:
                    i1, i2 = gt.split('/')
                else:
                    i1, i2 = '.', '.'

                if i1 == '.' or i2 == '.':
                    sample_gts[si].append(('0', '0'))
                else:
                    sample_gts[si].append((alleles[int(i1)], alleles[int(i2)]))

            n_vars += 1
            if n_vars % 50000 == 0:  # progress every 50k
                print(f"Processed {n_vars} variants...", file=sys.stderr)

    print("Total: {} variants".format(n_vars), file=sys.stderr)

    if not variants:
        print("ERROR: No variants found in VCF", file=sys.stderr)
        sys.exit(1)

    target_chrom = variants[0][0]
    print(f"Loading genetic map for {target_chrom}...", file=sys.stderr)
    gmap_pos, gmap_cm = load_genetic_map(map_file, target_chrom)
    print(f"Genetic map: {len(gmap_pos)} positions loaded", file=sys.stderr)

    # write .map
    print(f"Writing MAP file...", file=sys.stderr)
    with open(out_prefix + '.map', 'w') as f:
        for chrom, pos, snp_id in variants:
            cm = interpolate_cm(pos, gmap_pos, gmap_cm)
            f.write(f'{chrom}\t{snp_id}\t{cm:.6f}\t{pos}\n')

    # write .ped - FID IID PAT MAT SEX PHENO then tab-separated alleles
    print(f"Writing PED file...", file=sys.stderr)
    with open(out_prefix + '.ped', 'w') as f:
        for si, sample in enumerate(samples):
            sex = sex_map.get(sample, '0')
            fields = [sample, sample, '0', '0', sex, '-9']

            for a1, a2 in sample_gts[si]:
                fields.append(a1)
                fields.append(a2)

            f.write(' '.join(fields) + '\n')

            if (si + 1) % 100 == 0:
                print(f"Written {si + 1}/{len(samples)} samples...", file=sys.stderr)

    print(f"Done: {len(samples)} samples x {n_vars} variants", file=sys.stderr)


if __name__ == '__main__':
    main()
