# faster vcf->ped/map converter, uses byte arrays instead of storing strings
# two-pass: first build map, then read genotypes into compact arrays
import gzip
import sys
from bisect import bisect_left


def load_genetic_map(map_file, chrom):
    # was using .split('\t') but whitespace works better for map files
    positions, cm_vals = [], []
    chrom_clean = chrom.replace("chr", "")
    with open(map_file) as f:
        for line in f:
            p = line.strip().split()
            if len(p) < 4:
                continue
            if p[0].replace("chr", "") != chrom_clean:
                continue
            positions.append(int(p[3]))
            cm_vals.append(float(p[2]))
    return positions, cm_vals

def interpolate_cm(pos, positions, cm_vals):
    if not positions:
        return 0.0
    idx = bisect_left(positions, pos)
    if idx == 0:
        return (cm_vals[0])
    if idx >= len(positions):
        return cm_vals[-1]
    bp1, bp2 = positions[idx - 1], positions[idx]
    cm1, cm2 = cm_vals[idx - 1], cm_vals[idx]
    if bp2 == bp1:
        return cm1
    return cm1 + (pos - bp1) / (bp2 - bp1) * (cm2 - cm1)


# loads gender info from panel file
def load_panel_sex(panel_file):
    sex_map = {}
    if not panel_file:
        return sex_map
    try:
        with open(panel_file) as f:
            for line in f:
                p = line.strip().split('\t')
                if len(p) >= 4 and p[0] != "sample":
                    sex_map[p[0]] = '1' if p[3].strip() == 'male' else (
                        '2' if p[3].strip() == 'female' else '0')
    except FileNotFoundError:
        pass
    return sex_map


def main():
    if len(sys.argv) < 4:
        print("usage: vcf_to_germline_fast.py <vcf.gz> <genetic_map> <out_prefix> [panel]",
              file=sys.stderr)
        sys.exit(1)

    vcf_file = sys.argv[1]
    map_file = sys.argv[2]
    out_prefix = sys.argv[3]
    panel_file = sys.argv[4] if len(sys.argv) > 4 else None

    sex_map = load_panel_sex(panel_file)
    opener = gzip.open if vcf_file.endswith('.gz') else open

    # pass 1: scan VCF for header + variant info
    print(f"Pass 1: building MAP from {vcf_file}", file=sys.stderr)
    samples = []
    variants = []  # (chrom, pos, snp_id, ref, alt)

    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                samples = parts[9:]
                print("Found %d samples" % len(samples), file=sys.stderr)
                continue
            parts = line.split('\t', 6)  # only need first few cols
            chrom = parts[0]
            pos = int(parts[1])
            snp_id = parts[2] if parts[2] != '.' else f"{chrom}:{pos}"
            variants.append((chrom, pos, snp_id, parts[3], parts[4]))

    n_vars = len(variants)
    n_samp = len(samples)
    print(f"Total: {n_vars} variants, {n_samp} samples", file=sys.stderr)

    if n_vars == 0:
        print("error: no variants found in VCF", file=sys.stderr)
        sys.exit(1)

    target_chrom = variants[0][0]
    print(f"Loading genetic map for {target_chrom}...", file=sys.stderr)
    gmap_pos, gmap_cm = load_genetic_map(map_file, target_chrom)
    print(f"Genetic map: {len(gmap_pos)} positions", file=sys.stderr)

    print(f"Writing MAP file...", file=sys.stderr)
    with open(out_prefix + ".map", "w") as f:
        for chrom, pos, snp_id, ref, alt in variants:
            cm = interpolate_cm(pos, gmap_pos, gmap_cm)
            f.write(f"{chrom}\t{snp_id}\t{cm:.6f}\t{pos}\n")

    # pass 2: read genotypes into byte arrays
    # 0=ref, 1=alt, 255=missing
    # TODO could probably use numpy here but array module is good enough
    print(f"Pass 2: reading genotypes...", file=sys.stderr)

    # alleles = [ref, alt]  # old approach before switching to byte indices
    alleles_list = [(v[3], v[4]) for v in variants]

    import array
    gt1 = [array.array('B', [0] * n_vars) for _ in range(n_samp)]
    gt2 = [array.array('B', [0] * n_vars) for _ in range(n_samp)]

    vi = 0
    # TODO this is probably slow for huge files, maybe use mmap
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            fields = line.rstrip('\n').split('\t')
            for si in range(n_samp):
                tmp = fields[9 + si]
                colon = tmp.find(':')
                gt = tmp[:colon] if colon > 0 else tmp

                if '|' in gt:
                    sep = '|'
                elif '/' in gt:
                    sep = '/'
                else:
                    gt1[si][vi] = 255
                    gt2[si][vi] = 255
                    continue

                a, b = gt.split(sep)
                # print(f"DEBUG gt: {gt} a={a} b={b}")
                if (a == '.' or b == '.'):
                    gt1[si][vi] = 255
                    gt2[si][vi] = 255
                else:
                    gt1[si][vi] = int(a)
                    gt2[si][vi] = int(b)

            vi += 1
            if vi % 25000 == 0:
                print(f"Read {vi}/{n_vars} variants...", file=sys.stderr)

    print(f"Writing PED file ({n_samp} samples)...", file=sys.stderr)
    with open(out_prefix + '.ped', 'w') as f:
        for si, sample in enumerate(samples):
            sex = sex_map.get(sample, '0')
            parts = [sample, sample, '0', '0', sex, '-9']

            for j in range(n_vars):
                a1 = gt1[si][j]
                a2 = gt2[si][j]
                ref, alt = alleles_list[j]
                if a1 == 255:
                    parts.append('0')
                    parts.append('0')
                else:
                    parts.append(ref if a1 == 0 else alt)
                    parts.append(ref if a2 == 0 else alt)

            f.write(' '.join(parts) + '\n')

            if (si + 1) % 200 == 0:
                print(f"Written {si + 1}/{n_samp} samples...", file=sys.stderr)

    print(f"Done: {n_samp} samples x {n_vars} variants", file=sys.stderr)


if __name__ == "__main__":
    main()
