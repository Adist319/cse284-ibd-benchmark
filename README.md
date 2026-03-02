# Benchmarking PLINK and GERMLINE for IBD-Based Relative Finding

**CSE 284 - Adrian Ong (A69033975)**

## Overview

This project benchmarks two widely-used Identity by Descent (IBD) detection tools - **PLINK** (`--genome`) and **GERMLINE** - for identifying genetic relatives in diverse human populations.

We use the [1000 Genomes Project 30x high-coverage dataset](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38) (3,202 samples including 602 confirmed parent-child trios) to evaluate:

- IBD detection accuracy against 1,204 known parent-child relationships
- Computational performance (runtime, peak memory) across cohort sizes
- Parameter sensitivity (minimum segment length, mismatch tolerance, PI_HAT thresholds)
- Population structure effects comparing admixed vs. homogeneous cohorts

The analysis is limited to chromosome 22 for computational tractability.

## Key Findings

**PLINK** is extremely fast (< 1 second for most cohorts) and classifies parent-child relationships with near-perfect accuracy on the trios cohort (precision=0.987, recall=1.0, F1=0.993). It only outputs summary IBD statistics (PI_HAT, Z0/Z1/Z2) - no segment-level information.

**GERMLINE** detects IBD segments with full position and length information. It correctly detects ~98.8% of parent-child pairs when using proportion-based thresholds calibrated for chr22. However, its fixed whole-genome segment-length thresholds completely fail on single-chromosome data (F1=0.01), which is an important limitation to be aware of. GERMLINE is also substantially slower - 28 minutes for the 1,793-sample trios cohort vs. under 1 second for PLINK.

**Runtime comparison (default parameters):**

| Cohort | N Samples | PLINK | GERMLINE | Speedup |
|--------|-----------|-------|----------|---------|
| Admixed | 504 | 0.09 sec | 27 sec | ~300x |
| Homogeneous | 297 | 0.07 sec | 12 sec | ~170x |
| Trios | 1,793 | 0.90 sec | 1,677 sec | ~1,860x |

**Accuracy on known parent-child pairs (trios cohort):**

| Tool | Precision | Recall | F1 |
|------|-----------|--------|----|
| PLINK (Z0/Z1/Z2 thresholds) | 0.987 | 1.000 | 0.993 |
| GERMLINE (proportion-based, chr22) | ~0.99 | ~0.988 | ~0.989 |
| GERMLINE (whole-genome thresholds) | 0.353 | 0.005 | 0.010 |

The key methodological takeaway: GERMLINE's default classification thresholds assume whole-genome input. On chr22 alone (~2% of the genome), expected IBD lengths are proportionally smaller, so those thresholds don't work.

## Project Structure

```
.
|-- data/
|   |-- raw/                    # raw 1000 Genomes VCFs
|   |-- processed/              # population-filtered, LD-pruned files
|   `-- reference/              # pedigree, population panel, genetic maps
|-- notebooks/
|   `-- ibd_benchmarking_analysis.ipynb   # main analysis notebook
|-- scripts/
|   |-- preprocessing/
|   |   |-- download_data.sh            # download 1000 Genomes data
|   |   |-- preprocess.sh               # filter, convert, LD-prune
|   |   |-- vcf_to_germline_ped.py      # phase-preserving VCF->PED converter
|   |   `-- vcf_to_germline_fast.py     # faster two-pass version
|   `-- analysis/
|       |-- run_plink_ibd.sh                    # PLINK --genome + parameter sweeps
|       |-- run_germline_ibd.sh                 # GERMLINE + parameter sweeps
|       |-- classify_plink_relationships.py     # IBD -> relationship classifier
|       |-- classify_germline_relationships.py  # segment -> relationship classifier
|       `-- compare_tools.py                    # head-to-head comparison + figures
|-- results/
|   |-- plink/                  # PLINK .genome outputs + benchmarks
|   |-- germline/               # GERMLINE .match outputs + benchmarks
|   |-- comparison/             # cross-tool comparison results
|   `-- figures/                # generated plots
|-- tools/                      # PLINK and GERMLINE binaries
|-- run_all.sh                  # master pipeline script
|-- requirements.txt
`-- README.md
```

## Dependencies

- **Python 3.10+**
- **bcftools** (`brew install bcftools` or conda)
- **C++ compiler** (for compiling GERMLINE from source)
- ~500 MB disk space for chr22 data

Python packages (installed via requirements.txt):
- pandas, numpy, matplotlib, seaborn, networkx

## Installation

```bash
# clone the repo
git clone https://github.com/Adist319/cse284-ibd-benchmark.git
cd cse284-ibd-benchmark

# set up Python environment
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# download PLINK binaries (macOS ARM example)
mkdir -p tools
cd tools
curl -sL "https://s3.amazonaws.com/plink1-assets/plink_mac_20231018.zip" -o plink19.zip
unzip plink19.zip -d plink19
curl -sL "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_arm64_20250116.zip" -o plink2.zip
unzip plink2.zip -d plink2
cd ..

# compile GERMLINE from source
git clone https://github.com/gusevlab/germline.git tools/germline
cd tools/germline && make all; cd ../..
```

For Linux, replace the macOS PLINK URLs with the appropriate Linux builds from [plink1.9](https://www.cog-genomics.org/plink/) and [plink2](https://www.cog-genomics.org/plink/2.0/).

## Running the Full Pipeline

```bash
# download data, preprocess, run all analyses
bash run_all.sh
```

Or run individual steps:

```bash
bash scripts/preprocessing/download_data.sh   # download chr22 VCF + pedigree
bash scripts/preprocessing/preprocess.sh       # filter, LD-prune, split cohorts
bash scripts/analysis/run_plink_ibd.sh         # PLINK IBD + parameter sweep
bash scripts/analysis/run_germline_ibd.sh      # GERMLINE IBD + parameter sweep
python scripts/analysis/classify_plink_relationships.py --cohort trios
python scripts/analysis/classify_germline_relationships.py --cohort trios
python scripts/analysis/compare_tools.py --cohort trios
```

Note: the GERMLINE VCF-to-PED conversion for the trios cohort (1,793 samples) takes roughly 2 hours. The admixed and homogeneous cohorts take about 5-10 minutes each.

## Quick Test Example

To run a quick test on the homogeneous cohort (smallest, ~297 samples). This assumes you've already run the preprocessing and IBD steps (or are using precomputed results):

```bash
source .venv/bin/activate

# classify PLINK relationships for the homogeneous cohort
python scripts/analysis/classify_plink_relationships.py \
    --cohorts homogeneous \
    --suffix default

# view the classified pairs (each row is a pair with predicted relationship)
head -5 results/plink/homogeneous_classified.tsv
```

For GERMLINE on the same cohort:

```bash
python scripts/analysis/classify_germline_relationships.py \
    --match-files homogeneous:results/germline/homogeneous_default.match \
    --output-dir results/germline/

# view segment-level output
cat results/germline/homogeneous_pairs_summary.tsv | head -10
```

Note: the `run_plink_ibd.sh` and `run_germline_ibd.sh` scripts run all three cohorts together and rewrite the benchmarks file, so run them as part of the full pipeline rather than individually.

## Dataset

**1000 Genomes Project 30x High-Coverage** (Byrska-Bishop et al., Cell 2022):
- 3,202 samples from 26 populations
- 602 confirmed parent-child trios (1,204 known relationships)
- Phased VCFs on GRCh38, chr22 used here

### Population Cohorts

| Cohort | Populations | N Samples | N SNPs (chr22) | N LD-pruned SNPs |
|--------|------------|-----------|----------------|------------------|
| Admixed | PUR, CLM, MXL, PEL, ASW, ACB | 504 | ~104K | ~9K |
| Homogeneous | CEU, GBR, TSI | 297 | ~90K | ~7K |
| Trios | all samples with known relatives | 1,793 | ~106K | ~9.5K |

## Methods

### PLINK IBD (`--genome`)

Method-of-moments estimator of pairwise IBD proportions. Outputs Z0, Z1, Z2 (probability of sharing 0, 1, or 2 alleles IBD) and PI_HAT (weighted sum). Requires LD-pruned input. Runs in O(n^2 * m) time but the constant is very small.

### GERMLINE

Seed-and-extend hashing algorithm for detecting shared IBD segments. Outputs segment coordinates, genetic length in cM, and mismatch count. Requires phased input - do not run through PLINK first since that strips phase information. Roughly O(n) per chromosome via hashing, but slower in practice on large cohorts.

### Relationship Classification

**PLINK** - thresholds on PI_HAT + Z scores:

| Relationship | Criteria |
|---|---|
| Parent-child | PI_HAT > 0.4, Z0 < 0.15 |
| Full sibling | PI_HAT > 0.35, Z2 > 0.1 |
| Second-degree | 0.17 < PI_HAT < 0.4 |
| Third-degree | 0.08 < PI_HAT < 0.17 |

**GERMLINE** - for chr22 specifically, use proportion-based thresholds (total IBD / 55 cM expected for parent-child on chr22). Fixed segment-length thresholds designed for whole-genome data will not work on a single chromosome.

## Remaining Work and Open Questions

Things I still want to get to in the last week, and stuff I'd appreciate feedback on:

- **Extend beyond chr22.** Right now everything is single-chromosome, which is the biggest limitation. GERMLINE's thresholds were designed for whole-genome data, and running on just chr22 (~2% of the genome) makes their default classification basically useless. I want to at least try chr1 + chr22 together to see if the thresholds start working better.
- **Sibling detection.** The ground truth only has parent-child relationships (from the 1000 Genomes trios), so I can't really validate sibling or second-degree detection yet. I'm not sure where to get reliable ground truth for those -- would appreciate suggestions.
- **GERMLINE VCF-to-PED speed.** The conversion for the trios cohort (1793 samples) takes ~2 hours with my Python script. There's probably a faster way to do this, maybe with bcftools directly or by chunking the VCF.
- **Better parameter sweep analysis.** I swept min-segment-length and mismatch tolerance for both tools, but haven't done a proper grid search or plotted ROC curves yet. The current parameter sweep plots are kind of basic.
- **KING comparison.** A few people have mentioned KING as another IBD tool worth benchmarking. Haven't had time to set it up yet but it would round out the comparison.

## References

- Byrska-Bishop, M., et al. (2022). High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell*, 185(18), 3426-3440.
- Purcell, S., et al. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. *AJHG*, 81(3), 559-575.
- Gusev, A., et al. (2009). Whole population, genome-wide mapping of hidden relatedness. *Genome Research*, 19(2), 318-326.
- Manichaikul, A., et al. (2010). Robust relationship inference in genome-wide association studies. *Bioinformatics*, 26(22), 2867-2873.
