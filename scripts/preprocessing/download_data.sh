#!/bin/bash
set -euo pipefail

DATA_DIR="$(cd "$(dirname "$0")/../../data/raw" && pwd)"
REF_DIR="$(cd "$(dirname "$0")/../../data/reference" && pwd)"
mkdir -p "$DATA_DIR" "$REF_DIR"

BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp"

download() {
  local url="$1"
  local dest="$2"
  local fname
  fname=$(basename "$url")
  if [ -f "$dest/$fname" ]; then
    echo "  already exists: $fname"
  else
    echo "  downloading $fname..."
    # wget "$url" -O "$dest/$fname"  # switched to curl, wget not on macOS by default
    curl -sL --retry 3 -o "$dest/$fname" "$url"
    echo "  done: $(ls -lh "$dest/$fname" | awk '{print $5}')"
  fi
}

echo "downloading pedigree + sample metadata..."

download "${BASE_URL}/release/20130502/integrated_call_samples.20130502.ALL.ped" "$REF_DIR"
download "${BASE_URL}/release/20130502/20140625_related_individuals.txt" "$REF_DIR"
download "${BASE_URL}/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" "${REF_DIR}"

echo "downloading 1kG 30x phased VCF (chr22)..."

PHASED_URL="${BASE_URL}/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
VCF="1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
download "${PHASED_URL}/${VCF}" "$DATA_DIR"
download "${PHASED_URL}/${VCF}.tbi" "$DATA_DIR"

echo 'downloading genetic map...'
GMAP_URL="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"
if [ ! -f "$REF_DIR/plink.chr22.GRCh38.map" ]; then
  curl -sL --retry 3 -o "$REF_DIR/plink.GRCh38.map.zip" "$GMAP_URL"
  cd "$REF_DIR" && unzip -o plink.GRCh38.map.zip && cd -
fi

# TODO should probably verify checksums
echo ""
echo "done. contents:"
ls -lh "$DATA_DIR"/
echo ""
ls -lh ${REF_DIR}/
