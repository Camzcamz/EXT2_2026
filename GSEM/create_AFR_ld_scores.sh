#!/bin/bash
################################################################################
# FORMAT AFR 1KG REF & GET AFR LD SCORES FOR Genomic SEM
################################################################################

# Change your directory path
BASE=~/path/to/your/directory

# Set PLINK tool path (download plink 2 from https://www.cog-genomics.org/plink/)
export PLINK_DIR=$BASE/tools
export PATH="$PATH:$PLINK_DIR"

# Reference data directory
REF_DIR=$BASE/ref

# Decompress .zst 1KG files
plink2 --zst-decompress "$REF_DIR/all_phase3.pgen.zst" "$REF_DIR/all_phase3.pgen"
plink2 --zst-decompress "$REF_DIR/all_phase3.pvar.zst" "$REF_DIR/all_phase3.pvar"

# File paths
PGEN_RAW=$REF_DIR/all_phase3
PSAM_RAW=$REF_DIR/all_phase3.psam
REL_IDS=$REF_DIR/deg2_phase3.king.cutoff.out.id
TEMP_REF=$REF_DIR/temp
mkdir -p "$TEMP_REF"

#####################################
# AFR Filtering  #####
#####################################
# Subset to atgc SNPs
# Calculate allele frequencies
# Extract useful position info from .pvar for downstream joining
# Combine AFREQ and position data
# Final formatting to generate unified frequency file
# RSID file (RSID, CHROM, POS)
# ALLELE file (ChrPosID, REF, ALT, AF_ALT_AFR_1000G)
# GET BFILE FOR POPCORN AND LDSCORES

AFR_MAF=0.005
AFR_MAF_NAME=maf005

grep -vFf "$REL_IDS" "$PSAM_RAW" | awk '$5 == "AFR" {print $1}' > "$TEMP_REF/all_phase3_king2d_afr.txt"

plink2 --pfile "$PGEN_RAW" \
  --keep "$TEMP_REF/all_phase3_king2d_afr.txt" \
  --snps-only just-acgt \
  --maf "$AFR_MAF" \
  --make-pgen \
  --out "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}"

plink2 --pfile "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}" \
  --freq \
  --out "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}"

awk 'BEGIN {OFS="\t"} !/^##/{print $2, $3}' "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}.pvar" > "$TEMP_REF/positions_afr.txt"

paste "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}.afreq" "$TEMP_REF/positions_afr.txt" > "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}_freq.txt"

awk 'BEGIN {OFS="\t"} \
	NR==1 {print "ChrPosID", "Chr", "RSID", "Pos", "REF", "ALT", "AF_ALT_AFR_1000G"} \
	NR>1 {print $1":"$7, $1, $2, $7, $3, $4, $5}' \
	"$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}_freq.txt" > "$TEMP_REF/1KG_AFR_AF_${AFR_MAF_NAME}.txt"


for CHR in {1..22}; do
  plink2 --pfile "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}" \
    --chr $CHR \
    --max-alleles 2 \
    --make-bed \
    --out "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}_chr${CHR}"
done

cd "$REF_DIR"

# Set environment variables
export TEMP_REF=$REF_DIR/temp
export LDSCORE_OUT_DIR=$BASE/ref/1KG.ldscore.AFR
export AFR_MAF_NAME=maf005
mkdir -p "$LDSCORE_OUT_DIR"

# Generate .bed/.bim/.fam files per chromosome
for CHR in {1..22}; do
  plink2 --pfile "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}" \
    --chr $CHR \
    --max-alleles 2 \
    --make-bed \
    --out "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}_chr${CHR}"
done

# Run LDSC on each chromosome
cd "$BASE/tools/ldsc/" || exit
source activate ldsc
for CHR in {1..22}; do
  python ldsc.py \
    --bfile "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}_chr${CHR}" \
    --l2 \
    --ld-wind-kb 1000 \
    --out "$LDSCORE_OUT_DIR/chr${CHR}"
done

# Rename LDSC output files for Genomic SEM
cd "$LDSCORE_OUT_DIR" || exit
for file in *; do
  new_name=$(echo "$file" | sed 's/^chr//')
  mv "$file" "$new_name"
done
