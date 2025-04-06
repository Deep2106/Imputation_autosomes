#!/bin/bash

# Define input and output directories
VCF_DIR="<PATH>/1000G/GRCh38_VCF"  # Update this path to your VCF folder
OUT_DIR="<PATH>/1000G/GRCh38_VCF/LEGEND_HAP"  # Output folder for hap/legend files

# Create output directory if it doesn't exist
mkdir -p ${OUT_DIR}

# Chromosomes: Autosomes (1-22) + X
CHROMS=($(seq 1 22) X)

# Load required modules (if running on HPC)
# module load bcftools   # Uncomment if using an HPC environment

# Process each chromosome
for CHR in ${CHROMS[@]}; do
    echo "Processing Chromosome: ${CHR}"

    # Define input VCF and output file paths
    VCF_FILE="${VCF_DIR}/ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    CLEANED_VCF="${OUT_DIR}/chr${CHR}_biallelic_snps.vcf.gz"
    HAP_FILE="${OUT_DIR}/chr${CHR}.hap.gz"
    LEGEND_FILE="${OUT_DIR}/chr${CHR}.legend.gz"

    # Check if VCF file exists
    if [[ ! -f ${VCF_FILE} ]]; then
        echo "Warning: VCF file for chr${CHR} not found! Skipping..."
        continue
    fi

    # Step 1: Remove multi-allelic sites (split into biallelic)
    echo "Step 1: Splitting multi-allelic sites for chr${CHR} (8 threads)... change accordingly"
    bcftools norm -m -both --threads 8 -Oz -o ${OUT_DIR}/chr${CHR}_biallelic.vcf.gz ${VCF_FILE}
    bcftools index ${OUT_DIR}/chr${CHR}_biallelic.vcf.gz  # No --threads here

    # Step 2: (Optional) Remove indels, keep only SNPs
    echo "Step 2: Removing indels for chr${CHR} (8 threads)...change accordingly"
    bcftools view -v snps --threads 8 -Oz -o ${CLEANED_VCF} ${OUT_DIR}/chr${CHR}_biallelic.vcf.gz
    bcftools index ${CLEANED_VCF}  # No --threads here

    # Step 3: Generate HAP file
    echo "Step 3: Generating HAP file for chr${CHR}..."
    bcftools query -f '%CHROM %POS %ID %REF %ALT [ %GT]\n' ${CLEANED_VCF} | gzip > ${HAP_FILE}

    # Step 4: Generate LEGEND file
    echo "Step 4: Generating LEGEND file for chr${CHR}..."
    bcftools query -f '%ID %CHROM %POS %REF %ALT\n' ${CLEANED_VCF} | gzip > ${LEGEND_FILE}

    echo "Completed processing for Chromosome: ${CHR}"
done

echo "Conversion completed for all chromosomes!"
