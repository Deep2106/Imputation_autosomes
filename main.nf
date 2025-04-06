// =====================
// main.nf
// =====================

nextflow.enable.dsl=2

params.input_dir      = '/path/to/input_vcfs'
params.genetic_map    = '/path/to/genetic_maps'
params.reference_dir  = '/path/to/reference_panel'
params.output_dir     = 'imputation_results'
params.chroms         = (1..22)

workflow {
  Channel.from(params.chroms)
    .filter { chrom ->
      !file("${params.output_dir}/imputed_chr${chrom}.vcf.gz").exists()
    }
    | phase_shapeit5
    | bcf_to_vcf
    | generate_chunks
    | impute_chunks
    | merge_chunks
    | cleanup_chunks
    | qc_bcftools_stats

  qc_bcftools_stats.out.collect() | multiqc_report
}

process phase_shapeit5 {
  tag "chr${chrom}"
  publishDir "${params.output_dir}", mode: 'copy'

  input:
  val chrom

  output:
  tuple val(chrom), path("phased_chr${chrom}.bcf")

  script:
  """
  module load shapeit5
  shapeit5_phase_common \\
    --input ${params.input_dir}/input_chr${chrom}.vcf.gz \\
    --map ${params.genetic_map}/chr${chrom}.b38.gmap.gz \\
    --region ${chrom} \\
    --output phased_chr${chrom}.bcf \\
    --thread 8
  """
}

process bcf_to_vcf {
  tag "chr${chrom}"
  publishDir "${params.output_dir}", mode: 'copy'

  input:
  tuple val(chrom), path(phased_bcf)

  output:
  tuple val(chrom), path("phased_chr${chrom}.vcf.gz")

  script:
  """
  module load bcftools
  bcftools view ${phased_bcf} -Oz -o phased_chr${chrom}.vcf.gz --threads 8
  tabix -p vcf phased_chr${chrom}.vcf.gz
  """
}

process generate_chunks {
  tag "chr${chrom}"
  publishDir "${params.output_dir}", mode: 'copy'

  input:
  tuple val(chrom), path(vcf)

  output:
  tuple val(chrom), path("coordinates_chr${chrom}.txt")

  script:
  """
  module load imp5Chunker
  imp5Chunker \\
    --h ${params.reference_dir}/chr${chrom}_biallelic_snps.vcf.gz \\
    --g ${vcf} \\
    --r ${chrom} \\
    --o coordinates_chr${chrom}.txt
  """
}

process impute_chunks {
  tag "chr${chrom}"
  publishDir "${params.output_dir}", mode: 'copy'

  input:
  tuple val(chrom), path(vcf)
  path coords

  output:
  tuple val(chrom), path("imputed_chr${chrom}_*.vcf.gz")

  script:
  """
  module load impute5
  while read CHUNK_ID CHR FULL_REGION IMPUTE_REGION SIZE TARGET_SNP REF_SNP; do
    impute5 \\
      --g ${vcf} \\
      --m ${params.genetic_map}/chr${chrom}.b38.gmap.gz \\
      --h ${params.reference_dir}/chr${chrom}_biallelic_snps.vcf.gz \\
      --r ${IMPUTE_REGION} \\
      --buffer-region ${FULL_REGION} \\
      --o imputed_chr${chrom}_\${CHUNK_ID}.vcf.gz \\
      --threads 8
  done < ${coords}
  """
}

process merge_chunks {
  tag "chr${chrom}"
  publishDir "${params.output_dir}", mode: 'copy'

  input:
  tuple val(chrom), path(imputed_vcfs)

  output:
  tuple val(chrom), path("imputed_chr${chrom}.vcf.gz")

  script:
  """
  module load bcftools
  ls imputed_chr${chrom}_*.vcf.gz | sort -V > imputed_chunks_chr${chrom}.txt
  bcftools concat -n -f imputed_chunks_chr${chrom}.txt -Oz -o imputed_chr${chrom}.vcf.gz
  """
}

process cleanup_chunks {
  tag "chr${chrom}"

  input:
  val chrom

  script:
  """
  rm -f ${params.output_dir}/imputed_chr${chrom}_*.vcf.gz
  rm -f ${params.output_dir}/imputed_chunks_chr${chrom}.txt
  """
}

process qc_bcftools_stats {
  tag "chr${chrom}"
  publishDir "${params.output_dir}/qc", mode: 'copy'

  input:
  tuple val(chrom), path(vcf)

  output:
  path "bcftools_stats_chr${chrom}.txt"

  script:
  """
  module load bcftools
  bcftools stats ${vcf} > bcftools_stats_chr${chrom}.txt
  """
}

process multiqc_report {
  publishDir "${params.output_dir}/qc", mode: 'copy'

  input:
  path(stats_files)

  output:
  path "multiqc_report.html"

  script:
  """
  module load multiqc
  multiqc . --filename multiqc_report.html
  """
}
