# Imputation_autosomes
These repository uses shapeit5 and impute5 for phasing and imputation under nextflow 

# Genome-Wide Imputation Pipeline (Nextflow)

## 🛠 Requirements:
- All tools (`shapeit5_phase_common`, `bcftools`, `imp5Chunker`, `impute5`) must be available via module system.
- Input VCFs must be bgzipped (.vcf.gz) and tabix-indexed (.tbi).
- Genetic map (hg38) from https://github.com/odelaneau/shapeit4/blob/master/maps/genetic_maps.b38.tar.gz
- legend hap file can be created for hg38 1000G WGS dataset. Example bash script provided as extra-scripts

## 🚀 Running the Pipeline:
```bash
nextflow run main.nf -c nextflow.config
```

## ⚙ Parameters (set in nextflow.config):
- `input_dir`: Path to input VCFs (per chromosome)
- `genetic_map`: Directory containing SHAPEIT5-compatible genetic maps
- `reference_dir`: Directory containing reference haplotype panels
- `output_dir`: Where to store outputs
- `chroms`: List of chromosomes to process (default: 1-22 )

## 🔁 Behavior:
- The pipeline is **resumable**.
- Each chromosome is processed independently.
- If the final imputed VCF for a chromosome exists, it is **skipped**.
- Intermediate files are kept during processing and deleted after merge.

Happy Imputing! 🧬


## 🧪 Dry Run (Simulation Mode)
To simulate the pipeline without running any commands (great for debugging):

```bash
nextflow run main.nf -c nextflow.config -resume -dry-run
```

This will:
- Validate all inputs, channels, and processes
- Show which processes will run or be skipped
- NOT execute any tool (safe to test!)



## 📊 Full Monitoring & Reporting Options

Nextflow allows you to generate useful monitoring outputs with every run.

### ✅ Generate All Reports (Recommended)
This will produce an HTML report, timeline chart, and detailed resource log:
```bash
nextflow run main.nf -c nextflow.config -resume \
  --with-report execution_report.html \
  --with-trace trace.tsv \
  --with-timeline timeline.html
```

### 📈 Output Files:
- `execution_report.html`: High-level HTML summary of the entire run
- `trace.tsv`: Detailed tabular log of all process executions
- `timeline.html`: Interactive timeline showing when each process ran

### 🔁 Combine With Dry Run:
You can also simulate and preview what would run without executing anything:
```bash
nextflow run main.nf -c nextflow.config -resume -dry-run
```

These features are fully supported and help track, debug, and optimize your workflow effectively.


## 🧬 Post-Imputation Quality Control (QC)

### 🔎 Per-Chromosome QC with BCFtools
Each final imputed VCF file is evaluated using:
```bash
bcftools stats output_imputed_<CHR>.vcf.gz > bcftools_stats_<CHR>.txt
```

These stats include:
- Number of SNPs, indels
- INFO fields and missingness
- Imputation quality summary

---

### 📊 Aggregated QC Report with MultiQC
All BCFtools stats are compiled into a single interactive HTML report:
```bash
multiqc . --filename multiqc_report.html
```

Output:
- `multiqc_report.html`: Viewable in any browser
- Found under: `imputation_results/qc/`

---

## ✅ Full Workflow Output Structure

```
imputation_results/
├── TCD_322_phased.chr<CHR>.bcf
├── TCD_322_phased.chr<CHR>.vcf.gz
├── output_imputed_<CHR>.vcf.gz
├── qc/
│   ├── bcftools_stats_<CHR>.txt
│   └── multiqc_report.html
```

This ensures your pipeline is production-ready, reproducible, and QC-verified.
