# Pig Metagenomics - Quality Control and Host Removal Pipeline

This repository documents the preprocessing workflow for shotgun metagenomic sequencing data from pig gut samples using KneadData.

## Project Overview

- **Samples**: 30 pig gut metagenomic samples
- **Sequencing**: Paired-end Illumina (Nextera library prep)
- **Raw data size**: ~15 GB
- **Goal**: Remove low-quality reads, adapters, and pig host contamination

## Directory Structure

```
Complex_Diet/
├── data/
│   ├── raw/
│   │   └── raw_data/           # Raw FASTQ files (*_R1_001.fastq.gz, *_R2_001.fastq.gz)
│   ├── clean/                  # KneadData output (cleaned FASTQ files)
│   │   ├── *_cleaned_paired_1.fastq
│   │   ├── *_cleaned_paired_2.fastq
│   │   ├── *_cleaned_unmatched_1.fastq
│   │   ├── *_cleaned_unmatched_2.fastq
│   │   ├── *_pig_host_bowtie2_contam*.fastq  # Host reads (removed)
│   │   ├── fastqc/             # FastQC reports
│   │   └── kneaddata_summary.tsv
│   └── refs/
│       └── pig_host_index/     # Bowtie2 index for pig genome
```

## Environment Setup

### Conda Environment

```bash
# Create environment with mamba
mamba create -n pig_metagenomics python=3.9
mamba activate pig_metagenomics

# Install required tools
mamba install -c bioconda kneaddata bowtie2 trimmomatic fastqc
```

### Interactive Session (for testing)

```bash
srun -p msismall -N 1 -n 1 --cpus-per-task=16 --mem=64G --time=24:00:00 --pty bash
```

## Pipeline: KneadData

KneadData performs two main steps:
1. **Quality trimming & adapter removal** (Trimmomatic)
2. **Host decontamination** (Bowtie2 alignment to pig genome)

### SLURM Batch Script

```bash
#!/bin/bash
#SBATCH --job-name=knead
#SBATCH --output=kneaddata_%A_%a.out
#SBATCH --error=kneaddata_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --array=1-30

# Add environment to PATH
export PATH="/users/7/ssilvaju/.local/share/mamba/envs/pig_metagenomics/bin:$PATH"

# Fix Java memory for Trimmomatic
export _JAVA_OPTIONS="-Xmx32g"

# Define paths
RAW_DIR="/projects/standard/gomeza/ssilvaju/Complex_Diet/data/raw/raw_data"
OUT_DIR="/projects/standard/gomeza/ssilvaju/Complex_Diet/data/clean"
REF="/projects/standard/gomeza/ssilvaju/Complex_Diet/data/refs/pig_host_index/pig_host"
TRIMMO="/users/7/ssilvaju/.local/share/mamba/envs/pig_metagenomics/share/trimmomatic"
ADAPTERS="${TRIMMO}/adapters/NexteraPE-PE.fa"
KNEAD="/users/7/ssilvaju/.local/share/mamba/envs/pig_metagenomics/bin/kneaddata"
BOWTIE2_DIR="/users/7/ssilvaju/.local/share/mamba/envs/pig_metagenomics/bin"
THREADS=8

mkdir -p $OUT_DIR

# Get sample from array task
R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${RAW_DIR}/samples_R1.txt)
R2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${RAW_DIR}/samples_R2.txt)
SAMPLE=$(echo "$R1" | cut -d'_' -f1)

echo "Running Kneaddata on ${SAMPLE}"

$KNEAD \
    --input1 ${RAW_DIR}/${R1} \
    --input2 ${RAW_DIR}/${R2} \
    --reference-db ${REF} \
    --output ${OUT_DIR} \
    --output-prefix ${SAMPLE}_cleaned \
    --threads ${THREADS} \
    --trimmomatic ${TRIMMO} \
    --trimmomatic-options "ILLUMINACLIP:${ADAPTERS}:2:30:10 MINLEN:50" \
    --bowtie2 ${BOWTIE2_DIR} \
    --bowtie2-options "--very-sensitive --phred33" \
    --run-fastqc-end \
    --remove-intermediate-output \
    --serial

echo "Completed ${SAMPLE}"
```

### Key Parameters Explained

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `--reference-db` | pig_host | Bowtie2 index for host removal |
| `--trimmomatic-options` | ILLUMINACLIP:...:2:30:10 MINLEN:50 | Remove Nextera adapters, discard reads <50bp |
| `--bowtie2-options` | --very-sensitive | Stringent host detection |
| `--remove-intermediate-output` | - | Save disk space |
| `_JAVA_OPTIONS="-Xmx32g"` | - | Prevent Trimmomatic memory errors |

## Results Summary

### Generate Summary Table

```bash
kneaddata_read_count_table \
    --input /projects/standard/gomeza/ssilvaju/Complex_Diet/data/clean/ \
    --output /projects/standard/gomeza/ssilvaju/Complex_Diet/data/clean/kneaddata_summary.tsv
```

### Calculate Retention Rates

```bash
awk -F'\t' 'NR>1 && $2!="NA" {
    raw=$2; final=$12;
    trim_ret=$4/raw*100;
    host_ret=$8/$4*100;
    total_ret=final/raw*100;
    printf "%s\tTrim:%.1f%%\tHost:%.1f%%\tTotal:%.1f%%\n", $1, trim_ret, host_ret, total_ret
}' kneaddata_summary.tsv
```

### Our Results

| Metric | Range | Interpretation |
|--------|-------|----------------|
| Trimming retention | 32-72% | Variable quality across samples |
| Host removal retention | 92-97% | Low host contamination (good for fecal) |
| Total retention | 30-70% | Acceptable for downstream analysis |

#### Sample-level Results

| Sample | Trimming | Host Removal | Total |
|--------|----------|--------------|-------|
| 2827 | 60.1% | 93.6% | 56.2% |
| 2828 | 72.0% | 96.2% | 69.3% |
| 2830 | 68.3% | 94.9% | 64.9% |
| 2839 | 71.3% | 96.3% | 68.7% |
| 2859 | 72.3% | 96.1% | 69.5% |
| 2870 | 32.6% | 92.7% | **30.2%** ⚠️ |

### Acceptable Thresholds

| Step | Acceptable | Concern |
|------|------------|---------|
| After trimming | >50% | <40% |
| After host removal | >90% | <80% |
| Overall | >40% | <30% |

## Output Files

### Per Sample

| File | Description |
|------|-------------|
| `*_cleaned_paired_1.fastq` | Clean forward reads (paired) |
| `*_cleaned_paired_2.fastq` | Clean reverse reads (paired) |
| `*_cleaned_unmatched_1.fastq` | Orphan forward reads |
| `*_cleaned_unmatched_2.fastq` | Orphan reverse reads |
| `*_pig_host_bowtie2_*_contam.fastq` | Removed host reads |
| `*_cleaned.log` | Processing log |

### Summary

| File | Description |
|------|-------------|
| `kneaddata_summary.tsv` | Read counts at each step |
| `fastqc/` | Quality reports |

## Troubleshooting

### Common Errors

**1. Java OutOfMemoryError**
```
Caused by: java.lang.OutOfMemoryError: Java heap space
```
Solution: Add `export _JAVA_OPTIONS="-Xmx32g"` to script

**2. Bowtie2 not found**
```
ERROR: Unable to find bowtie2
```
Solution: Use `--bowtie2 /path/to/bin/` (directory, not executable)

**3. Disk quota exceeded**
```
OSError: [Errno 122] Disk quota exceeded
```
Solution: Move data to `/projects/` instead of home directory

### Check for Failed Samples

```bash
# Count successful samples
ls *_cleaned_paired_1.fastq | wc -l

# Check error logs
cat kneaddata_*.err | grep -i error
```

## Quality Control Checks

### Verify Host Removal Worked

```bash
# Check for contaminant files (should exist)
ls *_contam*.fastq | head -5

# Contaminant files contain host reads that were removed
```

### Check Sample Quality

```bash
# View FastQC summary
cat fastqc/*_fastqc/summary.txt
```

## Next Steps

After KneadData preprocessing:

1. **Taxonomic profiling**: MetaPhlAn, Kraken2, or similar
2. **Functional profiling**: HUMAnN, SUPER-FOCUS
3. **Assembly**: MEGAHIT, metaSPAdes
4. **Statistical analysis**: Use compositional methods (ALDEx2, ANCOM-II)

## References

- [KneadData Documentation](https://huttenhower.sph.harvard.edu/kneaddata/)
- [Trimmomatic Manual](http://www.usadellab.org/cms/?page=trimmomatic)
- [Bowtie2 Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

## Notes
