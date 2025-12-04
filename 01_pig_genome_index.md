# Building Pig Genome Bowtie2 Index for Host Removal

This guide documents how to download the pig reference genome (Sus scrofa) and build a Bowtie2 index for host decontamination in metagenomic analyses.

## Overview

- **Genome**: Sus scrofa (pig) assembly Sscrofa11.1
- **Source**: Ensembl Release 115
- **Purpose**: Remove host contamination from gut metagenomics data
- **Total size**: ~1.3-1.5 GB compressed

## Directory Structure

```
Complex_Diet/
└── data/
    └── refs/
        ├── pig_genome.fa                    # Combined genome FASTA
        └── pig_host_index/                  # Bowtie2 index files
            ├── pig_host.1.bt2
            ├── pig_host.2.bt2
            ├── pig_host.3.bt2
            ├── pig_host.4.bt2
            ├── pig_host.rev.1.bt2
            └── pig_host.rev.2.bt2
```

## Step 1: Download Pig Genome

### Navigate to Reference Directory

```bash
cd ~/Complex_Diet/data/refs
mkdir -p pig_host_index
```

### Download All Chromosomes

Downloads chromosomes 1-18, X, Y, and mitochondrial DNA from Ensembl:

```bash
for f in {1..18} X Y MT; do 
    wget https://ftp.ensembl.org/pub/release-115/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.primary_assembly.$f.fa.gz
done
```

### Verify Download

```bash
ls -lh Sus_scrofa.Sscrofa11.1.dna.primary_assembly*.fa.gz
```

Expected output: 21 files totaling ~1.3-1.5 GB

### Check File Integrity

```bash
# Count downloaded files (should be 21)
ls Sus_scrofa.Sscrofa11.1.dna.primary_assembly*.fa.gz | wc -l

# Check for any incomplete downloads (should show no errors)
for f in Sus_scrofa.Sscrofa11.1.dna.primary_assembly*.fa.gz; do
    gzip -t "$f" && echo "$f OK" || echo "$f CORRUPTED"
done
```

## Step 2: Combine and Decompress

### Concatenate All Chromosome Files

```bash
cat Sus_scrofa.Sscrofa11.1.dna.primary_assembly*.fa.gz > pig_genome.fa.gz
```

### Decompress

```bash
gunzip pig_genome.fa.gz
```

### Verify Combined Genome

```bash
# Check file size (should be ~2.5 GB uncompressed)
ls -lh pig_genome.fa

# Count sequences (should show 21 chromosome entries)
grep -c "^>" pig_genome.fa
```

## Step 3: Build Bowtie2 Index

### Option A: Interactive Session (for testing)

```bash
# Start interactive session
srun -p msismall -N 1 -n 1 --cpus-per-task=8 --mem=32G --time=24:00:00 --pty bash

# Activate environment
conda activate pig_metagenomics

# Build index
cd ~/Complex_Diet/data/refs
bowtie2-build pig_genome.fa pig_host_index/pig_host
```

### Option B: SLURM Batch Job (recommended)

Create the batch script:

```bash
cat > ~/Complex_Diet/data/refs/build_bowtie2_index.sbatch << 'EOF'
#!/bin/bash
#SBATCH --job-name=bowtie2_pig_index
#SBATCH --output=bowtie2_pig_index.out
#SBATCH --error=bowtie2_pig_index.err
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email@umn.edu

echo "Starting Bowtie2 index build at: $(date)"

# Add environment to PATH
export PATH="/users/7/ssilvaju/.local/share/mamba/envs/pig_metagenomics/bin:$PATH"

# Move to reference folder
cd ~/Complex_Diet/data/refs

# Create output directory
mkdir -p pig_host_index

# Build the Bowtie2 index
bowtie2-build --threads 8 pig_genome.fa pig_host_index/pig_host

echo "Index build finished at: $(date)"
EOF
```

Submit the job:

```bash
sbatch ~/Complex_Diet/data/refs/build_bowtie2_index.sbatch
```

Monitor progress:

```bash
squeue -u $USER
tail -f ~/Complex_Diet/data/refs/bowtie2_pig_index.out
```

## Step 4: Verify Index

### Check Index Files Exist

```bash
ls -lh ~/Complex_Diet/data/refs/pig_host_index/
```

Expected output:

```
pig_host.1.bt2
pig_host.2.bt2
pig_host.3.bt2
pig_host.4.bt2
pig_host.rev.1.bt2
pig_host.rev.2.bt2
```

### Test Index with Bowtie2

```bash
bowtie2-inspect -s ~/Complex_Diet/data/refs/pig_host_index/pig_host | head -20
```

## Step 5: Clean Up (Optional)

After confirming the index works, remove intermediate files to save space:

```bash
cd ~/Complex_Diet/data/refs

# Remove individual chromosome files
rm Sus_scrofa.Sscrofa11.1.dna.primary_assembly*.fa.gz

# Optionally compress the combined genome
gzip pig_genome.fa
```

## Using the Index with KneadData

Reference the index in your KneadData command:

```bash
kneaddata \
    --input1 sample_R1.fastq.gz \
    --input2 sample_R2.fastq.gz \
    --reference-db ~/Complex_Diet/data/refs/pig_host_index/pig_host \
    --output output_dir \
    ...
```

## Troubleshooting

### Download Fails

If wget fails, try with different options:

```bash
wget --no-check-certificate https://ftp.ensembl.org/pub/release-115/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.primary_assembly.1.fa.gz
```

Or use curl:

```bash
for f in {1..18} X Y MT; do 
    curl -O https://ftp.ensembl.org/pub/release-115/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.primary_assembly.$f.fa.gz
done
```

### Index Build Runs Out of Memory

Increase memory allocation:

```bash
#SBATCH --mem=64G
```

### Bowtie2 Not Found

Ensure the environment is activated or use full path:

```bash
/users/7/ssilvaju/.local/share/mamba/envs/pig_metagenomics/bin/bowtie2-build pig_genome.fa pig_host
```

## Resource Requirements

| Step | Time | Memory | Disk Space |
|------|------|--------|------------|
| Download | ~30 min | Minimal | ~1.5 GB |
| Decompress | ~5 min | Minimal | ~2.5 GB |
| Build index | 2-4 hours | 16-32 GB | ~4 GB |
| **Total** | ~5 hours | 32 GB | ~8 GB |

## References

- [Ensembl Sus scrofa](https://www.ensembl.org/Sus_scrofa/Info/Index)
- [Bowtie2 Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [Sscrofa11.1 Assembly](https://www.ncbi.nlm.nih.gov/assembly/GCF_000003025.6/)

## Version Information

```bash
# Record versions for reproducibility
bowtie2 --version
echo "Genome: Sus_scrofa.Sscrofa11.1 (Ensembl release 115)"
date
```
