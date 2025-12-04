# Kraken2 Taxonomic Profiling for Pig Metagenomics

This guide documents the taxonomic profiling workflow using Kraken2 for shotgun metagenomic data from pig gut samples.

## Overview

Kraken2 is a taxonomic sequence classifier that assigns taxonomic labels to DNA sequences by examining k-mers and querying a database that maps k-mers to the lowest common ancestor (LCA) of all genomes containing that k-mer.

### Pipeline Summary

| Step | Tool | Purpose |
|------|------|---------|
| 1. Database Setup | kraken2-build | Download/build standard database |
| 2. Classification | kraken2 | Assign taxonomy to reads |
| 3. Reporting | kraken2 --report | Generate abundance tables |
| 4. Downstream | Bracken (optional) | Refine abundance estimates |

---

## Requirements

### System Requirements

- **Disk Space**: ~100 GB for standard database
- **Memory**: ~32-64 GB RAM (database loaded into memory)
- **Time**: Database build ~2-4 hours; classification ~5-15 min/sample

### Software

```bash
# Install via conda
mamba install -c bioconda kraken2 bracken
```

---

## Directory Structure

```
Complex_Diet/
├── data/
│   ├── clean/                      # KneadData output (input for Kraken2)
│   │   ├── *_cleaned_paired_1.fastq
│   │   └── *_cleaned_paired_2.fastq
│   ├── kraken2/                    # Kraken2 output
│   │   ├── outputs/                # Per-sample classification
│   │   ├── reports/                # Per-sample reports
│   │   └── merged/                 # Combined results
│   └── refs/
│       └── kraken2_db/             # Kraken2 database
│           ├── hash.k2d
│           ├── opts.k2d
│           └── taxo.k2d
```

---

## Step 1: Build Kraken2 Database

### Option A: Download Pre-built Database (Recommended)

Download from: https://benlangmead.github.io/aws-indexes/k2

```bash
# Create database directory
mkdir -p ~/Complex_Diet/data/refs/kraken2_db
cd ~/Complex_Diet/data/refs/kraken2_db

# Download standard database (~50 GB compressed)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20231009.tar.gz

# Extract
tar -xzvf k2_standard_20231009.tar.gz
rm k2_standard_20231009.tar.gz
```

### Option B: Build Standard Database from Scratch

```bash
# Set database location
DBNAME=~/Complex_Diet/data/refs/kraken2_db

# Build standard database (archaea, bacteria, viral, plasmid, human, UniVec_Core)
kraken2-build --standard --threads 16 --db $DBNAME
```

### Option C: Build Custom Database

```bash
DBNAME=~/Complex_Diet/data/refs/kraken2_db

# Step 1: Download taxonomy
kraken2-build --download-taxonomy --db $DBNAME

# Step 2: Download specific libraries
kraken2-build --download-library bacteria --db $DBNAME
kraken2-build --download-library archaea --db $DBNAME
kraken2-build --download-library viral --db $DBNAME
kraken2-build --download-library plasmid --db $DBNAME
kraken2-build --download-library UniVec_Core --db $DBNAME

# Step 3: Build database
kraken2-build --build --threads 16 --db $DBNAME

# Step 4: Clean up intermediate files (optional)
kraken2-build --clean --db $DBNAME
```

### Verify Database

```bash
# Inspect database contents
kraken2-inspect --db $DBNAME | head -20
```

---

## Step 2: Run Kraken2 Classification

### Single Sample (Testing)

```bash
DB=~/Complex_Diet/data/refs/kraken2_db
CLEAN_DIR=~/Complex_Diet/data/clean
OUT_DIR=~/Complex_Diet/data/kraken2

mkdir -p $OUT_DIR/outputs $OUT_DIR/reports

# Run Kraken2 on paired-end reads
kraken2 --db $DB \
    --paired \
    --threads 8 \
    --output $OUT_DIR/outputs/2827_kraken2.out \
    --report $OUT_DIR/reports/2827_kraken2.report \
    --report-minimizer-data \
    --use-mpa-style \
    --use-names \
    --minimum-hit-groups 2 \
    $CLEAN_DIR/2827_cleaned_paired_1.fastq \
    $CLEAN_DIR/2827_cleaned_paired_2.fastq
```

### Batch Processing (SLURM Array Job)

Create the batch script:

```bash
cat > /projects/standard/gomeza/ssilvaju/Complex_Diet/scripts/kraken2_batch.sbatch << 'EOF'
#!/bin/bash
#SBATCH --job-name=kraken2
#SBATCH --output=kraken2_%A_%a.out
#SBATCH --error=kraken2_%A_%a.err
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=8:00:00
#SBATCH --array=1-30

# Add environment to PATH
export PATH="/users/7/ssilvaju/.local/share/mamba/envs/pig_metagenomics/bin:$PATH"

# Define paths
DB="/projects/standard/gomeza/ssilvaju/Complex_Diet/data/refs/kraken2_db"
CLEAN_DIR="/projects/standard/gomeza/ssilvaju/Complex_Diet/data/clean"
OUT_DIR="/projects/standard/gomeza/ssilvaju/Complex_Diet/data/kraken2"
SAMPLE_LIST="/projects/standard/gomeza/ssilvaju/Complex_Diet/config/samples.txt"

mkdir -p $OUT_DIR/outputs $OUT_DIR/reports

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST)

echo "Processing: ${SAMPLE}"
echo "Start: $(date)"

R1="${CLEAN_DIR}/${SAMPLE}_cleaned_paired_1.fastq"
R2="${CLEAN_DIR}/${SAMPLE}_cleaned_paired_2.fastq"

# Run Kraken2
kraken2 --db $DB \
    --paired \
    --threads 16 \
    --output $OUT_DIR/outputs/${SAMPLE}_kraken2.out \
    --report $OUT_DIR/reports/${SAMPLE}_kraken2.report \
    --minimum-hit-groups 2 \
    $R1 $R2

echo "Completed: ${SAMPLE}"
echo "End: $(date)"
EOF
```

### Create Sample List

```bash
# Generate sample list from cleaned files
ls /projects/standard/gomeza/ssilvaju/Complex_Diet/data/clean/*_cleaned_paired_1.fastq | \
    xargs -n1 basename | \
    sed 's/_cleaned_paired_1.fastq//' > \
    /projects/standard/gomeza/ssilvaju/Complex_Diet/config/samples.txt

# Verify
cat /projects/standard/gomeza/ssilvaju/Complex_Diet/config/samples.txt
```

### Submit Job

```bash
sbatch ~/Complex_Diet/scripts/kraken2_batch.sbatch
```

---

## Step 3: Kraken2 Parameters Explained

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--db` | path | Path to Kraken2 database |
| `--paired` | - | Input files are paired-end reads |
| `--threads` | 8 | Number of CPU threads |
| `--output` | file | Per-read classification output |
| `--report` | file | Sample-level taxonomy report |
| `--report-minimizer-data` | - | Include minimizer counts in report |
| `--minimum-hit-groups` | 2 | Require ≥2 hit groups for classification |
| `--confidence` | 0.1 | Confidence score threshold (0-1) |

### Confidence Scoring

The confidence score is calculated as:
```
C/Q = (k-mers mapped to clade) / (total k-mers queried)
```

A threshold of 0.1 means at least 10% of k-mers must support the classification.

---

## Step 4: Output Files

### Per-Read Classification (`*_kraken2.out`)

Tab-delimited with 5 columns:

| Column | Description |
|--------|-------------|
| 1 | C (classified) or U (unclassified) |
| 2 | Sequence ID |
| 3 | Taxonomy ID (0 if unclassified) |
| 4 | Sequence length |
| 5 | LCA mapping of each k-mer |

Example:
```
C	read_001	562	150|148	562:13 561:4 A:31 0:1 562:3
U	read_002	0	145|142	0:50 A:20
```

### Sample Report (`*_kraken2.report`)

Tab-delimited with 6-8 columns:

| Column | Description |
|--------|-------------|
| 1 | Percentage of reads in clade |
| 2 | Number of reads in clade |
| 3 | Number of reads assigned directly |
| 4 | Minimizers in clade (with --report-minimizer-data) |
| 5 | Distinct minimizers (with --report-minimizer-data) |
| 6 | Rank code (D, P, C, O, F, G, S) |
| 7 | NCBI Taxonomy ID |
| 8 | Scientific name |

Example:
```
 36.40	182	182	1688	18	S	211044	  Escherichia coli
```

---

## Step 5: Merge Reports and make CSV files

### Combine All Reports into a Single Table

```bash
cat > /projects/standard/gomeza/ssilvaju/Complex_Diet/scripts/kraken2_to_phyloseq.py << 'EOF'
#!/usr/bin/env python3
"""
Convert Kraken2 MPA-style reports to phyloseq-compatible CSV files.
Creates: otu_table.csv (counts) and tax_table.csv (taxonomy)
For reports generated with: --use-mpa-style --use-names
"""
import os
import sys
import re
from collections import defaultdict

# Regex to clean multiple spaces in taxon names
RE_SPACES = re.compile(r" +")

# MPA rank prefixes to full names
RANK_MAP = {
    'd': 'Domain',
    'k': 'Kingdom', 
    'p': 'Phylum',
    'c': 'Class',
    'o': 'Order',
    'f': 'Family',
    'g': 'Genus',
    's': 'Species'
}

RANKS = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

def parse_mpa_report(filepath):
    """Parse MPA-style Kraken2 report, extract species-level counts."""
    species_counts = {}
    
    with open(filepath) as f:
        for line in f:
            line = RE_SPACES.sub("_", line.strip())
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            
            taxon_path = parts[0]
            count = int(parts[1])
            
            # Only keep species level (contains 's__')
            if '|s__' in taxon_path:
                species_counts[taxon_path] = count
    
    return species_counts

def parse_taxonomy(taxon_path):
    """Parse MPA taxonomy string into rank dictionary."""
    taxonomy = {rank: 'NA' for rank in RANKS}
    
    levels = taxon_path.split('|')
    for level in levels:
        if '__' in level:
            prefix, name = level.split('__', 1)
            if prefix == 'd':
                taxonomy['Domain'] = name
            elif prefix == 'p':
                taxonomy['Phylum'] = name
            elif prefix == 'c':
                taxonomy['Class'] = name
            elif prefix == 'o':
                taxonomy['Order'] = name
            elif prefix == 'f':
                taxonomy['Family'] = name
            elif prefix == 'g':
                taxonomy['Genus'] = name
            elif prefix == 's':
                taxonomy['Species'] = name
    
    return taxonomy

def create_otu_id(taxon_path):
    """Create a clean OTU ID from taxonomy path."""
    # Extract species name for ID
    for level in taxon_path.split('|'):
        if level.startswith('s__'):
            species = level.replace('s__', '')
            # Clean special characters
            species = re.sub(r'[^a-zA-Z0-9_]', '_', species)
            return species
    return taxon_path.replace('|', '_')

def simplify_sample_name(filename):
    """Clean sample name from filename."""
    name = filename.replace("_kraken2.report", "")
    name = name.replace(".report", "")
    name = name.replace("_cleaned", "")
    return name

def main(report_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    # Collect all data
    all_taxa = set()
    sample_counts = {}
    
    report_files = sorted([f for f in os.listdir(report_dir) if f.endswith('.report')])
    
    print(f"Found {len(report_files)} report files")
    
    for filename in report_files:
        sample = simplify_sample_name(filename)
        filepath = os.path.join(report_dir, filename)
        
        counts = parse_mpa_report(filepath)
        sample_counts[sample] = counts
        all_taxa.update(counts.keys())
        
        print(f"  {sample}: {len(counts)} species, {sum(counts.values()):,} reads")
    
    # Sort taxa and samples
    taxa_list = sorted(all_taxa)
    samples = sorted(sample_counts.keys())
    
    print(f"\nTotal unique species: {len(taxa_list)}")
    print(f"Total samples: {len(samples)}")
    
    # Create OTU table
    print("\nWriting otu_table.csv...")
    with open(os.path.join(output_dir, 'otu_table.csv'), 'w') as f:
        # Header
        f.write("OTU_ID," + ",".join(samples) + "\n")
        
        # Data rows
        for taxon in taxa_list:
            otu_id = create_otu_id(taxon)
            counts = [str(sample_counts[s].get(taxon, 0)) for s in samples]
            f.write(f"{otu_id}," + ",".join(counts) + "\n")
    
    # Create taxonomy table
    print("Writing tax_table.csv...")
    with open(os.path.join(output_dir, 'tax_table.csv'), 'w') as f:
        # Header
        f.write("OTU_ID," + ",".join(RANKS) + "\n")
        
        # Data rows
        for taxon in taxa_list:
            otu_id = create_otu_id(taxon)
            taxonomy = parse_taxonomy(taxon)
            tax_values = [taxonomy[rank] for rank in RANKS]
            f.write(f"{otu_id}," + ",".join(tax_values) + "\n")
    
    # Create sample metadata template
    print("Writing sample_data.csv (template)...")
    with open(os.path.join(output_dir, 'sample_data.csv'), 'w') as f:
        f.write("SampleID,Group,Treatment,Timepoint\n")
        for sample in samples:
            f.write(f"{sample},,,\n")
    
    # Summary statistics
    print("\n" + "="*60)
    print("OUTPUT FILES CREATED:")
    print(f"  {output_dir}/otu_table.csv   - Count matrix ({len(taxa_list)} taxa x {len(samples)} samples)")
    print(f"  {output_dir}/tax_table.csv   - Taxonomy table ({len(taxa_list)} taxa x 7 ranks)")
    print(f"  {output_dir}/sample_data.csv - Sample metadata template (fill in manually)")
    print("="*60)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python kraken2_to_phyloseq.py <report_dir> <output_dir>")
        print("Example: python kraken2_to_phyloseq.py ./reports ./phyloseq_input")
        sys.exit(1)
    
    report_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    main(report_dir, output_dir)
EOF

chmod +x /projects/standard/gomeza/ssilvaju/Complex_Diet/scripts/kraken2_to_phyloseq.py
```

### Run Merge Script

```bash
python scripts/kraken2_to_phyloseq.py \
    data/kraken2/reports \
    data/kraken2/phyloseq_input
```

## Step 6: Quality Control

### Calculate Classification Rates

```bash
cat > ~/Complex_Diet/scripts/kraken2_summary.sh << 'EOF'
#!/bin/bash
# Summarize Kraken2 classification rates

REPORT_DIR=$1

echo -e "Sample\tTotal_Reads\tClassified\tUnclassified\tPct_Classified"

for report in ${REPORT_DIR}/*_kraken2.report; do
    sample=$(basename $report _kraken2.report)
    
    # Get unclassified count (first line, column 2)
    unclassified=$(head -1 $report | awk '{print $2}')
    
    # Get root count (classified reads)
    classified=$(grep -P "\tR\t" $report | head -1 | awk '{print $2}')
    
    # Calculate total and percentage
    total=$((classified + unclassified))
    if [ $total -gt 0 ]; then
        pct=$(echo "scale=2; $classified * 100 / $total" | bc)
    else
        pct=0
    fi
    
    echo -e "${sample}\t${total}\t${classified}\t${unclassified}\t${pct}%"
done
EOF

chmod +x ~/Complex_Diet/scripts/kraken2_summary.sh
```

### Run Summary

```bash
bash ~/Complex_Diet/scripts/kraken2_summary.sh \
    /projects/standard/gomeza/ssilvaju/Complex_Diet/data/kraken2/reports \
    > /projects/standard/gomeza/ssilvaju/Complex_Diet/data/kraken2/classification_summary.tsv
```

---

## Step 7: Filter Unidentified Reads (Following Paper Methods)

Based on the referenced paper methodology, filter out:
- Unclassified reads
- Reads only classified to Root
- Reads only classified to "Cellular organisms"

```bash
cat > ~/Complex_Diet/scripts/filter_kraken2.py << 'EOF'
#!/usr/bin/env python3
"""
Filter Kraken2 results following paper methodology:
- Remove unclassified reads
- Remove reads classified only to Root (taxid 1)
- Remove reads classified only to Cellular organisms (taxid 131567)
"""

import sys
import pandas as pd

def filter_kraken2_table(input_file, output_file):
    """Filter abundance table to remove uninformative classifications."""
    
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    
    # Taxa to remove (uninformative)
    remove_taxa = [
        'unclassified',
        'root',
        'cellular organisms',
        'Bacteria',  # Domain level only
        'Archaea',   # Domain level only
        'Viruses',   # Domain level only
    ]
    
    # Filter (case-insensitive)
    df_filtered = df[~df.index.str.lower().isin([t.lower() for t in remove_taxa])]
    
    # Remove taxa with zero counts across all samples
    df_filtered = df_filtered.loc[df_filtered.sum(axis=1) > 0]
    
    # Save
    df_filtered.to_csv(output_file, sep='\t')
    
    print(f"Original taxa: {len(df)}")
    print(f"Filtered taxa: {len(df_filtered)}")
    print(f"Removed: {len(df) - len(df_filtered)} taxa")

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    filter_kraken2_table(input_file, output_file)
EOF
```

---

## Step 8: Normalize to RPM

```bash
cat > ~/Complex_Diet/scripts/normalize_rpm.py << 'EOF'
#!/usr/bin/env python3
"""
Normalize Kraken2 counts to Reads Per Million (RPM).
"""

import sys
import pandas as pd

def normalize_to_rpm(input_file, output_file):
    """Normalize counts to RPM."""
    
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    
    # Calculate RPM: (count / total reads) * 1,000,000
    total_reads = df.sum(axis=0)
    df_rpm = df.div(total_reads, axis=1) * 1e6
    
    # Round to 2 decimal places
    df_rpm = df_rpm.round(2)
    
    # Save
    df_rpm.to_csv(output_file, sep='\t')
    
    print(f"Normalized {len(df.columns)} samples to RPM")
    print(f"Output: {output_file}")

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    normalize_to_rpm(input_file, output_file)
EOF
```

---
## Step 9: Create separate abundance tables for each domain at multiple taxonomic levels.

```

cat > /projects/standard/gomeza/ssilvaju/Complex_Diet/scripts/split_kraken2_tables.py << 'EOF'
#!/usr/bin/env python3
"""
Create separate abundance tables for each domain at multiple taxonomic levels.

Output structure:
    by_domain/
    ├── bacteria/
    │   ├── phylum.tsv
    │   ├── class.tsv
    │   ├── order.tsv
    │   ├── family.tsv
    │   ├── genus.tsv
    │   └── species.tsv
    ├── archaea/
    │   └── ...
    └── viruses/
        └── ...
"""

import os
import sys
import argparse
from collections import defaultdict

# Taxonomic level codes
LEVEL_NAMES = {
    'D': 'domain',
    'P': 'phylum',
    'C': 'class',
    'O': 'order',
    'F': 'family',
    'G': 'genus',
    'S': 'species'
}

def parse_report_by_domain(filepath):
    """
    Parse Kraken2 report and organize by domain and taxonomic level.
    
    Returns:
        dict: {domain: {level: {taxon: count}}}
    """
    data = defaultdict(lambda: defaultdict(dict))
    current_domain = None
    
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            
            # Parse columns
            if len(parts) == 8:
                pct, clade_reads, taxon_reads, mini, dist_mini, rank, taxid, name = parts
            else:
                pct, clade_reads, taxon_reads, rank, taxid, name = parts
            
            name = name.strip()
            rank = rank.strip()
            count = int(clade_reads)
            
            # Track current domain
            if rank == 'D':
                current_domain = name.lower().replace(' ', '_')
            
            # Get base level (S1, S2 -> S)
            base_rank = rank[0] if rank else ''
            
            # Store data
            if current_domain and base_rank in LEVEL_NAMES:
                level = LEVEL_NAMES[base_rank]
                data[current_domain][level][name] = count
    
    return data


def process_all_reports(report_dir, output_dir):
    """
    Process all reports and create domain/level-specific tables.
    """
    # Collect all data
    # Structure: {domain: {level: {taxon: {sample: count}}}}
    all_data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    samples = []
    
    report_files = sorted([f for f in os.listdir(report_dir) 
                          if f.endswith('.report')])
    
    print(f"Processing {len(report_files)} reports...")
    
    for filename in report_files:
        sample = filename.replace('_kraken2.report', '').replace('.report', '')
        samples.append(sample)
        filepath = os.path.join(report_dir, filename)
        
        sample_data = parse_report_by_domain(filepath)
        
        for domain, levels in sample_data.items():
            for level, taxa in levels.items():
                for taxon, count in taxa.items():
                    all_data[domain][level][taxon][sample] = count
    
    samples = sorted(set(samples))
    
    # Write output tables
    print(f"\nWriting output tables to: {output_dir}")
    
    for domain, levels in sorted(all_data.items()):
        domain_dir = os.path.join(output_dir, domain)
        os.makedirs(domain_dir, exist_ok=True)
        
        print(f"\n{domain.upper()}:")
        
        for level, taxa_dict in sorted(levels.items()):
            output_file = os.path.join(domain_dir, f"{level}.tsv")
            
            with open(output_file, 'w') as f:
                # Header
                f.write("Taxon\t" + "\t".join(samples) + "\n")
                
                # Data rows
                for taxon in sorted(taxa_dict.keys()):
                    counts = [str(taxa_dict[taxon].get(s, 0)) for s in samples]
                    f.write(f"{taxon}\t" + "\t".join(counts) + "\n")
            
            print(f"  {level}: {len(taxa_dict)} taxa")
    
    # Create summary
    summary_file = os.path.join(output_dir, "summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Domain\tLevel\tNumber_of_Taxa\n")
        for domain, levels in sorted(all_data.items()):
            for level, taxa_dict in sorted(levels.items()):
                f.write(f"{domain}\t{level}\t{len(taxa_dict)}\n")
    
    print(f"\nSummary written to: {summary_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Split Kraken2 results by domain and taxonomic level'
    )
    parser.add_argument('report_dir', help='Directory with Kraken2 reports')
    parser.add_argument('output_dir', help='Output directory')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.report_dir):
        print(f"ERROR: Directory not found: {args.report_dir}")
        sys.exit(1)
    
    process_all_reports(args.report_dir, args.output_dir)
EOF

chmod +x /projects/standard/gomeza/ssilvaju/Complex_Diet/scripts/split_kraken2_tables.py

```

Run it 

```
python /projects/standard/gomeza/ssilvaju/Complex_Diet/scripts/split_kraken2_tables.py \
    /projects/standard/gomeza/ssilvaju/Complex_Diet/data/kraken2/reports \
    /projects/standard/gomeza/ssilvaju/Complex_Diet/data/kraken2/by_domain
```

## Troubleshooting

### Database Loading Errors

```bash
# Use memory mapping if RAM is limited
kraken2 --db $DB --memory-mapping --paired ...
```

### Low Classification Rate

- Check database completeness
- Lower confidence threshold: `--confidence 0.05`
- Reduce minimum hit groups: `--minimum-hit-groups 1`

### Slow Performance

- Ensure database is on local/fast storage (not NFS)
- Increase threads: `--threads 16`
- Use SSD storage if available

---

## References

- [Kraken2 Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)
- [Kraken2 Manual](https://github.com/DerrickWood/kraken2/wiki/Manual)
- [Pre-built Databases](https://benlangmead.github.io/aws-indexes/k2)

---

## Citation

If using Kraken2, please cite:

```bibtex
@article{wood2019improved,
  title={Improved metagenomic analysis with Kraken 2},
  author={Wood, Derrick E and Lu, Jennifer and Langmead, Ben},
  journal={Genome biology},
  volume={20},
  number={1},
  pages={1--13},
  year={2019}
}
```
