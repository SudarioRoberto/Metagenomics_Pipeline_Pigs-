# 05 - Statistical Analysis

This folder contains the complete statistical analysis workflow for taxonomic abundance data generated from the Kraken2/Bracken pipeline.

## 🔗 Interactive Report

**[View the Full Interactive Analysis Report](https://sudarioroberto.github.io/Metagenomics_Pipeline_Pigs_Complex_Diets/)**

---

## Overview

The analysis follows best practices for microbiome data, including compositional-aware methods and multiple testing correction as recommended by recent guidelines.

### Analyses Performed

| Analysis | Methods | Description |
|----------|---------|-------------|
| **Alpha Diversity** | Observed, Shannon, Simpson, Chao1, Pielou | Within-sample diversity metrics |
| **Beta Diversity** | Bray-Curtis, Jaccard, Aitchison | Between-sample dissimilarity |
| **Ordination** | PCoA | Visualize community structure |
| **PERMANOVA** | adonis2 | Test for treatment effects |
| **Dispersion** | betadisper | Test homogeneity of dispersions |
| **Differential Abundance** | ALDEx2, ANCOM-BC | Identify treatment-associated taxa |
| **Visualization** | ggplot2, pheatmap | Publication-ready figures |

---

## Directory Contents

```
05_Statistical_Analysis/
│
├── 📄 README.md                      # This file
├── 📄 complex_diets.qmd              # Quarto source file (main analysis)
├── 📄 complex_diets.html             # Rendered HTML report
│
├── 📄 meta.csv                       # Sample metadata
│
├── 📊 Abundance Tables (Bracken output)
│   ├── merged_D_abundance.txt        # Domain level
│   ├── merged_P_abundance.txt        # Phylum level
│   ├── merged_C_abundance.txt        # Class level
│   ├── merged_O_abundance.txt        # Order level
│   ├── merged_F_abundance.txt        # Family level
│   ├── merged_G_abundance.txt        # Genus level
│   └── merged_S_abundance.txt        # Species level
│
└── 📂 figures/                       # Exported figures (optional)
```

---

## Experimental Design

| Factor | Levels | Description |
|--------|--------|-------------|
| **Complex** | low, medium, high | Diet complexity level |
| **Protease** | cont, prot | Protease supplementation (control vs treatment) |
| **Interaction** | 6 groups | Full factorial design |

---

## Statistical Methods

### Alpha Diversity

- **Metrics**: Observed richness, Shannon, Simpson, Inverse Simpson, Chao1, Pielou's evenness
- **Tests**: Kruskal-Wallis (3+ groups), Wilcoxon (2 groups)
- **Correction**: Benjamini-Hochberg FDR

### Beta Diversity

- **Distances**: 
  - Bray-Curtis (abundance-weighted)
  - Jaccard (presence/absence)
  - Aitchison (compositional, CLR-transformed)
- **Ordination**: Principal Coordinates Analysis (PCoA)
- **Tests**: PERMANOVA with 999 permutations
- **Pairwise**: pairwise.adonis with BH correction

### Differential Abundance

Following the consensus approach recommended by Nearing et al. (2022):

1. **ALDEx2**: CLR transformation with Monte Carlo sampling
   - Welch's t-test and Wilcoxon for 2 groups
   - Kruskal-Wallis for 3+ groups
   - Effect size estimation

2. **ANCOM-BC**: Bias-corrected log-linear model
   - Structural zero detection
   - Sample-specific normalization

3. **Consensus**: Taxa significant in both methods

---

## Running the Analysis

### Prerequisites

```r
# Required R packages
install.packages(c("tidyverse", "vegan", "ggpubr", "patchwork", "pheatmap"))

# Bioconductor packages
BiocManager::install(c("phyloseq", "ALDEx2", "ANCOMBC", "mia"))

# Additional packages
install.packages("pairwiseAdonis", repos = "https://github.com/pmartinezarbizu/pairwiseAdonis")
```

### Render the Report

```bash
# Using Quarto CLI
quarto render complex_diets.qmd

# Or in RStudio: Click "Render" button
```

---

## Input Data Format

### Abundance Tables

Tab-delimited files with taxa as rows and samples as columns. Each sample has two columns: `*_num` (read counts) and `*_frac` (relative abundance).

```
taxon                     sample1.bracken_num  sample1.bracken_frac  ...
Lactobacillus johnsonii   15234                0.0523                ...
Prevotella copri          8921                 0.0306                ...
```

### Metadata (meta.csv)

```csv
,complex,protease,interaction
2827,low,cont,low_cont
2828,low,prot,low_prot
...
```

---

## Key Outputs

### Figures Generated

| Figure | Description |
|--------|-------------|
| Alpha diversity boxplots | By Complex, Protease, and Interaction |
| PCoA ordinations | Bray-Curtis and Aitchison distances |
| Volcano plots | ALDEx2 and ANCOM-BC results |
| Heatmaps | Significant taxa with interactions |
| Barplots | Taxonomic composition at multiple levels |

### Statistical Tables

| Table | Description |
|-------|-------------|
| PERMANOVA summary | R² and p-values for all distance metrics |
| Pairwise comparisons | Post-hoc tests between groups |
| Differential abundance | Significant taxa from ALDEx2 and ANCOM-BC |
| Consensus taxa | Taxa significant in both methods |

---

## Important Notes

### FDR Correction

Results are reported both with raw p-values and BH-adjusted p-values. When adjusted p-values are not significant but raw p-values are, this is noted for biological interpretation.

### Compositional Data

Microbiome data is compositional (relative abundances sum to 1). This analysis uses:
- CLR transformation for ordination
- ALDEx2 (compositional-aware)
- ANCOM-BC (bias-corrected)

### Zero Handling

- Pseudocounts added for CLR transformation
- Prevalence filtering (taxa in ≥10% of samples)
- Zero-variance taxa removed for interaction analysis

---

## References

- Nearing JT et al. (2022). Microbiome differential abundance methods produce different results across 38 datasets. *Nature Communications*.
- Gloor GB et al. (2017). Microbiome datasets are compositional: and this is not optional. *Frontiers in Microbiology*.
- Lin H & Peddada SD (2020). Analysis of compositions of microbiomes with bias correction. *Nature Communications*.

---

## Contact

For questions about this analysis, please open an issue on the GitHub repository.

