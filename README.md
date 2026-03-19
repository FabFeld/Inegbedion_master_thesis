# Inegbedion Master Thesis - RNA-seq Analysis Pipelines

R analysis pipelines (DESeq2, ORA, GSEA) for RNA-seq data, used for the master thesis of Rachel Inegbedion at University of Bayreuth.

## Project Overview

Comprehensive R scripts for differential gene expression analysis and functional enrichment of RNA-seq data:
- **DESeq2** – Differential expression analysis
- **ORA** – Overrepresentation analysis
- **GSEA** – Gene Set Enrichment Analysis

---

## Requirements

- **R version**: 4.4.1 (validated)
- **Packages**:

| Package | Purpose |
|---------|---------|
| `DESeq2` | Differential expression analysis |
| `clusterProfiler` | ORA and GSEA |
| `tidyverse` | Data wrangling and visualization |
| `org.Dmagna.eg.db` | *Daphnia magna* annotations |
| `org.Dr.eg.db` | *Danio rerio* annotations |

Additional dependencies are listed in individual scripts.

---

## Architecture
```
├── DESeq2_pairwise_dmagna.R
├── DESeq2_pairwise_drerio.R
├── ORA_GSEA_for_DESeq2_output.R
├── ORA+GSEA_for_orgDb.Dmagna_SAC.R
└── fun/                 # Reusable helper functions
```

---

## Raw Data

Raw count matrices and metadata for this analysis are available on **Zenodo**:

🔗 **[Zenodo Repository](10.5281/zenodo.19112223)**

**Format**:
- Count matrices: TSV/CSV (genes × samples)
- Metadata (coldata): Sample conditions, replicates, treatment groups

---

## Usage

### 1. Clone the repository

```bash
git clone https://github.com/FabFeld/Inegbedion_master_thesis.git
cd Inegbedion_master_thesis
```

### 2. Setup

1. Ensure `fun/` is in your working directory (or adjust the path in scripts)
2. Download and extract raw data from Zenodo
3. Set working directory to the target dataset:

```r
### Example: zebrafish / sulfamic acid
setwd("Drerio/Sulfamic_acid")
```

### 3. Run analysis

```r
source("DESeq2_pairwise_drerio.R")
```

> **Note**: Each dataset must be analyzed separately.

---


## Script Descriptions

| Script | Purpose | Organism | Output |
|--------|---------|----------|--------|
| `DESeq2_pairwise_dmagna.R` | Differential expression analysis | *Daphnia magna* | DE genes, normalized counts, plots |
| `DESeq2_pairwise_drerio.R` | Differential expression analysis | *Danio rerio* | DE genes, normalized counts, plots |
| `ORA_GSEA_for_DESeq2_output.R` | Functional enrichment on DE results | *Danio rerio*| GO terms, pathway enrichment |
| `ORA+GSEA_for_orgDb.Dmagna_SAC.R` | Custom organism DB enrichment | *D. magna* | Enriched pathways, visualization |


