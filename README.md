# Inegbedion Master Thesis - RNA-seq Analysis Pipelines

R analysis pipelines (DESeq2, ORA, GSEA) for RNA-seq data, for the master thesis of Rachel Inegbedion at University of Bayreuth.

## Project Overview

This repository contains comprehensive R scripts for differential gene expression analysis and functional enrichment analysis of RNA-seq data. The project includes differential expression analysis using DESeq2, as well as overrepresentation analysis (ORA) and gene set enrichment analysis (GSEA) for biological interpretation.

---

## Requirements

### System Requirements
- **R version**: validated only for 4.4.1

### R Packages
The scripts depend on the following R packages:
- `DESeq2` - Differential expression analysis
- `dplyr` - Data manipulation
- `tidyverse` - Data wrangling and visualization
- `clusterProfiler` - Overrepresentation analysis (ORA) and GSEA
- `org.Dmagna.eg.db` - Daphnia magna gene annotations
- `org.Dr.eg.db` - Zebrafish gene annotations
- `ggplot2` - Data visualization
- Additional dependencies as specified in individual scripts


---

## Architecture

The project follows a modular pipeline structure organized into three main analysis stages:

### 1. **Differential Expression Analysis (DESeq2)**
   - `DESeq2_pairwise_dmagna.R` - DESeq2 pairwise comparisons for *Daphnia magna*
   - `DESeq2_pairwise_drerio.R` - DESeq2 pairwise comparisons for *Danio rerio* (zebrafish)
   
Identifies significantly differentially expressed genes between experimental conditions

### 2. **Functional Enrichment Analysis**
   - `ORA_GSEA_for_DESeq2_output.R` - Overrepresentation analysis and GSEA on DESeq2 results
   - `ORA+GSEA_for_orgDb.Dmagna_SAC.R` - ORA and GSEA specifically for *D. magna* using custom organism database
   
Identifies enriched biological pathways and gene ontology terms among significantly expressed genes

### 3. **Helper Functions (fun/ folder)**
   - Contains reusable R functions for common operations
   - Functions are sourced by the main analysis scripts
   - Must be available in the working directory or referenced via proper file paths

**Data Flow**:
```
Raw RNA-seq Data (Count matrices)
        ↓
  [DESeq2 Analysis]
        ↓
 Normalized counts & DE genes
        ↓
  [ORA/GSEA Analysis]
        ↓
 Enriched biological pathways
```

---

## Raw Data Location

Raw count matrices and metadata for this analysis are available on **Zenodo**:

🔗 **[Zenodo Repository](https://zenodo.org/)**

### Data Requirements
- **Count matrix format**: Tab-delimited or CSV files with gene IDs as rows and samples as columns
- **Metadata file**: Sample information including experimental conditions, replicates, and treatment groups
- Place raw data files in your working directory or specify the full path in scripts

---

## Using the Scripts

### Setup Instructions

1. **Clone the repository**:
   ```bash
   git clone https://github.com/FabFeld/Inegbedion_master_thesis.git
   cd Inegbedion_master_thesis
   ```

2. **Prepare your working directory**:
   ```bash
   # Ensure the fun/ folder is in your working directory
   ls fun/  # Verify helper functions are present
   ```

3. **Load the helper functions** (required):
   At the beginning of each analysis script, the `fun/` folder functions are sourced. Ensure this path is correct:
   ```r
   # In the main scripts:
   source("fun/function_name.R")  # Adjust path if needed
   ```

4. **Modify input/output paths**:
   - Update file paths in scripts to point to your data location
   - Specify output directories for results

5. **Run the analysis**:
   ```r
   # For DESeq2 analysis
   source("DESeq2_pairwise_dmagna.R")
   # or
   source("DESeq2_pairwise_drerio.R")
   ```

### Important Notes
- The `fun/` folder **must be available** in your working directory or the location specified in the scripts
- Update data paths and parameters in each script before running
- Ensure all required R packages are installed before execution

---

## Script Descriptions

| Script | Purpose | Organism | Output |
|--------|---------|----------|--------|
| `DESeq2_pairwise_dmagna.R` | Differential expression analysis | *Daphnia magna* | DE genes, normalized counts, plots |
| `DESeq2_pairwise_drerio.R` | Differential expression analysis | *Danio rerio* | DE genes, normalized counts, plots |
| `ORA_GSEA_for_DESeq2_output.R` | Functional enrichment on DE results | *Danio rerio*| GO terms, pathway enrichment |
| `ORA+GSEA_for_orgDb.Dmagna_SAC.R` | Custom organism DB enrichment | *D. magna* | Enriched pathways, visualization |


