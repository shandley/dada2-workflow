# DADA2 Workflow for 16S rRNA Sequencing Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A generalized workflow for processing 16S rRNA gene amplicon data using the [DADA2](https://benjjneb.github.io/dada2/) package in R. This workflow identifies exact amplicon sequence variants (ASVs) with higher resolution than traditional OTU-based methods.

## Overview

This repository contains an RMarkdown workflow that implements the complete DADA2 pipeline:

1. **Quality filtering** - Remove low-quality reads
2. **Error rate learning** - Learn error rates from the data
3. **Dereplication** - Combine identical reads
4. **Sample inference** - Apply the DADA2 algorithm
5. **Merger of paired-end reads** - Align and merge paired reads
6. **Chimera removal** - Remove artifactual sequences
7. **Taxonomic assignment** - Classify ASVs
8. **Output generation** - Create tables and visualizations

## Requirements

- R ≥ 4.0.0
- Required R packages:
  - dada2
  - ggplot2
  - phyloseq
  - Biostrings
  - ShortRead
  - tidyverse

## Getting Started

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/dada2-workflow.git
   cd dada2-workflow
   ```

2. Install required R packages:
   ```r
   install.packages(c("ggplot2", "tidyverse"))
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install(c("dada2", "phyloseq", "Biostrings", "ShortRead"))
   ```

3. Reference database:
   - The workflow uses the Silva database (v138.1) for taxonomy assignments
   - The database is either accessed from the DADA2 package or downloaded from Zenodo
   - Both genus-level and species-level assignments are performed

4. Open `dada2_workflow.Rmd` in RStudio and update the following:
   - Path to your sequence files (supports both `_R1_001.fastq.gz` and `_R1.fastq.gz` naming conventions)
   - Read trimming parameters based on your data quality
   - Path to your sample metadata file

5. Execute the workflow by running all code chunks in the RMarkdown document.

6. View results in the interactive dashboard:
   ```r
   # Install required packages if needed
   install.packages(c("shiny", "shinydashboard", "plotly", "DT", "vegan", "viridis"))
   
   # Make sure phyloseq is installed
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install("phyloseq")
   
   # Run the dashboard (will open in your browser)
   source("run_dashboard.R")
   
   # If the dashboard doesn't open automatically, you can access it at:
   # http://127.0.0.1:4321
   
   # Note: If you encounter a greyed-out dashboard, try restarting R
   # and running the dashboard again
   ```

## Customizing the Workflow

The main parameters to adjust in the workflow:

- **Read trimming** (`truncLen`, `maxEE`): Modify based on your sequencing quality
- **Taxonomy thresholds**: Adjust confidence thresholds for taxonomic assignments
- **Visualization options**: Adjust plotting parameters as needed

## Output Files

The workflow produces several output files in the `results/` directory:

- `seqtab_nochim.csv`: ASV count table
- `taxonomy.csv`: Taxonomic assignments for each ASV (with species-level assignments)
- `phyloseq_object.rds`: R object for downstream analysis
- `ASVs.fasta`: FASTA file containing ASV sequences

## Interactive Dashboard

The included Shiny dashboard (`dashboard.R`) provides interactive visualization of DADA2 results:

- **Overview**: Quick summary statistics and sample metrics
- **Sample Quality**: Read count distributions and filtering statistics
- **Alpha Diversity**: Shannon, Simpson, Observed, and Chao1 diversity metrics
- **Beta Diversity**: PCoA, NMDS, t-SNE, and UMAP ordinations with PERMANOVA
- **Taxonomy**: Interactive bar plots and heatmaps of taxonomic composition
- **ASV Table**: Browse ASV sequences and abundance data
- **Differential Abundance**: Identify taxa that differ between sample groups using DESeq2 or ALDEx2

Launch the dashboard by running:
```r
source("run_dashboard.R")
```

## Contributing

Contributions to improve this workflow are welcome. Please feel free to submit a pull request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## References

- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016). "DADA2: High-resolution sample inference from Illumina amplicon data." Nature Methods, 13, 581-583. doi: 10.1038/nmeth.3869
- McMurdie PJ, Holmes S (2013). "phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data." PLoS ONE, 8(4):e61217
- Quast C, et al. (2013). "The SILVA ribosomal RNA gene database project: improved data processing and web-based tools." Nucleic Acids Research, 41(D1), D590-D596. doi: 10.1093/nar/gks1219
