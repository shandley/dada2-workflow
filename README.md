# DADA2 Optimized Workflow for 16S rRNA Sequencing Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An optimized workflow for processing 16S rRNA gene amplicon data using the [DADA2](https://benjjneb.github.io/dada2/) package in R. This workflow identifies exact amplicon sequence variants (ASVs) with higher resolution than traditional OTU-based methods while implementing advanced optimizations for improved performance, reliability, and insights.

## Overview

This repository contains two main components:

1. **dada2_workflow_optimize.Rmd** - An enhanced, performance-optimized RMarkdown workflow that implements the complete DADA2 pipeline
2. **dashboard.R** - An interactive Shiny dashboard for visualizing and exploring results

*Note: A basic DADA2 workflow (dada2_workflow.Rmd) is also included for users who want the simplest possible implementation.*

## Key Features of the Optimized Workflow

- **Automatic sequencing platform detection** - Identifies platform based on read length and quality patterns
- **Parameter optimization** - Tunes filtering parameters to your specific sequence data characteristics
- **Adaptive truncation lengths** - Dynamically adjusts to optimize read quality and overlap
- **Expected error threshold optimization** - Balances quality control and read retention
- **Primer detection** - Automatically identifies primer sequences for amplicon size calculations
- **Memory-optimized processing** - Efficient batched operations with automatic memory management
- **Enhanced checkpointing system** - Robust recovery from interruptions with comprehensive tracking
- **Parallelized execution** - Automatic multi-core utilization with adaptive worker allocation
- **Reference-based taxonomy confidence scoring** - Bootstrap confidence values for all taxonomic assignments
- **Multi-method taxonomy assignment** - Combines results from multiple classifiers for improved accuracy
- **Phylogenetic tree construction** - Integrates phylogenetic information using optimized alignment methods
- **Rarefaction analysis** - Depth optimization with saturation detection for proper diversity comparisons
- **Detailed quality visualization** - Enhanced plots with quality interpretation zones
- **Multi-run support** - Process and integrate data from multiple sequencing runs with batch effect analysis
- **Comprehensive reporting** - Generate detailed HTML/PDF reports with code hiding option

## Dashboard Features

The interactive dashboard (`dashboard.R`) provides advanced visualization and analysis capabilities:

- **Overview Panel** - Sample metrics, read statistics, and ASV summaries
- **Quality Control** - Filtering performance, read tracking, and quality metric distributions
- **Alpha Diversity** - Multiple indices with statistical comparisons between groups
- **Beta Diversity** - Multiple ordination methods (PCoA, NMDS, t-SNE, UMAP) with statistical tests
- **Taxonomy Explorer** - Interactive hierarchical visualization of taxonomic composition
- **ASV Browser** - Searchable ASV table with sequence information and abundance patterns
- **Differential Abundance** - Multiple testing methods for identifying biomarkers between groups
- **Batch Effect Analysis** - For multi-run studies, quantifies and visualizes run effects
- **Normalization Methods** - Compare various count normalization approaches
- **Export Options** - Download plots, tables, and processed data in multiple formats

## Requirements

- R ≥ 4.0.0
- Required R packages:
  - dada2
  - ggplot2
  - phyloseq
  - Biostrings
  - ShortRead
  - tidyverse
  - future (for parallelization)
  - DECIPHER (for improved taxonomy & phylogeny)
  - vegan (for diversity analyses)

## Getting Started

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/dada2-workflow.git
   cd dada2-workflow
   ```

2. Install required R packages:
   ```r
   install.packages(c("ggplot2", "tidyverse", "argparse", "future", "future.apply"))
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install(c("dada2", "phyloseq", "Biostrings", "ShortRead", "DECIPHER"))
   ```

3. Run the optimized workflow using one of the following methods:

   **Method 1: RStudio**
   - Open `dada2_workflow_optimize.Rmd` in RStudio 
   - Update parameters in the YAML header if needed
   - Execute by running all code chunks

   **Method 2: Command Line**
   - For a single sequencing run:
     ```bash
     Rscript run_dada2_workflow.R
     ```
   
   - For multi-run analysis with batch effect correction:
     ```bash
     Rscript run_dada2_workflow.R --multi-run --run-dir path/to/run_directory
     ```

4. View results in the interactive dashboard:
   ```bash
   Rscript run_dashboard.R
   ```
   Or for advanced options:
   ```bash
   Rscript run_dashboard.R --optimize --cores 4 --multi-run
   ```

## Directory Structure

### Single Run Mode
Place your fastq files in the `data/` directory:
```
data/
  ├── sample1_R1.fastq.gz
  ├── sample1_R2.fastq.gz
  ├── sample2_R1.fastq.gz
  ├── sample2_R2.fastq.gz
  └── ...
```

### Multi-Run Mode
Organize your data with each run in a separate subdirectory:
```
data/
  ├── run1/
  │   ├── sample1_R1.fastq.gz
  │   ├── sample1_R2.fastq.gz
  │   └── ...
  ├── run2/
  │   ├── sampleA_R1.fastq.gz
  │   ├── sampleA_R2.fastq.gz
  │   └── ...
  └── run3/
      ├── sampleX_R1.fastq.gz
      ├── sampleX_R2.fastq.gz
      └── ...
```

## Command-Line Options

The `run_dada2_workflow.R` script provides the following options:

```
usage: run_dada2_workflow.R [-h] [-m] [-d RUN_DIR] [-b] [-r] [-f FORMAT]
                            [-o OUTPUT_DIR] [-n OUTPUT_FILE] [--cores CORES]
                            [--optimize]

Run optimized DADA2 workflow for 16S rRNA amplicon sequence processing

optional arguments:
  -h, --help            show this help message and exit
  -m, --multi-run       Enable multi-run processing mode
  -d RUN_DIR, --run-dir RUN_DIR
                        Directory containing run subdirectories (for multi-run
                        mode)
  -b, --big-data        Enable big data mode with optimized memory management
  -r, --report          Generate HTML report
  -f FORMAT, --format FORMAT
                        Output format for report (e.g., html_document,
                        pdf_document)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory for reports
  -n OUTPUT_FILE, --output-file OUTPUT_FILE
                        Base output filename for reports
  --cores CORES         Number of CPU cores to use for parallelization
  --optimize            Enable additional performance optimizations
```

## Multi-Run Processing and Batch Effect Analysis

For studies with samples across multiple sequencing runs, the workflow:

1. Processes each run separately through the sample inference step with run-specific error models
2. Merges sequence tables from all runs while preserving run information
3. Performs chimera removal and taxonomic assignment on the combined data
4. Provides batch effect detection and correction methods:
   - PERMANOVA to test for significant run effects
   - Beta dispersion analysis to check homogeneity across runs
   - Batch effect visualization with ordination methods
   - Optional normalization methods specifically for batch correction

This approach gives you:
- More accurate error models specific to each sequencing run
- Detection of potential batch effects that could bias results
- Methods to correct or account for batch effects in downstream analyses
- Better integration of data from different sequencing platforms or centers

## Output Files

The workflow produces comprehensive output files in the `results/` directory:

- `seqtab_nochim.csv`: ASV count table
- `taxonomy.csv`: Taxonomic assignments for each ASV
- `taxonomy_with_confidence.csv`: Taxonomy with bootstrap confidence scores
- `phyloseq_object.rds`: R object for downstream analysis
- `ASVs.fasta`: FASTA file containing ASV sequences
- `filter_summary.csv`: Quality filtering statistics
- `chimera_summary.csv`: Chimera detection statistics
- `read_tracking_detailed.csv`: Read counts through each pipeline step
- `rarefaction_curves.rds`: Data for rarefaction analysis
- `workflow_summary.rds`: Complete statistics about the analysis run

## Contributing

Contributions to improve this workflow are welcome. Please feel free to submit a pull request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## References

- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016). "DADA2: High-resolution sample inference from Illumina amplicon data." Nature Methods, 13, 581-583. doi: 10.1038/nmeth.3869
- McMurdie PJ, Holmes S (2013). "phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data." PLoS ONE, 8(4):e61217
- Murali A, Bhargava A, Wright ES (2018). "IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences." Microbiome, 6, 140. doi: 10.1186/s40168-018-0521-5
