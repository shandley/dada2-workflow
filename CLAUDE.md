# CLAUDE.md for dada2-workflow

## Build/Run Commands
- Run R script: `Rscript script_name.R`
- Interactive R session: `R` or `Rscript -e "source('script_name.R')"`
- Install packages: `install.packages("package_name")`
- Load DADA2: `library(dada2)`
- Run optimized workflow: `Rscript -e "rmarkdown::render('dada2_workflow_optimize.Rmd')"`
- Run dashboard with batch effect features: `Rscript run_dashboard.R --multi-run`
- Run dashboard with advanced performance options: `Rscript run_dashboard.R --optimize --cores 4`
- Launch workflow app: `Rscript launch_workflow_app.R`

## Optimization Features
- Automatic parallelization based on system capabilities
- Batched processing for memory-intensive operations
- Chunked processing for large datasets
- Efficient garbage collection for long-running operations
- Multi-run batch effect detection and correction

## Parallelization Settings
- Use `future` package for parallel processing
- Automatically detect cores: `parallel::detectCores()`
- Set parallelization plan: `plan(multisession, workers = optimal_workers)`
- Memory management: `options(future.globals.maxSize = 2 * 1024^3)`
- Batch size recommendations:
  - Filtering: 20 samples/batch
  - Error learning: 40 samples max
  - Sample inference: 10 samples/batch
  - Merging: 10 samples/batch
  - Chimera removal: 2000 ASVs/batch
  - Taxonomy assignment: 1000 ASVs/batch

## Batch Effect Analysis & Normalization Methods
- PERMANOVA (adonis2): Tests for differences in microbial composition between runs
- Beta dispersion: Tests for homogeneity of multivariate dispersions across runs
- RDA (Redundancy Analysis): Assesses variance explained by run factor
- Mantel test: Correlates run membership with community dissimilarities
- ANOSIM: Non-parametric test similar to PERMANOVA

- Normalization methods available in dashboard:
  - Relative abundance: Simple proportion scaling
  - CSS: Cumulative sum scaling from metagenomeSeq
  - RLE: Relative log expression from DESeq2
  - TMM: Trimmed mean of M-values from edgeR
  - CLR: Centered log-ratio transformation (compositional data analysis)
  - ALR: Additive log-ratio transformation (compositional data analysis)
  - VST: Variance stabilizing transformation from DESeq2

## Code Style Guidelines
- Function names: Use snake_case (e.g., `filter_and_trim`)
- Variable names: Use descriptive snake_case (e.g., `filtered_reads`)
- Indentation: 2 spaces (no tabs)
- Line length: Max 80 characters
- Documentation: Use roxygen2-style comments for functions
- Error handling: Use tryCatch() for graceful error handling
- Use pipes (%>%) for readable data transformations
- Separate concerns into distinct functions
- Comment complex algorithms or parameter choices

## Repo Structure
- R scripts in root directory
- Data in separate data/ directory (not committed to git)
- Results in separate results/ directory (not committed to git)
- Checkpoints in separate checkpoints/ directory (not committed to git)

## CRAN Mirror Management
- Automatic CRAN mirror detection and selection
- Regional mirror selection based on timezone (US, Europe, Asia, Australia)
- Fallback to reliable global CDN mirrors
- Connection testing with 5-second timeout
- Handles mirror errors gracefully with multiple fallback options