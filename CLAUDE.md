# CLAUDE.md for dada2-workflow

## Build/Run Commands
- Run R script: `Rscript script_name.R`
- Interactive R session: `R` or `Rscript -e "source('script_name.R')"`
- Install packages: `install.packages("package_name")`
- Load DADA2: `library(dada2)`

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