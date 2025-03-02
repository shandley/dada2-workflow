#!/usr/bin/env Rscript

# Command-line script to run the DADA2 workflow with various options
# Supports single run and multi-run processing modes

# Load required libraries
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
library(argparse)

# Create argument parser
parser <- ArgumentParser(description="Run DADA2 workflow for 16S rRNA amplicon sequence processing")

# Add command line arguments
parser$add_argument("-m", "--multi-run", 
                   action="store_true", default=FALSE,
                   help="Enable multi-run processing mode")

parser$add_argument("-d", "--run-dir", type="character", default=NULL,
                   help="Directory containing run subdirectories (for multi-run mode)")

parser$add_argument("-b", "--big-data", 
                   action="store_true", default=FALSE,
                   help="Enable big data mode with optimized memory management")

parser$add_argument("-r", "--report", 
                   action="store_true", default=FALSE,
                   help="Generate HTML report")

parser$add_argument("-f", "--format", type="character", default="html_document",
                   help="Output format for report (e.g., html_document, pdf_document)")

parser$add_argument("-o", "--output-dir", type="character", default="reports",
                   help="Output directory for reports")

parser$add_argument("-n", "--output-file", type="character", default="dada2_workflow_report",
                   help="Base output filename for reports")

# Parse arguments
args <- parser$parse_args()

# Function to run the workflow
run_workflow <- function(args) {
  # Include timestamp in output
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Build parameter list for rendering
  params <- list(
    generate_report = args$report,
    output_format = args$format,
    output_dir = args$output_dir,
    output_file = paste0(args$output_file, "_", timestamp),
    multi_run = args$multi_run,
    run_dir = args$run_dir,
    big_data = args$big_data
  )
  
  # Display run configuration
  cat("=== DADA2 Workflow Configuration ===\n")
  cat("Timestamp:", timestamp, "\n")
  cat("Multi-run mode:", ifelse(args$multi_run, "ENABLED", "DISABLED"), "\n")
  if(args$multi_run) {
    cat("Run directory:", ifelse(is.null(args$run_dir), "data (default)", args$run_dir), "\n")
  }
  cat("Big data mode:", ifelse(args$big_data, "ENABLED", "DISABLED"), "\n")
  cat("Generate report:", ifelse(args$report, "YES", "NO"), "\n")
  if(args$report) {
    cat("Report format:", args$format, "\n")
    cat("Report output directory:", args$output_dir, "\n")
    cat("Report output file:", params$output_file, "\n")
  }
  cat("=====================================\n\n")
  
  # Check if rmarkdown is installed
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    install.packages("rmarkdown")
  }
  
  # Render the workflow
  cat("Starting DADA2 workflow...\n")
  rmarkdown::render(
    input = "dada2_workflow_optimize.Rmd",
    output_format = args$format,
    output_file = paste0(params$output_file, ".", ifelse(args$format == "html_document", "html", "pdf")),
    output_dir = args$output_dir,
    params = params,
    envir = new.env(parent = globalenv())
  )
}

# Check if output directory exists, create if needed
if (args$report && !dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
  cat("Created output directory:", args$output_dir, "\n")
}

# Run the workflow
run_workflow(args)