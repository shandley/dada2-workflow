#!/usr/bin/env Rscript

# Create demo data for DADA2 multi-run dashboard testing
# This script creates sample data with multiple "runs" for testing the batch effects visualization

# Libraries
library(phyloseq)
library(dada2)

# Create sample data for testing multi-run functionality
create_multi_run_demo_data <- function() {
  # Create directory for sample data
  if (!dir.exists("results")) dir.create("results")
  
  # Create parameters for 3 "runs" with slightly different characteristics
  n_samples_per_run <- 10
  n_taxa <- 200
  n_runs <- 3
  
  # Run names
  run_names <- paste0("run", 1:n_runs)
  
  # Create sample names for each run
  sample_names <- c()
  for (run in run_names) {
    sample_names <- c(sample_names, paste0(run, "_sample", 1:n_samples_per_run))
  }
  
  # Parameters for the different runs
  # Each run will have slightly different abundance distributions to simulate batch effects
  run_means <- list(
    run1 = c(0.3, 0.2, 0.15, 0.1, 0.05), # Run 1 has higher abundance of first taxa
    run2 = c(0.1, 0.3, 0.25, 0.2, 0.05), # Run 2 has higher abundance of middle taxa
    run3 = c(0.05, 0.1, 0.15, 0.2, 0.4)  # Run 3 has higher abundance of last taxa
  )
  
  # Create different abundance distributions for each run
  otu_mat <- matrix(0, nrow = n_taxa, ncol = length(sample_names))
  rownames(otu_mat) <- paste0("ASV", 1:n_taxa)
  colnames(otu_mat) <- sample_names
  
  # Generate different abundance profiles per run
  set.seed(123) # For reproducibility
  
  sample_idx <- 1
  for (run_idx in 1:n_runs) {
    run_name <- run_names[run_idx]
    
    # Select which samples belong to this run
    run_samples <- grep(paste0("^", run_name), sample_names)
    
    # Create taxa with different abundance patterns per run
    # First 5 taxa will have specific abundance patterns per run
    for (i in 1:5) {
      otu_mat[i, run_samples] <- rnorm(length(run_samples), 
                                     mean = run_means[[run_name]][i] * 10000, 
                                     sd = 1000)
    }
    
    # Remaining taxa have more random distribution
    for (i in 6:n_taxa) {
      # Use exponential distribution for long tail of rare taxa
      if (i <= 50) {
        # More abundant taxa
        otu_mat[i, run_samples] <- rexp(length(run_samples), rate = 1/1000) +
          rnorm(length(run_samples), mean = 0, sd = 100)
      } else {
        # Less abundant taxa
        otu_mat[i, run_samples] <- rexp(length(run_samples), rate = 1/100) +
          rnorm(length(run_samples), mean = 0, sd = 20)
      }
    }
  }
  
  # Ensure all values are non-negative and round to integers
  otu_mat[otu_mat < 0] <- 0
  otu_mat <- round(otu_mat)
  
  # Create taxonomy table
  tax_mat <- matrix("", nrow = n_taxa, ncol = 7)
  colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rownames(tax_mat) <- rownames(otu_mat)
  
  # Fill taxonomy - create a realistic bacterial community
  # Phyla with different distributions across runs
  phyla <- c("Bacteroidetes", "Firmicutes", "Proteobacteria", "Actinobacteria", "Verrucomicrobia")
  phyla_probs <- list(
    run1 = c(0.4, 0.3, 0.15, 0.1, 0.05),  # Run 1: more Bacteroidetes
    run2 = c(0.2, 0.4, 0.2, 0.15, 0.05),  # Run 2: more Firmicutes
    run3 = c(0.1, 0.2, 0.4, 0.2, 0.1)     # Run 3: more Proteobacteria
  )
  
  # Classes per phylum
  classes <- list(
    Bacteroidetes = c("Bacteroidia", "Cytophagia", "Flavobacteriia", "Sphingobacteriia"),
    Firmicutes = c("Bacilli", "Clostridia", "Erysipelotrichia", "Negativicutes"),
    Proteobacteria = c("Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria"),
    Actinobacteria = c("Actinobacteria", "Coriobacteriia", "Thermoleophilia"),
    Verrucomicrobia = c("Verrucomicrobiae", "Opitutae", "Spartobacteria")
  )
  
  # Orders and families - simplified for demo
  orders_per_class <- 2
  families_per_order <- 3
  genera_per_family <- 2
  
  # Assign taxonomy
  for (i in 1:n_taxa) {
    asv <- rownames(otu_mat)[i]
    
    # Find which run this ASV is more prevalent in
    asv_counts <- otu_mat[i,]
    run_totals <- sapply(run_names, function(run) {
      sum(asv_counts[grep(paste0("^", run), sample_names)])
    })
    dominant_run <- run_names[which.max(run_totals)]
    
    # Pick phylum based on run-specific probabilities
    tax_mat[i, "Kingdom"] <- "Bacteria"
    phylum <- sample(phyla, 1, prob = phyla_probs[[dominant_run]])
    tax_mat[i, "Phylum"] <- phylum
    
    # Pick class
    class_options <- classes[[phylum]]
    tax_mat[i, "Class"] <- sample(class_options, 1)
    
    # Generate order, family, genus
    tax_mat[i, "Order"] <- paste0(substring(tax_mat[i, "Class"], 1, 4), "ales_", sample(1:orders_per_class, 1))
    tax_mat[i, "Family"] <- paste0(substring(tax_mat[i, "Order"], 1, 4), "aceae_", sample(1:families_per_order, 1))
    tax_mat[i, "Genus"] <- paste0(substring(tax_mat[i, "Family"], 1, 4), "_", sample(1:genera_per_family, 1))
    
    # Species - just use ASV ID 
    tax_mat[i, "Species"] <- paste0("sp_", substring(asv, 4))
  }
  
  # Create sample metadata with run information
  sample_data_df <- data.frame(
    Sample_ID = sample_names,
    Run = gsub("_sample.*", "", sample_names),
    Group = sample(c("Control", "Treatment"), length(sample_names), replace = TRUE),
    Batch = rep(1:3, each = n_samples_per_run * n_runs/3),
    ReadCount = colSums(otu_mat)
  )
  rownames(sample_data_df) <- sample_names
  
  # Create phyloseq object
  otu <- otu_table(otu_mat, taxa_are_rows = TRUE)
  tax <- tax_table(tax_mat)
  samples <- sample_data(sample_data_df)
  
  ps <- phyloseq(otu, tax, samples)
  
  # Save the phyloseq object
  saveRDS(ps, "results/phyloseq_object.rds")
  
  # Generate tracking data
  tracking_df <- data.frame(
    Sample = sample_names,
    Input = round(colSums(otu_mat) * runif(length(sample_names), 1.3, 1.5)),
    Filtered = round(colSums(otu_mat) * runif(length(sample_names), 1.1, 1.3)),
    Denoised = round(colSums(otu_mat) * runif(length(sample_names), 1.0, 1.1)),
    Merged = round(colSums(otu_mat) * runif(length(sample_names), 0.9, 1.0)),
    NonChimeric = colSums(otu_mat),
    Run = gsub("_sample.*", "", sample_names)
  )
  
  saveRDS(tracking_df, "results/read_tracking.rds")
  
  cat("Multi-run demo data created successfully!\n")
  cat("- Created", n_runs, "simulated runs with", n_samples_per_run, "samples each\n")
  cat("- Total of", n_taxa, "ASVs with different distributions across runs\n")
  cat("- Saved as results/phyloseq_object.rds and results/read_tracking.rds\n")
}

# If this script is run directly (not sourced), create the demo data
if (!interactive()) {
  create_multi_run_demo_data()
}