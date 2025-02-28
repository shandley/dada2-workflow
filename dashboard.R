## DADA2 Results Dashboard
## Interactive visualization of DADA2 workflow outputs

# Load required libraries
library(shiny)           # Web application framework for R
library(shinydashboard)  # Dashboard interface for Shiny apps
library(phyloseq)        # Handling and analysis of microbiome census data
library(plotly)          # Interactive web-based data visualization
library(DT)              # R interface to the DataTables JavaScript library
library(tidyverse)       # Collection of packages for data manipulation and visualization
library(vegan)           # Community ecology package for diversity analysis
library(ggplot2)         # Data visualization package for creating graphics
library(viridis)         # Colorblind-friendly color palettes for visualizations

# UI
ui <- dashboardPage(
  dashboardHeader(title = "DADA2 Results Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Sample Quality", tabName = "quality", icon = icon("chart-line")),
      menuItem("Alpha Diversity", tabName = "alpha", icon = icon("chart-bar")),
      menuItem("Beta Diversity", tabName = "beta", icon = icon("project-diagram")),
      menuItem("Taxonomy", tabName = "taxonomy", icon = icon("sitemap")),
      menuItem("ASV Table", tabName = "asv", icon = icon("table"))
    ),
    
    hr(),
    
    # Filtering controls
    conditionalPanel(
      condition = "input.tab == 'alpha' || input.tab == 'beta' || input.tab == 'taxonomy' || input.tab == 'asv'",
      sliderInput("minReads", "Minimum Read Depth:", 
                  min = 1, max = 10000, value = 1000),
      sliderInput("prevalence", "ASV Prevalence (%):", 
                  min = 0, max = 100, value = 5)
    )
  ),
  
  dashboardBody(
    tabItems(
      # Overview tab
      tabItem(tabName = "overview",
              fluidRow(
                box(width = 12, title = "DADA2 Pipeline Summary", status = "primary",
                    "This dashboard provides interactive visualizations of the DADA2 analysis results."),
                
                valueBoxOutput("totalSamples", width = 3),
                valueBoxOutput("totalASVs", width = 3),
                valueBoxOutput("totalReads", width = 3),
                valueBoxOutput("medianReads", width = 3)
              ),
              
              fluidRow(
                box(width = 6, title = "Sample Read Counts", status = "info",
                    plotlyOutput("sampleReads")),
                box(width = 6, title = "Read Processing Stats", status = "info",
                    plotlyOutput("readStats"))
              ),
              
              fluidRow(
                box(width = 12, title = "Taxonomic Composition Overview", status = "info",
                    plotlyOutput("taxaOverview"))
              )
      ),
      
      # Sample Quality tab
      tabItem(tabName = "quality",
              fluidRow(
                box(width = 12, title = "Sample Read Quality", status = "info",
                    "This tab shows the quality metrics for samples after DADA2 processing.")
              ),
              
              # Quality Profiles section
              fluidRow(
                box(width = 12, title = "Quality Profiles", status = "info",
                    "Examine the aggregated quality scores of reads to assess data quality and determine optimal trimming parameters.")
              ),
              
              fluidRow(
                box(width = 6, title = "Forward Read Quality Profiles", status = "info",
                    plotOutput("forwardQualityProfilePlot", height = "500px")),
                box(width = 6, title = "Reverse Read Quality Profiles", status = "info",
                    plotOutput("reverseQualityProfilePlot", height = "500px"))
              ),
              
              fluidRow(
                box(width = 12, status = "info",
                    p("Aggregated quality score distribution across the first 10 samples. The blue line shows the mean quality score at each position, 
                      and the gray ribbon represents the range of quality scores observed across samples."))
              ),
              
              # Original sections
              fluidRow(
                box(width = 6, title = "Read Count Distribution", status = "info",
                    plotlyOutput("readCountHist")),
                box(width = 6, title = "Tracking Reads Through Pipeline", status = "info",
                    plotlyOutput("readTracking"))
              ),
              
              fluidRow(
                box(width = 12, title = "Sample Read Counts", status = "info",
                    DT::dataTableOutput("sampleTable"))
              )
      ),
      
      # Alpha Diversity tab
      tabItem(tabName = "alpha",
              fluidRow(
                box(width = 12, title = "Alpha Diversity Metrics", status = "info",
                    "This tab shows different alpha diversity metrics for samples.")
              ),
              
              fluidRow(
                box(width = 4, 
                    selectInput("alphaMethod", "Diversity Metric:", 
                                choices = c("Shannon", "Simpson", "Observed", "Chao1"),
                                selected = "Shannon")),
                box(width = 8, 
                    checkboxInput("logScale", "Log Scale", value = FALSE))
              ),
              
              fluidRow(
                box(width = 12, title = "Alpha Diversity by Group", status = "info",
                    plotlyOutput("alphaPlot"))
              )
      ),
      
      # Beta Diversity tab
      tabItem(tabName = "beta",
              fluidRow(
                box(width = 12, title = "Beta Diversity Analysis", status = "info",
                    "This tab shows beta diversity metrics and ordination plots.")
              ),
              
              fluidRow(
                box(width = 4,
                    selectInput("betaMethod", "Distance Method:",
                                choices = c("Bray-Curtis", "Jaccard", "UniFrac", "Weighted UniFrac"),
                                selected = "Bray-Curtis")),
                box(width = 4,
                    selectInput("ordMethod", "Ordination Method:",
                                choices = c("PCoA", "NMDS", "t-SNE", "UMAP"),
                                selected = "PCoA"))
              ),
              
              fluidRow(
                box(width = 12, title = "Ordination Plot", status = "info",
                    plotlyOutput("betaPlot"))
              ),
              
              fluidRow(
                box(width = 12, title = "Sample Distance Heatmap", status = "info",
                    plotlyOutput("distHeatmap"))
              )
      ),
      
      # Taxonomy tab
      tabItem(tabName = "taxonomy",
              fluidRow(
                box(width = 12, title = "Taxonomic Composition", status = "info",
                    "This tab shows the taxonomic composition of samples.")
              ),
              
              fluidRow(
                box(width = 3,
                    selectInput("taxLevel", "Taxonomic Level:", 
                                choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                selected = "Phylum")),
                box(width = 3,
                    selectInput("plotType", "Plot Type:",
                                choices = c("Bar Plot", "Heatmap"),
                                selected = "Bar Plot")),
                box(width = 3,
                    numericInput("topTaxa", "Top N Taxa:", value = 10, min = 1, max = 100)),
                box(width = 3,
                    checkboxInput("otherCat", "Group Remaining as 'Other'", value = TRUE))
              ),
              
              fluidRow(
                box(width = 12, title = "Taxonomic Composition Plot", status = "info",
                    plotlyOutput("taxaPlot"))
              ),
              
              fluidRow(
                box(width = 12, title = "Taxonomy Table", status = "info",
                    DT::dataTableOutput("taxaTable"))
              )
      ),
      
      # ASV Table tab
      tabItem(tabName = "asv",
              fluidRow(
                box(width = 12, title = "ASV Table", status = "info",
                    "This tab shows the ASV counts and sequence information.")
              ),
              
              fluidRow(
                box(width = 4,
                    selectInput("taxonomyDisplayLevel", "Show taxonomy at level:", 
                                choices = c("All Levels", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                selected = "All Levels"))
              ),
              
              fluidRow(
                box(width = 12, title = "ASV Data Table", status = "info",
                    DT::dataTableOutput("asvTable"))
              ),
              
              fluidRow(
                box(width = 12, title = "ASV Sequence Information", status = "info",
                    verbatimTextOutput("seqInfo"))
              )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Load phyloseq object
  phyloseq_obj <- reactive({
    # Path to saved phyloseq object
    ps_path <- "results/phyloseq_object.rds"
    
    # Check if file exists
    if (!file.exists(ps_path)) {
      return(NULL)
    }
    
    # Try to safely load and process the phyloseq object
    tryCatch({
      # Load the phyloseq object
      ps <- readRDS(ps_path)
      
      # Safe check if taxa table exists
      has_tax_table <- tryCatch({
        !is.null(phyloseq::tax_table(ps))
      }, error = function(e) {
        FALSE
      })
      
      # Safe check for taxa orientation
      if (has_tax_table) {
        # Get orientation info safely
        tax_rows <- tryCatch({
          is.matrix(phyloseq::tax_table(ps)) && phyloseq::taxa_are_rows(phyloseq::tax_table(ps))
        }, error = function(e) {
          NULL
        })
        
        otu_rows <- tryCatch({
          phyloseq::taxa_are_rows(phyloseq::otu_table(ps))
        }, error = function(e) {
          NULL
        })
        
        # Only proceed with orientation fix if we have valid information
        if (!is.null(tax_rows) && !is.null(otu_rows) && tax_rows != otu_rows) {
          warning("Fixing orientation mismatch between OTU and tax tables")
          # Transpose tax table if needed
          phyloseq::tax_table(ps) <- t(phyloseq::tax_table(ps))
        }
      }
      
      return(ps)
    }, error = function(e) {
      # If there's an error, try to create a basic phyloseq object
      warning("Error processing phyloseq object: ", conditionMessage(e))
      warning("Creating simplified phyloseq object")
      
      tryCatch({
        # Try to extract and reconstruct just the OTU table
        otu_data <- as.matrix(readRDS(ps_path))
        if (is.matrix(otu_data)) {
          # Create a simple phyloseq object with just the OTU table
          simple_otu <- phyloseq::otu_table(otu_data, taxa_are_rows = FALSE)
          return(phyloseq::phyloseq(simple_otu))
        } else {
          NULL
        }
      }, error = function(e2) {
        warning("Could not create simplified phyloseq object: ", conditionMessage(e2))
        return(NULL)
      })
    })
  })
  
  # Load tracking data
  tracking_data <- reactive({
    # Path to read tracking file (we'll assume it's saved separately)
    track_path <- "results/read_tracking.rds"
    
    # If file doesn't exist, return NULL
    if (!file.exists(track_path)) {
      return(NULL)
    }
    
    # Load tracking data with error handling
    tryCatch({
      track_data <- readRDS(track_path)
      
      # Make sure it's properly formatted with required column
      if (!("Sample" %in% colnames(track_data))) {
        # If no Sample column, try to add it
        if (is.data.frame(track_data) && nrow(track_data) > 0) {
          track_data$Sample <- paste0("Sample", 1:nrow(track_data))
        }
      }
      
      return(track_data)
    }, error = function(e) {
      warning("Error loading tracking data: ", conditionMessage(e))
      return(NULL)
    })
  })
  
  # No need to update UI elements based on sample data
  
  # Filter phyloseq object based on user inputs
  filtered_ps <- reactive({
    ps <- phyloseq_obj()
    if (is.null(ps)) return(NULL)
    
    # Extra validation to ensure phyloseq object is valid
    valid_phyloseq <- tryCatch({
      # Basic validation: check if it has an OTU table
      has_otu <- !is.null(phyloseq::otu_table(ps))
      
      # Return result
      has_otu
    }, error = function(e) FALSE)
    
    if (!valid_phyloseq) {
      warning("Invalid phyloseq object detected. Creating a basic one.")
      return(NULL)
    }
    
    # Try all filtering operations with error handling
    tryCatch({
      # Filter by read depth
      if (!is.null(input$minReads)) {
        tryCatch({
          sample_sums_val <- phyloseq::sample_sums(ps)
          ps <- phyloseq::prune_samples(sample_sums_val >= input$minReads, ps)
        }, error = function(e) {
          warning("Error in prune_samples: ", conditionMessage(e))
        })
      }
      
      # Filter by prevalence - with careful validation
      if (!is.null(input$prevalence)) {
        tryCatch({
          # Get OTU table with correct orientation
          otu <- phyloseq::otu_table(ps)
          
          # Check orientation to calculate prevalence correctly
          if (phyloseq::taxa_are_rows(otu)) {
            prevalence <- apply(otu, 1, function(x) sum(x > 0) / length(x) * 100)
          } else {
            prevalence <- apply(otu, 2, function(x) sum(x > 0) / length(x) * 100)
          }
          
          # Make sure prevalence vector has names matching taxa
          taxa_names <- phyloseq::taxa_names(ps)
          if (length(prevalence) == length(taxa_names)) {
            names(prevalence) <- taxa_names
            # Use prune_taxa with pre-calculated logical vector
            to_keep <- prevalence >= input$prevalence
            ps <- phyloseq::prune_taxa(to_keep, ps)
          }
        }, error = function(e) {
          warning("Error in prevalence filtering: ", conditionMessage(e))
        })
      }
      
      return(ps)
    }, error = function(e) {
      warning("Error in filtered_ps: ", conditionMessage(e))
      return(ps)  # Return original if filtering fails
    })
  })
  
  # Overview tab outputs
  output$totalSamples <- renderValueBox({
    ps <- filtered_ps()
    if (is.null(ps)) return(valueBox(0, "Samples", icon = icon("vials"), color = "blue"))
    
    # Safely get sample count
    n_samples <- tryCatch({
      phyloseq::nsamples(ps)
    }, error = function(e) {
      # Fallback if nsamples fails
      length(phyloseq::sample_names(ps))
    })
    
    valueBox(n_samples, "Samples", icon = icon("vials"), color = "blue")
  })
  
  output$totalASVs <- renderValueBox({
    ps <- filtered_ps()
    if (is.null(ps)) return(valueBox(0, "ASVs", icon = icon("dna"), color = "green"))
    
    # Safely get taxa count
    n_taxa <- tryCatch({
      phyloseq::ntaxa(ps)
    }, error = function(e) {
      # Fallback if ntaxa fails
      length(phyloseq::taxa_names(ps))
    })
    
    valueBox(n_taxa, "ASVs", icon = icon("dna"), color = "green")
  })
  
  output$totalReads <- renderValueBox({
    ps <- filtered_ps()
    if (is.null(ps)) return(valueBox(0, "Total Reads", icon = icon("list"), color = "yellow"))
    
    # Safely get total reads
    total_reads <- tryCatch({
      sum(phyloseq::sample_sums(ps))
    }, error = function(e) {
      # Fallback to OTU table sum if sample_sums fails
      sum(phyloseq::otu_table(ps))
    })
    
    valueBox(format(total_reads, big.mark = ","), 
             "Total Reads", icon = icon("list"), color = "yellow")
  })
  
  output$medianReads <- renderValueBox({
    ps <- filtered_ps()
    if (is.null(ps)) return(valueBox(0, "Median Reads/Sample", icon = icon("tachometer-alt"), color = "red"))
    
    # Safely get median reads
    median_reads <- tryCatch({
      median(phyloseq::sample_sums(ps))
    }, error = function(e) {
      # Fallback to mean if median fails
      mean(rowSums(phyloseq::otu_table(ps)))
    })
    
    valueBox(format(median_reads, big.mark = ","), 
             "Median Reads/Sample", icon = icon("tachometer-alt"), color = "red")
  })
  
  output$sampleReads <- renderPlotly({
    ps <- filtered_ps()
    if (is.null(ps)) return(NULL)
    
    # Create data frame with sample sums
    df <- data.frame(
      Sample = sample_names(ps),
      Reads = sample_sums(ps)
    )
    
    # Create plot
    p <- ggplot(df, aes(x = reorder(Sample, -Reads), y = Reads)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Sample", y = "Read Count") +
      scale_y_continuous(labels = scales::comma)
    
    ggplotly(p)
  })
  
  output$readStats <- renderPlotly({
    track <- tracking_data()
    if (is.null(track)) return(NULL)
    
    # Make sure track is a data frame with required columns
    if (!is.data.frame(track)) {
      return(NULL)
    }
    
    # Make sure the Sample column exists
    if (!"Sample" %in% colnames(track)) {
      # Try to add Sample column
      track$Sample <- paste0("Sample", 1:nrow(track))
    }
    
    # Try to reshape data safely
    tryCatch({
      # Keep only numeric columns besides Sample
      numeric_cols <- sapply(track, is.numeric)
      id_vars <- c("Sample")
      measure_vars <- setdiff(names(track)[numeric_cols], id_vars)
      
      # Only proceed if we have measure variables
      if (length(measure_vars) == 0) {
        return(NULL)
      }
      
      # Create melted data safely with base R
      track_melt <- data.frame(
        Sample = rep(track$Sample, length(measure_vars)),
        Step = rep(measure_vars, each = nrow(track)),
        Reads = unlist(lapply(measure_vars, function(var) track[[var]]))
      )
      
      # Create plot
      p <- ggplot(track_melt, aes(x = Step, y = Reads, group = Sample, color = Sample)) +
        geom_line() +
        geom_point() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Read Count") +
        scale_y_continuous(labels = scales::comma)
      
      ggplotly(p)
    }, error = function(e) {
      # Return a blank plot with error message
      warning("Error in read stats plot: ", conditionMessage(e))
      plot_ly() %>% 
        layout(title = "Error: Could not process tracking data",
               annotations = list(
                 x = 0.5, y = 0.5, 
                 text = paste("Error:", conditionMessage(e)), 
                 showarrow = FALSE
               ))
    })
  })
  
  output$taxaOverview <- renderPlotly({
    ps <- filtered_ps()
    if (is.null(ps)) return(NULL)
    
    # Safely get taxonomy data
    has_tax_table <- tryCatch({
      !is.null(phyloseq::tax_table(ps))
    }, error = function(e) {
      FALSE
    })
    
    # If we have a taxonomy table with Phylum information
    if (has_tax_table && "Phylum" %in% colnames(phyloseq::tax_table(ps))) {
      # Try to agglomerate at phylum level
      tryCatch({
        ps_phylum <- phyloseq::tax_glom(ps, taxrank = "Phylum")
        
        # Transform to relative abundance
        ps_phylum_rel <- phyloseq::transform_sample_counts(ps_phylum, function(x) x / sum(x))
        
        # Create data frame for plotting
        phylum_data <- phyloseq::psmelt(ps_phylum_rel)
        
        # Handle NA values in Phylum
        phylum_data$Phylum <- as.character(phylum_data$Phylum)
        phylum_data$Phylum[is.na(phylum_data$Phylum)] <- "Unknown"
        
        # Get taxa sums safely
        taxa_sums_safe <- tryCatch({
          phyloseq::taxa_sums(ps_phylum_rel)
        }, error = function(e) {
          # Create a named vector with OTU abundances
          otu_abundances <- rowSums(phyloseq::otu_table(ps_phylum_rel))
          names(otu_abundances) <- phyloseq::taxa_names(ps_phylum_rel)
          otu_abundances
        })
        
        # Get top phyla (up to 10)
        n_phyla <- min(10, length(unique(phylum_data$Phylum)))
        top_phyla_names <- names(sort(taxa_sums_safe, decreasing = TRUE)[1:n_phyla])
        
        # If we have valid top phyla, use them
        if (length(top_phyla_names) > 0) {
          phylum_data$Phylum <- ifelse(phylum_data$OTU %in% top_phyla_names, 
                                    phylum_data$Phylum, "Other")
        }
        
        # Create plot
        p <- ggplot(phylum_data, aes(x = Sample, y = Abundance, fill = Phylum)) +
          geom_bar(stat = "identity") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Sample", y = "Relative Abundance") +
          scale_fill_viridis(discrete = TRUE)
        
        return(ggplotly(p))
      }, error = function(e) {
        # Fallback to simple abundance plot if agglomeration fails
        warning("Error in taxonomy overview: ", conditionMessage(e))
        
        # Get sample data
        sample_counts <- phyloseq::sample_sums(ps)
        sample_data <- data.frame(
          Sample = names(sample_counts),
          Reads = sample_counts
        )
        
        # Create basic read count plot
        p <- ggplot(sample_data, aes(x = Sample, y = Reads)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Sample", y = "Read Count") +
          ggtitle("Sample Read Counts (Taxonomy overview failed)")
        
        return(ggplotly(p))
      })
    } else {
      # Create simple abundance plot if no taxonomy
      # Get sample data
      sample_counts <- phyloseq::sample_sums(ps)
      sample_data <- data.frame(
        Sample = names(sample_counts),
        Reads = sample_counts
      )
      
      # Create basic read count plot
      p <- ggplot(sample_data, aes(x = Sample, y = Reads)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Sample", y = "Read Count") +
        ggtitle("Sample Read Counts (No taxonomy data available)")
      
      ggplotly(p)
    }
  })
  
  # Load raw fastq paths
  raw_fastq_paths <- reactive({
    # Look for fastq files in the data directory
    fastq_dir <- "data"
    
    # Check if directory exists
    if (!dir.exists(fastq_dir)) {
      return(NULL)
    }
    
    # Try to find fastq files
    tryCatch({
      # Look for fastq files (support multiple extensions)
      fastq_files <- list.files(path = fastq_dir, 
                               pattern = "\\.(fastq|fastq\\.gz|fq|fq\\.gz)$", 
                               full.names = TRUE)
      
      if (length(fastq_files) == 0) {
        return(NULL)
      }
      
      # Separate forward and reverse reads
      forward_pattern <- "_R1_|_1\\.|_F\\."
      reverse_pattern <- "_R2_|_2\\.|_R\\."
      
      forward_files <- fastq_files[grepl(forward_pattern, fastq_files)]
      reverse_files <- fastq_files[grepl(reverse_pattern, fastq_files)]
      
      # Only proceed if we found both forward and reverse reads
      if (length(forward_files) == 0 && length(reverse_files) == 0) {
        # Try to be more flexible in patterns
        forward_files <- fastq_files[grepl("1|forward|f", tolower(basename(fastq_files)))]
        reverse_files <- fastq_files[grepl("2|reverse|r", tolower(basename(fastq_files)))]
      }
      
      # Create a named list to store the paths
      list(forward = forward_files, reverse = reverse_files)
    }, error = function(e) {
      warning("Error loading fastq files: ", conditionMessage(e))
      return(NULL)
    })
  })

  # No need for sample selection update, since we're always aggregating

  # Render forward read quality profile plot
  output$forwardQualityProfilePlot <- renderPlot({
    # Get fastq paths
    fastq_paths <- raw_fastq_paths()
    
    # Check if we have data
    if (is.null(fastq_paths) || length(fastq_paths[["forward"]]) == 0) {
      # Try to create a dummy quality plot for forward reads
      return(create_dummy_quality_plot("forward"))
    }
    
    # Use all available forward read files
    read_files <- fastq_paths[["forward"]]
    
    # Limit to first 10 files for aggregated plot
    if (length(read_files) > 10) {
      read_files <- read_files[1:10]
    }
    
    # Try to generate quality profiles
    tryCatch({
      # Load the dada2 package
      if (!requireNamespace("dada2", quietly = TRUE)) {
        stop("DADA2 package is required for quality profiles")
      }
      
      # Generate and return the plot with aggregate=TRUE to show aggregated profiles
      dada2::plotQualityProfile(read_files, aggregate = TRUE)
    }, error = function(e) {
      warning("Error generating forward quality profile: ", conditionMessage(e))
      create_dummy_quality_plot("forward")
    })
  })
  
  # Render reverse read quality profile plot
  output$reverseQualityProfilePlot <- renderPlot({
    # Get fastq paths
    fastq_paths <- raw_fastq_paths()
    
    # Check if we have data
    if (is.null(fastq_paths) || length(fastq_paths[["reverse"]]) == 0) {
      # Try to create a dummy quality plot for reverse reads
      return(create_dummy_quality_plot("reverse"))
    }
    
    # Use all available reverse read files
    read_files <- fastq_paths[["reverse"]]
    
    # Limit to first 10 files for aggregated plot
    if (length(read_files) > 10) {
      read_files <- read_files[1:10]
    }
    
    # Try to generate quality profiles
    tryCatch({
      # Load the dada2 package
      if (!requireNamespace("dada2", quietly = TRUE)) {
        stop("DADA2 package is required for quality profiles")
      }
      
      # Generate and return the plot with aggregate=TRUE to show aggregated profiles
      dada2::plotQualityProfile(read_files, aggregate = TRUE)
    }, error = function(e) {
      warning("Error generating reverse quality profile: ", conditionMessage(e))
      create_dummy_quality_plot("reverse")
    })
  })

  # Helper function to create a dummy aggregated quality profile plot 
  # This is used when real fastq files aren't available
  create_dummy_quality_plot <- function(read_type) {
    # Create dummy data
    set.seed(123)
    cycle <- 1:100
    mean_quality <- c(35 + rnorm(30, 0, 1), 
                      35 + rnorm(30, 0, 2) - seq(0, 6, length.out = 30), 
                      28 + rnorm(40, 0, 3) - seq(6, 15, length.out = 40))
    
    # Different pattern for reverse reads
    if (read_type == "reverse") {
      mean_quality <- c(30 + rnorm(20, 0, 2), 
                        30 + rnorm(40, 0, 2) - seq(0, 8, length.out = 40), 
                        22 + rnorm(40, 0, 3) - seq(8, 15, length.out = 40))
    }
    
    # Create data for error bars (to simulate aggregated plot)
    lower_q <- mean_quality - runif(length(mean_quality), 1, 3)
    upper_q <- mean_quality + runif(length(mean_quality), 1, 3)
    
    # Create a data frame for the plot
    qual_data <- data.frame(
      cycle = cycle,
      quality = mean_quality,
      lower = lower_q,
      upper = upper_q
    )
    
    # Create a quality profile plot similar to dada2's aggregated plotQualityProfile
    p <- ggplot(qual_data, aes(x = cycle, y = quality)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "darkgrey") +
      geom_line(size = 1, color = "blue") +
      labs(x = "Cycle", y = "Quality Score", 
           title = paste("Aggregated Quality Profile -", ifelse(read_type == "forward", "Forward", "Reverse"), "Reads"),
           subtitle = "No actual fastq files found. This is a simulated aggregated profile for demonstration.") +
      theme_minimal() +
      scale_y_continuous(limits = c(0, 40)) +
      annotate("rect", xmin = 0, xmax = 100, ymin = 0, ymax = 20, alpha = 0.2, fill = "red") +
      annotate("rect", xmin = 0, xmax = 100, ymin = 20, ymax = 30, alpha = 0.2, fill = "yellow") +
      annotate("rect", xmin = 0, xmax = 100, ymin = 30, ymax = 40, alpha = 0.2, fill = "green")
    
    return(p)
  }

  # Sample Quality tab outputs
  output$readCountHist <- renderPlotly({
    ps <- filtered_ps()
    if (is.null(ps)) return(NULL)
    
    # Create data frame with read counts
    read_counts <- data.frame(
      Reads = sample_sums(ps)
    )
    
    # Create histogram
    p <- ggplot(read_counts, aes(x = Reads)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "black") +
      theme_minimal() +
      labs(x = "Read Count", y = "Number of Samples") +
      scale_x_continuous(labels = scales::comma)
    
    ggplotly(p)
  })
  
  output$readTracking <- renderPlotly({
    track <- tracking_data()
    if (is.null(track)) return(NULL)
    
    # Extra validation to handle potential errors
    tryCatch({
      # Make sure track is a data frame
      if (!is.data.frame(track)) {
        stop("Tracking data is not a data frame")
      }
      
      # Make sure there are numeric columns to plot
      numeric_cols <- sapply(track, is.numeric)
      if (sum(numeric_cols) == 0) {
        stop("No numeric columns found in tracking data")
      }
      
      # Identify the Sample column (if any)
      sample_col <- NULL
      for (col_name in colnames(track)) {
        if (grepl("sample", tolower(col_name))) {
          sample_col <- col_name
          break
        }
      }
      
      # If we found a sample column, exclude it from calculations
      if (!is.null(sample_col)) {
        # Keep only numeric columns excluding the Sample column
        track_numeric <- track[, setdiff(names(which(numeric_cols)), sample_col), drop = FALSE]
      } else {
        # Just use all numeric columns
        track_numeric <- track[, numeric_cols, drop = FALSE]
      }
      
      # Convert any character columns to numeric
      track_numeric <- as.data.frame(sapply(track_numeric, function(x) {
        if (is.character(x)) {
          as.numeric(x)
        } else {
          x
        }
      }))
      
      # Calculate mean reads per step, handling NAs
      track_mean <- colMeans(track_numeric, na.rm = TRUE)
      
      # Create data frame for plotting
      track_df <- data.frame(
        Step = names(track_mean),
        Reads = as.numeric(track_mean)
      )
      
      # Extra check to make sure reads are numeric
      if (!is.numeric(track_df$Reads)) {
        track_df$Reads <- as.numeric(track_df$Reads)
      }
      
      # Remove any NaN or Inf values
      track_df <- track_df[is.finite(track_df$Reads), ]
      
      # If we have data to plot
      if (nrow(track_df) > 0) {
        # Create plot
        p <- ggplot(track_df, aes(x = Step, y = Reads, group = 1)) +
          geom_line(size = 1.2, color = "steelblue") +
          geom_point(size = 3, color = "steelblue") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(y = "Mean Read Count") +
          scale_y_continuous(labels = scales::comma)
        
        return(ggplotly(p))
      } else {
        # Create empty plot with message
        return(plot_ly() %>% 
               layout(title = "No valid tracking data to display",
                     annotations = list(
                       x = 0.5, y = 0.5, 
                       text = "No valid tracking data available", 
                       showarrow = FALSE
                     )))
      }
    }, error = function(e) {
      # Create empty plot with error message
      return(plot_ly() %>% 
             layout(title = "Error in read tracking plot",
                   annotations = list(
                     x = 0.5, y = 0.5, 
                     text = paste("Error:", conditionMessage(e)), 
                     showarrow = FALSE
                   )))
    })
  })
  
  output$sampleTable <- renderDataTable({
    ps <- filtered_ps()
    if (is.null(ps)) return(NULL)
    
    # Create data frame with just sample names and read counts
    sample_df <- data.frame(
      Sample = sample_names(ps),
      ReadCount = sample_sums(ps),
      row.names = NULL
    )
    
    # Format read count column
    sample_df$ReadCount <- format(sample_df$ReadCount, big.mark = ",")
    
    DT::datatable(sample_df, options = list(pageLength = 10))
  })
  
  # Alpha Diversity tab outputs
  output$alphaPlot <- renderPlotly({
    ps <- filtered_ps()
    if (is.null(ps) || is.null(input$alphaMethod)) return(NULL)
    
    # Calculate alpha diversity
    alpha_div <- estimate_richness(ps, measures = input$alphaMethod)
    
    # Create data frame for plotting
    alpha_df <- data.frame(
      Sample = sample_names(ps),
      Diversity = alpha_div[[input$alphaMethod]]
    )
    
    # Create plot - simple bar plot since we don't have group information
    p <- ggplot(alpha_df, aes(x = reorder(Sample, -Diversity), y = Diversity)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Sample", y = paste(input$alphaMethod, "Diversity"))
    
    # Apply log scale if selected
    if (input$logScale) {
      p <- p + scale_y_log10()
    }
    
    ggplotly(p)
  })
  
  # Beta Diversity tab outputs
  output$betaPlot <- renderPlotly({
    ps <- filtered_ps()
    if (is.null(ps) || is.null(input$betaMethod) || is.null(input$ordMethod)) return(NULL)
    
    # Transform to relative abundance
    ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
    
    # Calculate distance matrix
    dist_method <- switch(input$betaMethod,
                          "Bray-Curtis" = "bray",
                          "Jaccard" = "jaccard",
                          "UniFrac" = "unifrac",
                          "Weighted UniFrac" = "wunifrac")
    
    # For UniFrac methods, check if tree is available
    if (dist_method %in% c("unifrac", "wunifrac") && is.null(phy_tree(ps_rel))) {
      return(NULL)
    }
    
    dist_matrix <- phyloseq::distance(ps_rel, method = dist_method)
    
    # Perform ordination
    ord_method <- switch(input$ordMethod,
                         "PCoA" = "PCoA",
                         "NMDS" = "NMDS",
                         "t-SNE" = "tsne",
                         "UMAP" = "umap")
    
    ord <- ordinate(ps_rel, method = ord_method, distance = dist_matrix)
    
    # Get ordination coordinates
    if (ord_method == "PCoA") {
      ord_data <- data.frame(ord$vectors[,1:2])
      colnames(ord_data) <- c("Axis1", "Axis2")
      var_explained <- round(ord$values$Relative_eig[1:2] * 100, 1)
      axis_labels <- paste0("Axis ", 1:2, " (", var_explained, "%)")
    } else if (ord_method == "NMDS") {
      ord_data <- data.frame(scores(ord)[,1:2])
      colnames(ord_data) <- c("Axis1", "Axis2")
      axis_labels <- c("NMDS1", "NMDS2")
    } else {
      ord_data <- as.data.frame(ord$points)
      colnames(ord_data) <- c("Axis1", "Axis2")
      axis_labels <- c("Dim1", "Dim2")
    }
    
    # Add sample data
    ord_data$Sample <- rownames(ord_data)
    
    # Create plot
    p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, text = Sample)) +
      geom_point(size = 3, alpha = 0.7, color = "steelblue") +
      theme_minimal() +
      labs(x = axis_labels[1], y = axis_labels[2])
    
    ggplotly(p, tooltip = "text")
  })
  
  output$distHeatmap <- renderPlotly({
    ps <- filtered_ps()
    if (is.null(ps) || is.null(input$betaMethod)) return(NULL)
    
    # Transform to relative abundance
    ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
    
    # Calculate distance matrix
    dist_method <- switch(input$betaMethod,
                          "Bray-Curtis" = "bray",
                          "Jaccard" = "jaccard",
                          "UniFrac" = "unifrac",
                          "Weighted UniFrac" = "wunifrac")
    
    # For UniFrac methods, check if tree is available
    if (dist_method %in% c("unifrac", "wunifrac") && is.null(phy_tree(ps_rel))) {
      return(NULL)
    }
    
    dist_matrix <- phyloseq::distance(ps_rel, method = dist_method)
    
    # Convert distance matrix to regular matrix
    dist_mat <- as.matrix(dist_matrix)
    
    # Create heatmap
    p <- plot_ly(
      x = rownames(dist_mat),
      y = rownames(dist_mat),
      z = dist_mat,
      type = "heatmap",
      colorscale = "Viridis"
    ) %>%
      layout(
        title = paste(input$betaMethod, "Distance Heatmap"),
        xaxis = list(title = ""),
        yaxis = list(title = "")
      )
    
    p
  })
  
  # Taxonomy tab outputs
  output$taxaPlot <- renderPlotly({
    # Add tryCatch to handle all errors
    tryCatch({
      ps <- filtered_ps()
      if (is.null(ps) || is.null(input$taxLevel)) {
        return(plot_ly() %>% 
               layout(title = "No data available",
                     annotations = list(
                       x = 0.5, y = 0.5, 
                       text = "No valid data for plotting", 
                       showarrow = FALSE
                     )))
      }
      
      # Check for valid tax_table - critical for this function
      has_tax_table <- tryCatch({
        tt <- phyloseq::tax_table(ps)
        !is.null(tt) && ncol(tt) > 0 && nrow(tt) > 0
      }, error = function(e) FALSE)
      
      if (!has_tax_table) {
        return(plot_ly() %>% 
               layout(title = "No taxonomy data available",
                     annotations = list(
                       x = 0.5, y = 0.5, 
                       text = "Phyloseq object missing taxonomy table", 
                       showarrow = FALSE
                     )))
      }
      
      # Check if requested taxonomy level exists
      if (!input$taxLevel %in% colnames(phyloseq::tax_table(ps))) {
        return(plot_ly() %>% 
               layout(title = paste("No", input$taxLevel, "data available"),
                     annotations = list(
                       x = 0.5, y = 0.5, 
                       text = paste("Taxonomy table does not contain", input$taxLevel), 
                       showarrow = FALSE
                     )))
      }
      
      # Safely agglomerate at selected taxonomic level
      ps_safe <- ps
      
      # Try to agglomerate
      agglom_success <- FALSE
      tryCatch({
        ps_glom <- phyloseq::tax_glom(ps, taxrank = input$taxLevel)
        ps_safe <- ps_glom
        agglom_success <- TRUE
      }, error = function(e) {
        warning("Could not agglomerate taxa: ", conditionMessage(e))
      })
      
      # Transform to relative abundance
      ps_rel <- phyloseq::transform_sample_counts(ps_safe, function(x) x / sum(x))
      
      # Safely create data frame for plotting
      taxa_data <- tryCatch({
        phyloseq::psmelt(ps_rel)
      }, error = function(e) {
        warning("Error in psmelt: ", conditionMessage(e))
        # Create a basic data frame as fallback
        otu <- phyloseq::otu_table(ps_rel)
        tax <- phyloseq::tax_table(ps_rel)
        
        # Get dimensions
        if (phyloseq::taxa_are_rows(otu)) {
          n_taxa <- nrow(otu)
          n_samples <- ncol(otu)
          sample_names <- colnames(otu)
        } else {
          n_taxa <- ncol(otu)
          n_samples <- nrow(otu)
          sample_names <- rownames(otu)
        }
        
        # Build a simplified data frame
        data.frame(
          Sample = rep(sample_names, each = n_taxa),
          OTU = rep(phyloseq::taxa_names(ps_rel), times = n_samples),
          Abundance = as.vector(phyloseq::otu_table(ps_rel)),
          stringsAsFactors = FALSE
        )
      })
      
      # Make sure taxa_data has the taxonomy column we need
      if (!input$taxLevel %in% colnames(taxa_data)) {
        # Add the column from tax_table if possible
        tax <- phyloseq::tax_table(ps_rel)
        if (input$taxLevel %in% colnames(tax)) {
          taxa_data[[input$taxLevel]] <- tax[taxa_data$OTU, input$taxLevel]
        } else {
          # Create a dummy column as last resort
          taxa_data[[input$taxLevel]] <- "Unknown"
        }
      }
      
      # Replace NAs in taxonomy column
      taxa_data[[input$taxLevel]] <- as.character(taxa_data[[input$taxLevel]])
      taxa_data[[input$taxLevel]][is.na(taxa_data[[input$taxLevel]])] <- "Unknown"
      
      # Safely get top taxa
      n_top <- min(input$topTaxa, length(unique(taxa_data$OTU)))
      if (n_top < 1) n_top <- 1
      
      # Safe method to get top taxa
      top_taxa <- tryCatch({
        taxa_sums_val <- phyloseq::taxa_sums(ps_rel)
        names(sort(taxa_sums_val, decreasing = TRUE)[1:n_top])
      }, error = function(e) {
        # Fallback method to get top taxa
        abundance_by_taxa <- aggregate(Abundance ~ OTU, data = taxa_data, sum)
        abundance_by_taxa$OTU[order(abundance_by_taxa$Abundance, decreasing = TRUE)[1:n_top]]
      })
      
      # Group remaining taxa as "Other" if selected
      if (input$otherCat) {
        # Make sure the taxonomy column is character
        taxa_data[[input$taxLevel]] <- as.character(taxa_data[[input$taxLevel]])
        # Apply the grouping safely
        taxa_data[[input$taxLevel]] <- ifelse(taxa_data$OTU %in% top_taxa, 
                                         taxa_data[[input$taxLevel]], "Other")
      } else {
        # Filter to show only top taxa
        taxa_data <- taxa_data[taxa_data$OTU %in% top_taxa, ]
      }
      
      # If we have no data after filtering, show message
      if (nrow(taxa_data) == 0) {
        return(plot_ly() %>% 
               layout(title = "No data after filtering",
                     annotations = list(
                       x = 0.5, y = 0.5, 
                       text = "No taxa match the filter criteria", 
                       showarrow = FALSE
                     )))
      }
      
      # Create plot based on plot type
      if (input$plotType == "Bar Plot") {
        # Safe grouping without dplyr (which can fail with rlang errors)
        taxa_sum <- tryCatch({
          # Try with dplyr first
          taxa_data %>%
            dplyr::group_by(Sample, !!rlang::sym(input$taxLevel)) %>%
            dplyr::summarize(Abundance = sum(Abundance), .groups = "drop")
        }, error = function(e) {
          # Fallback to base R aggregation
          agg <- aggregate(Abundance ~ Sample + get(input$taxLevel), data = taxa_data, sum)
          colnames(agg)[2] <- input$taxLevel
          agg
        })
        
        # Create bar plot safely
        p <- tryCatch({
          # Create basic plot with taxonomy coloring
          ggplot(taxa_sum, aes(x = Sample, y = Abundance, fill = !!rlang::sym(input$taxLevel))) +
            geom_bar(stat = "identity") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = "Sample", y = "Relative Abundance", fill = input$taxLevel) +
            scale_fill_viridis_d()
        }, error = function(e) {
          # Fallback without rlang evaluation
          ggplot(taxa_sum, aes_string(x = "Sample", y = "Abundance", fill = input$taxLevel)) +
            geom_bar(stat = "identity") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = "Sample", y = "Relative Abundance", fill = input$taxLevel) +
            scale_fill_viridis_d()
        })
        
        return(ggplotly(p))
        
      } else {  # Heatmap
        # Safe creation of heatmap data
        taxa_mat <- tryCatch({
          # Try tidyr pivot_wider approach
          taxa_wide <- taxa_data %>%
            dplyr::group_by(Sample, !!rlang::sym(input$taxLevel)) %>%
            dplyr::summarize(Abundance = sum(Abundance), .groups = "drop") %>%
            tidyr::pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)
          
          # Convert to matrix
          mat <- as.matrix(taxa_wide[, -1])
          rownames(mat) <- taxa_wide[[input$taxLevel]]
          mat
        }, error = function(e) {
          # Fallback to using base R reshape
          taxa_agg <- aggregate(Abundance ~ Sample + get(input$taxLevel), data = taxa_data, sum)
          colnames(taxa_agg)[2] <- input$taxLevel
          
          # Reshape to wide format
          taxa_wide <- reshape(taxa_agg, 
                              idvar = input$taxLevel, 
                              timevar = "Sample", 
                              direction = "wide")
          
          # Fix column names
          colnames(taxa_wide) <- gsub("Abundance\\.", "", colnames(taxa_wide))
          
          # Convert to matrix
          mat <- as.matrix(taxa_wide[, -1])
          rownames(mat) <- taxa_wide[[input$taxLevel]]
          
          # Replace NAs with zeros
          mat[is.na(mat)] <- 0
          mat
        })
        
        # Check if we have matrix data
        if (!is.matrix(taxa_mat) || nrow(taxa_mat) == 0 || ncol(taxa_mat) == 0) {
          return(plot_ly() %>% 
                 layout(title = "Unable to create heatmap",
                       annotations = list(
                         x = 0.5, y = 0.5, 
                         text = "Could not convert data to matrix format", 
                         showarrow = FALSE
                       )))
        }
        
        # Create heatmap safely
        tryCatch({
          plot_ly(
            x = colnames(taxa_mat),
            y = rownames(taxa_mat),
            z = taxa_mat,
            type = "heatmap",
            colorscale = "Viridis"
          ) %>%
            layout(
              title = paste("Abundance of Top", n_top, input$taxLevel),
              xaxis = list(title = ""),
              yaxis = list(title = input$taxLevel)
            )
        }, error = function(e) {
          # Return error plot
          plot_ly() %>% 
            layout(title = "Error creating heatmap",
                  annotations = list(
                    x = 0.5, y = 0.5, 
                    text = paste("Error:", conditionMessage(e)), 
                    showarrow = FALSE
                  ))
        })
      }
    }, error = function(e) {
      # Global error handler
      plot_ly() %>% 
        layout(title = "Error in taxonomy plot",
              annotations = list(
                x = 0.5, y = 0.5, 
                text = paste("Error:", conditionMessage(e)), 
                showarrow = FALSE
              ))
    })
  })
  
  output$taxaTable <- renderDataTable({
    ps <- filtered_ps()
    if (is.null(ps) || is.null(input$taxLevel)) return(NULL)
    
    # Check if taxonomy table exists
    has_tax_table <- tryCatch({
      tt <- phyloseq::tax_table(ps)
      !is.null(tt) && ncol(tt) > 0 && nrow(tt) > 0
    }, error = function(e) FALSE)
    
    if (!has_tax_table) {
      return(NULL)
    }
    
    # Check if requested taxonomy level exists
    if (!input$taxLevel %in% colnames(phyloseq::tax_table(ps))) {
      return(NULL)
    }
    
    # Try to agglomerate
    tryCatch({
      # Agglomerate at selected taxonomic level
      ps_glom <- tax_glom(ps, taxrank = input$taxLevel)
      
      # Transform to relative abundance
      ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
      
      # Extract taxonomy table
      tax_table <- as.data.frame(tax_table(ps_rel))
      
      # Add abundance information
      tax_abund <- data.frame(
        Mean_Abundance = rowMeans(otu_table(ps_rel)),
        tax_table
      )
      
      # Sort by abundance
      tax_abund <- tax_abund[order(tax_abund$Mean_Abundance, decreasing = TRUE), ]
      
      # Format abundance column
      tax_abund$Mean_Abundance <- round(tax_abund$Mean_Abundance * 100, 2)
      
      # Add prevalence information
      prevalence <- apply(otu_table(ps_rel) > 0, 1, sum) / nsamples(ps_rel) * 100
      tax_abund$Prevalence <- round(prevalence[rownames(tax_abund)], 1)
      
      DT::datatable(tax_abund, 
                    options = list(pageLength = 10),
                    caption = paste("Taxonomy at", input$taxLevel, "level"))
    }, error = function(e) {
      # Return NULL if there was an error
      warning("Error creating taxonomy table: ", conditionMessage(e))
      return(NULL)
    })
  })
  
  # ASV Table tab outputs
  output$asvTable <- renderDataTable({
    ps <- filtered_ps()
    if (is.null(ps)) return(NULL)
    
    # Extract ASV table
    asv_table <- as.data.frame(t(otu_table(ps)))
    
    # Add taxonomy information safely
    tryCatch({
      tax_table <- as.data.frame(tax_table(ps))
      
      # Check if tax_table has rows for all ASVs
      missing_taxa <- setdiff(rownames(asv_table), rownames(tax_table))
      
      if (length(missing_taxa) > 0) {
        # Create empty taxonomy for missing ASVs
        empty_tax <- matrix(NA, 
                           nrow = length(missing_taxa), 
                           ncol = ncol(tax_table), 
                           dimnames = list(missing_taxa, colnames(tax_table)))
        
        # Bind with existing taxonomy
        tax_table <- rbind(tax_table, empty_tax)
      }
      
      # Add individual taxonomy columns to the ASV table
      for (level in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
        if (level %in% colnames(tax_table)) {
          asv_table[[level]] <- tax_table[rownames(asv_table), level]
        }
      }
      
      # Add combined taxonomy string
      tax_table$Taxonomy <- paste(
        tax_table$Kingdom, tax_table$Phylum, tax_table$Class,
        tax_table$Order, tax_table$Family, tax_table$Genus,
        tax_table$Species, sep = "; "
      )
      tax_table$Taxonomy <- gsub("NA", "", tax_table$Taxonomy)
      tax_table$Taxonomy <- gsub("; ; ", "; ", tax_table$Taxonomy)
      
      # Add taxonomy to ASV table (safely match rows)
      asv_table$Taxonomy <- tax_table[rownames(asv_table), "Taxonomy"]
    }, error = function(e) {
      # If taxonomy processing fails, add a placeholder
      warning("Error processing taxonomy: ", conditionMessage(e))
      asv_table$Taxonomy <- "Taxonomy unavailable"
    })
    
    # Add sequence length
    asv_table$Length <- nchar(rownames(asv_table))
    
    # Filter columns based on user-selected taxonomy level
    if (!is.null(input$taxonomyDisplayLevel) && input$taxonomyDisplayLevel != "All Levels") {
      # Show only the selected level
      taxonomy_cols <- c("Taxonomy", input$taxonomyDisplayLevel)
    } else {
      # Show all taxonomy levels
      taxonomy_cols <- c("Taxonomy", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    }
    
    # Order columns
    available_tax_cols <- intersect(taxonomy_cols, colnames(asv_table))
    col_order <- c(available_tax_cols, "Length", 
                   setdiff(colnames(asv_table), c(available_tax_cols, "Length")))
    asv_table <- asv_table[, col_order]
    
    # Create caption based on taxonomy level selection
    if (!is.null(input$taxonomyDisplayLevel) && input$taxonomyDisplayLevel != "All Levels") {
      caption_text <- paste("ASV Table with", input$taxonomyDisplayLevel, "Taxonomy")
    } else {
      caption_text <- "ASV Table with Complete Taxonomy"
    }
    
    DT::datatable(asv_table, 
                  options = list(pageLength = 10, scrollX = TRUE),
                  caption = caption_text)
  })
  
  output$seqInfo <- renderPrint({
    ps <- filtered_ps()
    if (is.null(ps)) return(NULL)
    
    # Get ASV sequences
    asv_seqs <- colnames(otu_table(ps))
    
    # Calculate sequence length distribution
    seq_lengths <- nchar(asv_seqs)
    
    cat("ASV Sequence Information:\n\n")
    cat("Total ASVs:", length(asv_seqs), "\n")
    cat("Min length:", min(seq_lengths), "\n")
    cat("Max length:", max(seq_lengths), "\n")
    cat("Mean length:", round(mean(seq_lengths), 1), "\n")
    cat("Median length:", median(seq_lengths), "\n\n")
    
    # Distribution of lengths
    cat("Sequence length distribution:\n")
    print(table(seq_lengths))
    
    # Show a few example sequences
    cat("\nExample sequences (first 5 ASVs):\n")
    for (i in 1:min(5, length(asv_seqs))) {
      cat("ASV", i, ":", asv_seqs[i], "\n")
    }
  })
}

# Create a simplified UI version to run first
ui_simple <- dashboardPage(
  dashboardHeader(title = "DADA2 Results (Simple)"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "overview",
              fluidRow(
                box(width = 12, title = "DADA2 Pipeline", status = "primary",
                    "This is a simplified version of the dashboard. Data is being loaded...")
              )
      )
    )
  )
)

# Run the application with simplified approach
app <- shinyApp(ui = ui, server = server)
runApp(app, port = 4321, launch.browser = TRUE)