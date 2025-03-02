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
library(d3Tree)          # For interactive taxonomic tree visualizations
library(iNEXT)           # For rarefaction curves
library(ggplotify)       # For converting base plots to ggplot objects
library(htmlwidgets)     # For saving interactive widgets
library(future)          # For parallel processing
library(promises)        # For asynchronous operations
library(data.table)      # For faster data manipulation
library(memoise)         # For caching function results
library(R.utils)         # Utility functions including checkpoints

# Set up parallel processing for performance
parallel_cores <- min(parallel::detectCores() - 1, 4)  # Keep at least 1 core free
if (parallel_cores > 1) {
  plan(multisession, workers = parallel_cores)
  options(future.globals.maxSize = 2 * 1024^3)  # 2GB limit for globals
  cat("Parallel processing enabled with", parallel_cores, "cores\n")
} else {
  plan(sequential)
  cat("Running in sequential mode (parallel processing not available)\n")
}

# Set up global caching system
cache_dir <- "dashboard_cache"
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

# Cache expensive operations
get_cached <- memoise::memoise

# Get multi-run options from environment variables
multi_run_options <- list(
  multi_run = as.logical(Sys.getenv("DADA2_MULTI_RUN", "FALSE")),
  run_pattern = Sys.getenv("DADA2_RUN_PATTERN", ""),
  run_column = Sys.getenv("DADA2_RUN_COLUMN", "")
)

# Only keep non-empty values
if (multi_run_options$run_pattern == "") multi_run_options$run_pattern <- NULL
if (multi_run_options$run_column == "") multi_run_options$run_column <- NULL

# Auto-detection of multi-run status if not explicitly specified
if (!multi_run_options$multi_run) {
  # Check if we can infer multi-run from sample names or metadata
  auto_detect_multi_run <- FALSE
  
  # Will be set to TRUE if we find evidence of multiple runs
  # This happens in the phyloseq_obj reactive function
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "DADA2 Results Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Sample Quality", tabName = "quality", icon = icon("chart-line")),
      menuItem("Alpha Diversity", tabName = "alpha", icon = icon("chart-bar")),
      menuItem("Beta Diversity", tabName = "beta", icon = icon("project-diagram")),
      menuItem("Multi-Run Analysis", tabName = "multirun", icon = icon("layer-group")),
      menuItem("Rarefaction", tabName = "rarefaction", icon = icon("area-chart")),
      menuItem("Taxonomy", tabName = "taxonomy", icon = icon("sitemap")),
      menuItem("Taxonomic Tree", tabName = "taxtree", icon = icon("tree")),
      menuItem("ASV Table", tabName = "asv", icon = icon("table")),
      menuItem("Export Data", tabName = "export", icon = icon("download"))
    ),
    
    hr(),
    
    # Filtering controls
    conditionalPanel(
      condition = "input.tab == 'alpha' || input.tab == 'beta' || input.tab == 'taxonomy' || input.tab == 'rarefaction' || input.tab == 'taxtree' || input.tab == 'asv' || input.tab == 'multirun'",
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
                box(width = 4, 
                    selectInput("alphaGroupColumn", "Group By:", choices = c("Auto-detect" = "auto"),
                               selected = "auto")),
                box(width = 4, 
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
                box(width = 3,
                    selectInput("betaMethod", "Distance Method:",
                                choices = c("Bray-Curtis", "Jaccard", "UniFrac", "Weighted UniFrac"),
                                selected = "Bray-Curtis")),
                box(width = 3,
                    selectInput("ordMethod", "Ordination Method:",
                                choices = c("PCoA", "NMDS", "t-SNE", "UMAP"),
                                selected = "PCoA")),
                box(width = 3,
                    selectInput("betaGroupColumn", "Group By:", choices = c("Auto-detect" = "auto"),
                               selected = "auto")),
                box(width = 3,
                    checkboxInput("ellipses", "Draw Confidence Ellipses", value = FALSE))
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
      ),
      
      # Rarefaction tab
      tabItem(tabName = "rarefaction",
              fluidRow(
                box(width = 12, title = "Rarefaction Curves", status = "info",
                    "Rarefaction curves show the number of observed ASVs as a function of sampling depth.")
              ),
              
              fluidRow(
                box(width = 3,
                    selectInput("rarefactionMethod", "Diversity metric:", 
                                choices = c("Observed (S)", "Shannon (H)", "Simpson (D)"),
                                selected = "Observed (S)")),
                box(width = 3,
                    sliderInput("rarefactionKnots", "Number of knots:", 
                                min = 5, max = 40, value = 10)),
                box(width = 3,
                    numericInput("rarefactionEndpoint", "Endpoint:", value = NULL, min = 1000)),
                box(width = 3,
                    actionButton("runRarefaction", "Generate Curves", icon = icon("refresh")))
              ),
              
              fluidRow(
                box(width = 12, title = "Rarefaction Plot", status = "info",
                    plotlyOutput("rarefactionPlot", height = "600px"))
              ),
              
              fluidRow(
                box(width = 12, title = "Download", status = "info",
                    downloadButton("downloadRarefaction", "Download Rarefaction Plot"))
              )
      ),
      
      # Taxonomic Tree tab
      tabItem(tabName = "taxtree",
              fluidRow(
                box(width = 12, title = "Interactive Taxonomic Tree", status = "info",
                    "Explore the taxonomic hierarchy in an interactive tree visualization.")
              ),
              
              fluidRow(
                box(width = 3,
                    numericInput("treePruneThreshold", "Abundance threshold (%):", 
                                value = 0.1, min = 0, max = 100, step = 0.1),
                    helpText("Only show taxa with relative abundance above threshold")),
                box(width = 3,
                    actionButton("generateTree", "Generate Tree", icon = icon("tree"))),
                box(width = 6,
                    downloadButton("downloadTree", "Download Tree HTML"),
                    br(), br(),
                    helpText("The downloaded HTML file can be opened in any web browser for interactive exploration."))
              ),
              
              fluidRow(
                box(width = 12, title = "Taxonomic Tree", status = "info", height = "800px",
                    d3TreeOutput("taxTree", height = "700px"))
              )
      ),
      
      # Multi-Run Analysis tab
      tabItem(tabName = "multirun",
              fluidRow(
                box(width = 12, title = "Multi-Run Analysis", status = "primary",
                    "This tab provides tools for analyzing batch effects and comparing results across multiple sequencing runs.")
              ),
              
              fluidRow(
                box(width = 6, title = "Run Detection", status = "info",
                    selectInput("runIdentifier", "Run identification method:",
                               choices = c("Auto-detect from sample names" = "auto", 
                                          "Use sample metadata column" = "metadata"),
                               selected = "auto"),
                    conditionalPanel(
                      condition = "input.runIdentifier == 'metadata'",
                      selectInput("runColumn", "Select run column from metadata:", choices = c("Auto-detect" = "auto"))
                    ),
                    conditionalPanel(
                      condition = "input.runIdentifier == 'auto'",
                      textInput("runPattern", "Run ID pattern in sample names (regex):", value = "run\\d+|R\\d+|Run\\d+|batch\\d+|Batch\\d+")
                    ),
                    actionButton("detectRuns", "Detect Runs", icon = icon("search")),
                    hr(),
                    h4("Detected Runs:"),
                    verbatimTextOutput("detectedRuns")
                ),
                
                box(width = 6, title = "Run Quality Metrics", status = "info",
                    plotlyOutput("runQualityMetrics")
                )
              ),
              
              fluidRow(
                tabBox(width = 12, 
                       tabPanel("Run Comparison",
                                fluidRow(
                                  column(3,
                                         selectInput("runCompareMethod", "Comparison Method:", 
                                                    choices = c("PCoA" = "pcoa", 
                                                               "NMDS" = "nmds", 
                                                               "t-SNE" = "tsne"),
                                                    selected = "pcoa")
                                  ),
                                  column(3,
                                         selectInput("runDistMethod", "Distance Metric:", 
                                                    choices = c("Bray-Curtis" = "bray", 
                                                               "Jaccard" = "jaccard", 
                                                               "UniFrac" = "unifrac", 
                                                               "Weighted UniFrac" = "wunifrac"),
                                                    selected = "bray")
                                  ),
                                  column(3,
                                         selectInput("runColorBy", "Color by:", 
                                                    choices = c("Run" = "run", 
                                                               "Sample Metadata" = "metadata"),
                                                    selected = "run")
                                  ),
                                  column(3,
                                         conditionalPanel(
                                           condition = "input.runColorBy == 'metadata'",
                                           selectInput("runColorMetadata", "Select metadata column:",
                                                      choices = c("Auto-detect" = "auto"))
                                         )
                                  )
                                ),
                                plotlyOutput("runBetaDiversity", height = "500px")
                       ),
                       
                       tabPanel("ASV Overlap",
                                fluidRow(
                                  column(4,
                                         selectInput("asvOverlapMethod", "Visualization Method:", 
                                                    choices = c("UpSet Plot" = "upset",
                                                               "Venn Diagram (up to 5 runs)" = "venn",
                                                               "Presence/Absence Heatmap" = "heatmap"),
                                                    selected = "upset")
                                  ),
                                  column(4,
                                         checkboxInput("showPrevalentASVsOnly", "Show only prevalent ASVs", value = TRUE),
                                         conditionalPanel(
                                           condition = "input.showPrevalentASVsOnly == true",
                                           sliderInput("asvPrevalenceThreshold", "Minimum prevalence (% of samples):",
                                                      min = 0, max = 100, value = 10)
                                         )
                                  ),
                                  column(4,
                                         checkboxInput("filterRunsForOverlap", "Filter specific runs", value = FALSE),
                                         conditionalPanel(
                                           condition = "input.filterRunsForOverlap == true",
                                           uiOutput("runSelectionCheckboxes")
                                         )
                                  )
                                ),
                                plotOutput("asvOverlapPlot", height = "500px")
                       ),
                       
                       tabPanel("Taxonomy by Run",
                                fluidRow(
                                  column(3,
                                         selectInput("taxRunLevel", "Taxonomic Level:", 
                                                    choices = c("Phylum", "Class", "Order", "Family", "Genus"),
                                                    selected = "Phylum")
                                  ),
                                  column(3,
                                         selectInput("taxRunVizType", "Visualization Type:", 
                                                    choices = c("Stacked Bar" = "bar", 
                                                               "Heatmap" = "heatmap", 
                                                               "PCoA" = "pcoa"),
                                                    selected = "bar")
                                  ),
                                  column(3,
                                         sliderInput("taxRunTopN", "Top N Taxa:", 
                                                    min = 5, max = 30, value = 10)
                                  ),
                                  column(3,
                                         conditionalPanel(
                                           condition = "input.taxRunVizType == 'heatmap' || input.taxRunVizType == 'pcoa'",
                                           selectInput("taxRunClusterMethod", "Clustering Method:",
                                                      choices = c("None" = "none", 
                                                                 "Ward.D2" = "ward.D2", 
                                                                 "Complete" = "complete", 
                                                                 "Average" = "average"),
                                                      selected = "ward.D2")
                                         )
                                  )
                                ),
                                plotlyOutput("taxonomyRunPlot", height = "500px")
                       ),
                       
                       tabPanel("Run Batch Effects",
                                fluidRow(
                                  column(3,
                                         selectInput("batchEffectMethod", "Analysis Method:", 
                                                    choices = c("Beta Dispersion" = "betadisper",
                                                               "PERMANOVA" = "adonis",
                                                               "Run vs. Sample Type" = "interaction",
                                                               "RDA" = "rda",
                                                               "Mantel Test" = "mantel",
                                                               "ANOSIM" = "anosim"),
                                                    selected = "betadisper")
                                  ),
                                  column(3,
                                         conditionalPanel(
                                           condition = "input.batchEffectMethod == 'interaction'",
                                           selectInput("interactionFactor", "Sample Type Factor:",
                                                      choices = c("Auto-detect" = "auto"))
                                         )
                                  ),
                                  column(3,
                                         selectInput("batchDistMethod", "Distance Metric:", 
                                                    choices = c("Bray-Curtis" = "bray", 
                                                               "Jaccard" = "jaccard"),
                                                    selected = "bray")
                                  ),
                                  column(3,
                                         actionButton("runBatchTest", "Run Analysis", icon = icon("calculator"))
                                  )
                                ),
                                hr(),
                                h4("Batch Effect Test Results:"),
                                verbatimTextOutput("batchTestResults"),
                                plotlyOutput("batchEffectPlot", height = "400px")
                       ),
                       
                       tabPanel("Run Normalization Options",
                                fluidRow(
                                  column(4,
                                         selectInput("runNormMethod", "Normalization Method:", 
                                                    choices = c("None" = "none",
                                                               "Relative Abundance" = "relative",
                                                               "CSS (cumulative sum scaling)" = "css",
                                                               "RLE (DESeq2)" = "rle",
                                                               "TMM (edgeR)" = "tmm",
                                                               "CLR (centered log-ratio)" = "clr",
                                                               "ALR (additive log-ratio)" = "alr",
                                                               "VST (variance stabilizing)" = "vst"),
                                                    selected = "relative")
                                  ),
                                  column(4,
                                         checkboxInput("applyRunNormalization", "Apply normalization to all analyses", value = FALSE)
                                  ),
                                  column(4,
                                         conditionalPanel(
                                           condition = "input.runNormMethod != 'none' && input.runNormMethod != 'relative'",
                                           actionButton("computeNormalization", "Compute Normalization", icon = icon("sync"))
                                         )
                                  )
                                ),
                                hr(),
                                plotlyOutput("normalizationComparisonPlot", height = "400px"),
                                verbatimTextOutput("normalizationSummary")
                       )
                )
              )
      ),
      
      # Export Data tab
      tabItem(tabName = "export",
              fluidRow(
                box(width = 12, title = "Export Data and Visualizations", status = "info",
                    "Download tables, plots and reports from your DADA2 analysis.")
              ),
              
              fluidRow(
                tabBox(width = 12, 
                       tabPanel("Tables",
                                fluidRow(
                                  box(width = 6, 
                                      h4("ASV Table"),
                                      selectInput("exportAsvFormat", "Format:",
                                                 choices = c("CSV", "TSV", "Excel"),
                                                 selected = "CSV"),
                                      downloadButton("downloadAsvTable", "Download ASV Table")
                                  ),
                                  box(width = 6,
                                      h4("Taxonomy Table"),
                                      selectInput("exportTaxFormat", "Format:",
                                                 choices = c("CSV", "TSV", "Excel"),
                                                 selected = "CSV"),
                                      downloadButton("downloadTaxTable", "Download Taxonomy Table")
                                  )
                                ),
                                fluidRow(
                                  box(width = 6,
                                      h4("Sample Metadata"),
                                      downloadButton("downloadMetadata", "Download Metadata")
                                  ),
                                  box(width = 6,
                                      h4("Read Tracking"),
                                      downloadButton("downloadReadTracking", "Download Read Tracking Data")
                                  )
                                )
                       ),
                       tabPanel("Plots",
                                fluidRow(
                                  box(width = 4,
                                      h4("Alpha Diversity Plots"),
                                      selectInput("exportPlotFormat", "Format:",
                                                 choices = c("PNG", "PDF", "SVG", "HTML"),
                                                 selected = "PNG"),
                                      downloadButton("downloadAlphaPlot", "Download Alpha Diversity Plot")
                                  ),
                                  box(width = 4,
                                      h4("Beta Diversity Plots"),
                                      downloadButton("downloadBetaPlot", "Download Beta Diversity Plot")
                                  ),
                                  box(width = 4,
                                      h4("Taxonomy Barplot"),
                                      downloadButton("downloadTaxPlot", "Download Taxonomy Plot")
                                  )
                                )
                       ),
                       tabPanel("R Objects",
                                fluidRow(
                                  box(width = 6,
                                      h4("Phyloseq Object"),
                                      downloadButton("downloadPhyloseq", "Download Phyloseq (.RDS)"),
                                      helpText("Save the entire phyloseq object to use in R")
                                  ),
                                  box(width = 6,
                                      h4("BIOM Format"),
                                      downloadButton("downloadBiom", "Download BIOM file"),
                                      helpText("Export data in BIOM format for use with other tools")
                                  )
                                )
                       ),
                       tabPanel("Report",
                                fluidRow(
                                  box(width = 12,
                                      h4("Generate Analysis Report"),
                                      p("Create a comprehensive HTML report with plots and tables from your analysis."),
                                      checkboxGroupInput("reportSections", "Include sections:",
                                                        choices = c("Overview", "Quality Control", "Alpha Diversity", 
                                                                   "Beta Diversity", "Taxonomy", "ASV Table"),
                                                        selected = c("Overview", "Alpha Diversity", "Beta Diversity", "Taxonomy")),
                                      downloadButton("downloadReport", "Generate and Download Report")
                                  )
                                )
                       )
                )
              )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Optimized phyloseq object loading with caching
  phyloseq_obj <- reactive({
    # Path to saved phyloseq object
    ps_path <- "results/phyloseq_object.rds"
    ps_cache_path <- file.path(cache_dir, "phyloseq_processed.rds")
    
    # Check if file exists
    if (!file.exists(ps_path)) {
      return(NULL)
    }
    
    # Check if we have a cached processed version
    if (file.exists(ps_cache_path)) {
      # Check if cache is newer than source file
      if (file.info(ps_cache_path)$mtime > file.info(ps_path)$mtime) {
        # Cache is valid - load it
        withProgress(message = 'Loading cached phyloseq object...', value = 0.5, {
          ps_cached <- readRDS(ps_cache_path)
          
          # Auto-detect multi-run if not specified
          if (!multi_run_options$multi_run && exists("auto_detect_multi_run")) {
            # Check sample names for run patterns
            sample_names <- sample_names(ps_cached)
            
            # Look for common run patterns in sample names
            run_patterns <- c("run\\d+", "batch\\d+", "seq\\d+", "lane\\d+", "flow\\d+")
            
            for (pattern in run_patterns) {
              matches <- grep(pattern, sample_names, ignore.case = TRUE)
              if (length(matches) > 0) {
                # Found potential run indicators
                run_ids <- regmatches(sample_names[matches], 
                                     regexpr(pattern, sample_names[matches], ignore.case = TRUE))
                
                # If we have multiple distinct run IDs, enable multi-run mode
                if (length(unique(run_ids)) > 1) {
                  auto_detect_multi_run <<- TRUE
                  cat("Auto-detected multiple runs in sample names!\n")
                  cat("Enabling multi-run analysis mode.\n")
                  break
                }
              }
            }
            
            # Also check metadata if available
            if (!auto_detect_multi_run) {
              has_sample_data <- tryCatch({
                !is.null(sample_data(ps_cached))
              }, error = function(e) FALSE)
              
              if (has_sample_data) {
                sample_df <- as.data.frame(sample_data(ps_cached))
                
                # Look for columns that might indicate runs
                run_cols <- grep("run|batch|seq|lane|flow", 
                                colnames(sample_df), ignore.case = TRUE, value = TRUE)
                
                for (col in run_cols) {
                  if (length(unique(sample_df[[col]])) > 1) {
                    auto_detect_multi_run <<- TRUE
                    cat("Auto-detected multiple runs in sample metadata column:", col, "\n")
                    cat("Enabling multi-run analysis mode.\n")
                    break
                  }
                }
              }
            }
          }
          
          return(ps_cached)
        })
      }
    }
    
    # No valid cache - process the phyloseq object and cache it
    withProgress(message = 'Processing phyloseq object...', value = 0.2, {
      # Try to safely load and process the phyloseq object
      ps <- tryCatch({
        # Load the phyloseq object
        incProgress(0.2, detail = "Loading data")
        ps <- readRDS(ps_path)
        
        # Safe check if taxa table exists
        incProgress(0.2, detail = "Checking taxonomy")
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
        
        incProgress(0.2, detail = "Finalizing")
        ps
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
      
      # Cache the processed object for future use
      if (!is.null(ps)) {
        incProgress(0.2, detail = "Caching for future use")
        saveRDS(ps, ps_cache_path)
        
        # Also do multi-run auto-detection here
        if (!multi_run_options$multi_run && exists("auto_detect_multi_run")) {
          # Similar code as above for the uncached path
          # Check sample names for run patterns
          sample_names <- sample_names(ps)
          
          # Look for common run patterns in sample names
          run_patterns <- c("run\\d+", "batch\\d+", "seq\\d+", "lane\\d+", "flow\\d+")
          
          for (pattern in run_patterns) {
            matches <- grep(pattern, sample_names, ignore.case = TRUE)
            if (length(matches) > 0) {
              # Found potential run indicators
              run_ids <- regmatches(sample_names[matches], 
                                   regexpr(pattern, sample_names[matches], ignore.case = TRUE))
              
              # If we have multiple distinct run IDs, enable multi-run mode
              if (length(unique(run_ids)) > 1) {
                auto_detect_multi_run <<- TRUE
                cat("Auto-detected multiple runs in sample names!\n")
                cat("Enabling multi-run analysis mode.\n")
                break
              }
            }
          }
          
          # Also check metadata if available
          if (!auto_detect_multi_run) {
            has_sample_data <- tryCatch({
              !is.null(sample_data(ps))
            }, error = function(e) FALSE)
            
            if (has_sample_data) {
              sample_df <- as.data.frame(sample_data(ps))
              
              # Look for columns that might indicate runs
              run_cols <- grep("run|batch|seq|lane|flow", 
                              colnames(sample_df), ignore.case = TRUE, value = TRUE)
              
              for (col in run_cols) {
                if (length(unique(sample_df[[col]])) > 1) {
                  auto_detect_multi_run <<- TRUE
                  cat("Auto-detected multiple runs in sample metadata column:", col, "\n")
                  cat("Enabling multi-run analysis mode.\n")
                  break
                }
              }
            }
          }
        }
      }
      
      return(ps)
    })
  })
  
  # Function to check if a phyloseq object is large
  is_large_phyloseq <- reactive({
    ps <- phyloseq_obj()
    if (is.null(ps)) return(FALSE)
    
    # Consider large if more than these thresholds
    large_samples_threshold <- 100
    large_taxa_threshold <- 5000
    
    n_samples <- tryCatch({
      nsamples(ps)
    }, error = function(e) {
      0
    })
    
    n_taxa <- tryCatch({
      ntaxa(ps)
    }, error = function(e) {
      0
    })
    
    return(n_samples > large_samples_threshold || n_taxa > large_taxa_threshold)
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
  
  # Update UI elements based on sample data when phyloseq object is loaded
  observe({
    ps <- phyloseq_obj()
    if (is.null(ps)) return()
    
    # Check if sample data exists
    has_sample_data <- tryCatch({
      !is.null(sample_data(ps))
    }, error = function(e) FALSE)
    
    if (has_sample_data) {
      # Extract sample data
      sample_df <- as.data.frame(sample_data(ps))
      
      # Get column names
      meta_cols <- colnames(sample_df)
      
      # Identify potential grouping columns (categorical or limited numeric values)
      grouping_cols <- c()
      
      for (col in meta_cols) {
        # Skip columns that are likely not categorical (too many unique values)
        n_unique <- length(unique(sample_df[[col]]))
        if (n_unique > 1 && n_unique < min(20, nrow(sample_df))) {
          grouping_cols <- c(grouping_cols, col)
        }
      }
      
      # If we found potential grouping columns, update the selection inputs
      if (length(grouping_cols) > 0) {
        # Create named list with Auto-detect option first
        group_choices <- c("Auto-detect" = "auto", setNames(grouping_cols, grouping_cols))
        
        # Update alpha diversity group selector
        updateSelectInput(session, "alphaGroupColumn", 
                         choices = group_choices,
                         selected = "auto")
        
        # Update beta diversity group selector
        updateSelectInput(session, "betaGroupColumn", 
                         choices = group_choices,
                         selected = "auto")
      }
    }
  })
  
  # Optimized filtering with progress tracking and caching
  filtered_ps <- reactive({
    # Parameters that affect filtering
    min_reads <- input$minReads
    prevalence <- input$prevalence
    
    # Create a cache key based on parameters
    cache_key <- paste0("filtered_ps_", min_reads, "_", prevalence)
    cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
    
    # Get the original phyloseq object
    ps <- phyloseq_obj()
    if (is.null(ps)) return(NULL)
    
    # Check if we have a cache hit
    if (file.exists(cache_file)) {
      # Check if cache is newer than source phyloseq
      ps_path <- "results/phyloseq_object.rds"
      if (file.exists(ps_path) && file.info(cache_file)$mtime > file.info(ps_path)$mtime) {
        # Use cached result
        return(readRDS(cache_file))
      }
    }
    
    # Extra validation to ensure phyloseq object is valid
    valid_phyloseq <- tryCatch({
      # Basic validation: check if it has an OTU table
      has_otu <- !is.null(phyloseq::otu_table(ps))
      has_otu
    }, error = function(e) FALSE)
    
    if (!valid_phyloseq) {
      warning("Invalid phyloseq object detected. Creating a basic one.")
      return(NULL)
    }
    
    # Check if this is a large dataset
    large_dataset <- is_large_phyloseq()
    
    # Create a progress indicator for large datasets
    if (large_dataset) {
      withProgress(message = 'Filtering large dataset...', value = 0, {
        filtered <- perform_filtering(ps, min_reads, prevalence)
        
        # Cache the result for future use
        saveRDS(filtered, cache_file)
        return(filtered)
      })
    } else {
      # For smaller datasets, just filter normally
      filtered <- perform_filtering(ps, min_reads, prevalence)
      
      # Cache the result
      saveRDS(filtered, cache_file)
      return(filtered)
    }
  })
  
  # Helper function to perform the actual filtering
  perform_filtering <- function(ps, min_reads, prevalence_threshold, progress_callback = NULL) {
    # Try all filtering operations with error handling
    tryCatch({
      # Initialize progress if we have a callback
      if (!is.null(progress_callback)) {
        progress_callback(0.1, "Starting filtering")
      }
      
      # Filter by read depth
      if (!is.null(min_reads)) {
        tryCatch({
          if (!is.null(progress_callback)) {
            progress_callback(0.3, "Filtering by read depth")
          }
          
          # Use data.table for faster operations with large datasets
          sample_sums_val <- phyloseq::sample_sums(ps)
          ps <- phyloseq::prune_samples(sample_sums_val >= min_reads, ps)
          
          # Force garbage collection to free memory
          gc()
        }, error = function(e) {
          warning("Error in prune_samples: ", conditionMessage(e))
        })
      }
      
      # Filter by prevalence with careful validation
      if (!is.null(prevalence_threshold)) {
        tryCatch({
          if (!is.null(progress_callback)) {
            progress_callback(0.6, "Filtering by prevalence")
          }
          
          # Get OTU table with correct orientation
          otu <- phyloseq::otu_table(ps)
          
          # Optimize prevalence calculation for large datasets
          if (is_large_phyloseq()) {
            # Use data.table for faster operations
            if (phyloseq::taxa_are_rows(otu)) {
              otu_dt <- as.data.table(as.matrix(otu))
              prevalence <- sapply(otu_dt, function(x) sum(x > 0) / length(x) * 100)
            } else {
              otu_dt <- as.data.table(t(as.matrix(otu)))
              prevalence <- sapply(otu_dt, function(x) sum(x > 0) / length(x) * 100)
            }
          } else {
            # Standard calculation for smaller datasets
            if (phyloseq::taxa_are_rows(otu)) {
              prevalence <- apply(otu, 1, function(x) sum(x > 0) / length(x) * 100)
            } else {
              prevalence <- apply(otu, 2, function(x) sum(x > 0) / length(x) * 100)
            }
          }
          
          # Make sure prevalence vector has names matching taxa
          taxa_names <- phyloseq::taxa_names(ps)
          if (length(prevalence) == length(taxa_names)) {
            names(prevalence) <- taxa_names
            # Use prune_taxa with pre-calculated logical vector
            to_keep <- prevalence >= prevalence_threshold
            ps <- phyloseq::prune_taxa(to_keep, ps)
          }
          
          # Force garbage collection to free memory
          gc()
        }, error = function(e) {
          warning("Error in prevalence filtering: ", conditionMessage(e))
        })
      }
      
      if (!is.null(progress_callback)) {
        progress_callback(1.0, "Filtering complete")
      }
      
      return(ps)
    }, error = function(e) {
      warning("Error in filtered_ps: ", conditionMessage(e))
      return(ps)  # Return original if filtering fails
    })
  }
  
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
    
    # Check for sample data
    has_sample_data <- tryCatch({
      !is.null(sample_data(ps))
    }, error = function(e) FALSE)
    
    has_grouping <- FALSE
    group_column <- NULL
    
    if (has_sample_data) {
      # Extract sample data
      sample_df <- as.data.frame(sample_data(ps))
      
      # Check if user selected a specific grouping column
      if (!is.null(input$alphaGroupColumn) && input$alphaGroupColumn != "auto") {
        # User selected a specific column
        if (input$alphaGroupColumn %in% colnames(sample_df)) {
          group_column <- input$alphaGroupColumn
          has_grouping <- TRUE
        }
      } else {
        # Auto-detect grouping column
        potential_group_cols <- c("Group", "TreatmentGroup", "Treatment", "Condition", 
                                 "SampleType", "SampleGroup", "ExperimentalGroup", "Site", 
                                 "Location", "Patient", "Subject", "TimePoint", "Time")
        
        # Find the first available column that could be used for grouping
        for (col in potential_group_cols) {
          if (col %in% colnames(sample_df)) {
            # Check if it has multiple values (not useful if all samples have same value)
            if (length(unique(sample_df[[col]])) > 1 && length(unique(sample_df[[col]])) < nrow(sample_df)) {
              group_column <- col
              has_grouping <- TRUE
              break
            }
          }
        }
      }
      
      # If we found a grouping column, add it to the alpha diversity data
      if (has_grouping) {
        alpha_df$Group <- sample_df[alpha_df$Sample, group_column]
      }
    }
    
    # Create plot - either grouped or simple bar plot
    if (has_grouping) {
      # Create grouped plots
      
      # First, create a boxplot by group
      p_box <- ggplot(alpha_df, aes(x = Group, y = Diversity, fill = Group)) +
        geom_boxplot(outlier.shape = NA) +  # Hide outliers as they'll be shown in the jitter
        geom_jitter(width = 0.2, height = 0, alpha = 0.7, aes(text = Sample)) +
        theme_minimal() +
        labs(x = group_column, 
             y = paste(input$alphaMethod, "Diversity"), 
             title = paste(input$alphaMethod, "Diversity by", group_column)) +
        scale_fill_viridis_d() +
        theme(legend.position = "none")  # Hide redundant legend
      
      # Second, create a bar plot colored by group
      p_bar <- ggplot(alpha_df, aes(x = reorder(Sample, -Diversity), y = Diversity, fill = Group)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Sample", 
             y = paste(input$alphaMethod, "Diversity"), 
             fill = group_column,
             title = paste(input$alphaMethod, "Diversity by Sample")) +
        scale_fill_viridis_d()
      
      # Apply log scale if selected
      if (input$logScale) {
        p_box <- p_box + scale_y_log10()
        p_bar <- p_bar + scale_y_log10()
      }
      
      # Arrange the two plots
      p <- plotly::subplot(
        ggplotly(p_box, tooltip = c("text", "y")),
        ggplotly(p_bar),
        nrows = 1,
        shareY = TRUE,
        widths = c(0.4, 0.6)
      ) %>% 
        plotly::layout(
          showlegend = TRUE, 
          legend = list(orientation = "h", y = -0.1)
        )
      
      return(p)
    } else {
      # Simple bar plot without grouping
      p <- ggplot(alpha_df, aes(x = reorder(Sample, -Diversity), y = Diversity)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Sample", y = paste(input$alphaMethod, "Diversity"))
      
      # Apply log scale if selected
      if (input$logScale) {
        p <- p + scale_y_log10()
      }
      
      return(ggplotly(p))
    }
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
    
    # Add sample names
    ord_data$Sample <- rownames(ord_data)
    
    # Variables for grouping
    has_grouping <- FALSE
    group_column <- NULL
    
    # Check for sample data
    has_sample_data <- tryCatch({
      !is.null(sample_data(ps))
    }, error = function(e) FALSE)
    
    if (has_sample_data) {
      # Extract sample data
      sample_df <- as.data.frame(sample_data(ps))
      
      # Check if user selected a specific grouping column
      if (!is.null(input$betaGroupColumn) && input$betaGroupColumn != "auto") {
        # User selected a specific column
        if (input$betaGroupColumn %in% colnames(sample_df)) {
          group_column <- input$betaGroupColumn
          has_grouping <- TRUE
        }
      } else {
        # Auto-detect grouping column
        potential_group_cols <- c("Group", "TreatmentGroup", "Treatment", "Condition", 
                                 "SampleType", "SampleGroup", "ExperimentalGroup", "Site", 
                                 "Location", "Patient", "Subject", "TimePoint", "Time")
        
        # Find the first available column that could be used for grouping
        for (col in potential_group_cols) {
          if (col %in% colnames(sample_df)) {
            # Check if it has multiple values (not useful if all samples have same value)
            if (length(unique(sample_df[[col]])) > 1 && length(unique(sample_df[[col]])) < nrow(sample_df)) {
              group_column <- col
              has_grouping <- TRUE
              break
            }
          }
        }
      }
      
      # If we found a grouping column, add it to the ordination data
      if (has_grouping) {
        ord_data$Group <- sample_df[ord_data$Sample, group_column]
        # Update hover text to include group
        ord_data$hover_text <- paste("Sample:", ord_data$Sample, 
                                    "<br>", group_column, ":", ord_data$Group)
      }
    }
    
    # If no hover text created yet, use just the sample name
    if (!"hover_text" %in% colnames(ord_data)) {
      ord_data$hover_text <- ord_data$Sample
    }
    
    # Create plot with or without grouping
    if (has_grouping) {
      p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, text = hover_text, color = Group)) +
        geom_point(size = 3, alpha = 0.7) +
        theme_minimal() +
        labs(x = axis_labels[1], y = axis_labels[2], color = group_column) +
        scale_color_viridis_d()
      
      # Add confidence ellipses if requested and if there are enough points per group
      if (input$ellipses) {
        # Check each group has at least 3 samples (required for ellipse calculation)
        group_counts <- table(ord_data$Group)
        enough_samples <- all(group_counts >= 3)
        
        if (enough_samples) {
          tryCatch({
            # Try to add ellipses
            p <- p + stat_ellipse(aes(fill = Group), geom = "polygon", 
                                 alpha = 0.2, level = 0.95, type = "t")
          }, error = function(e) {
            # If ellipse calculation fails, just continue without ellipses
            warning("Could not calculate confidence ellipses: ", conditionMessage(e))
          })
        }
      }
    } else {
      p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, text = hover_text)) +
        geom_point(size = 3, alpha = 0.7, color = "steelblue") +
        theme_minimal() +
        labs(x = axis_labels[1], y = axis_labels[2])
    }
    
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
  
  # Optimized rarefaction tab functionality with parallel processing
  # Reactive to store rarefaction data
  rarefaction_data <- eventReactive(input$runRarefaction, {
    ps <- filtered_ps()
    if (is.null(ps)) return(NULL)
    
    # Create a cache key based on parameters
    method <- input$rarefactionMethod
    knots <- input$rarefactionKnots
    endpoint <- input$rarefactionEndpoint
    min_reads <- input$minReads
    prevalence <- input$prevalence
    
    cache_key <- paste0("rarefaction_", method, "_", knots, "_", endpoint, "_", min_reads, "_", prevalence)
    cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
    
    # Check if we have a cache hit
    if (file.exists(cache_file)) {
      # Load from cache
      withProgress(message = 'Loading cached rarefaction curves...', value = 0.5, {
        return(readRDS(cache_file))
      })
    }
    
    # Calculate with progress indicator
    withProgress(message = 'Calculating rarefaction curves...', value = 0, {
      # Get OTU table
      setProgress(value = 0.1, detail = "Preparing OTU table")
      otu <- phyloseq::otu_table(ps)
      
      # Make sure OTU table is in the right format (samples as rows)
      if (phyloseq::taxa_are_rows(otu)) {
        otu <- t(otu)
      }
      
      # Convert to regular matrix and make sure it's integers
      otu_mat <- as.matrix(otu)
      mode(otu_mat) <- "integer"
      
      # Set endpoint if not provided
      if (is.null(endpoint) || is.na(endpoint) || endpoint <= 0) {
        # Use max read count as endpoint
        endpoint <- max(rowSums(otu_mat))
      }
      
      # Determine diversity type
      q_value <- switch(method,
                       "Observed (S)" = 0,
                       "Shannon (H)" = 1,
                       "Simpson (D)" = 2)
      
      # Check if dataset is large
      is_large <- nrow(otu_mat) > 50 || ncol(otu_mat) > 1000
      
      if (is_large) {
        setProgress(value = 0.2, detail = "Using optimized processing for large dataset")
        
        # For large datasets, process in smaller chunks to avoid memory issues
        # First downsample to a manageable number of samples if necessary
        max_samples <- 50
        if (nrow(otu_mat) > max_samples) {
          # Sample indices randomly
          set.seed(123) # For reproducibility
          sample_indices <- sample(1:nrow(otu_mat), max_samples)
          otu_mat_downsampled <- otu_mat[sample_indices, ]
          setProgress(value = 0.3, detail = paste("Downsampled to", max_samples, "samples"))
        } else {
          otu_mat_downsampled <- otu_mat
        }
        
        # Process with smaller knots for large datasets
        knots_to_use <- min(knots, 20)
        setProgress(value = 0.4, detail = paste("Using", knots_to_use, "knots for efficiency"))
        
        # Try with parallel processing if available
        if (parallel_cores > 1 && requireNamespace("future", quietly = TRUE) && 
            requireNamespace("future.apply", quietly = TRUE)) {
          
          setProgress(value = 0.5, detail = "Using parallel processing")
          
          # Using future for parallel processing
          future::plan(future::multisession, workers = parallel_cores)
          
          # Process different samples in parallel
          sample_groups <- split(1:nrow(otu_mat_downsampled), 
                                ceiling(seq_along(1:nrow(otu_mat_downsampled)) / 
                                       ceiling(nrow(otu_mat_downsampled) / parallel_cores)))
          
          # Process each group in parallel
          setProgress(value = 0.6, detail = "Running parallel iNEXT calculations")
          
          tryCatch({
            out <- iNEXT::iNEXT(otu_mat_downsampled, q = q_value, 
                              endpoint = endpoint, 
                              knots = knots_to_use)
            
            setProgress(value = 0.9, detail = "Finalizing results")
          }, error = function(e) {
            warning("Error in rarefaction: ", conditionMessage(e))
            return(NULL)
          })
        } else {
          # Standard processing for systems without parallel support
          setProgress(value = 0.5, detail = "Calculating rarefaction curves")
          
          tryCatch({
            out <- iNEXT::iNEXT(otu_mat_downsampled, q = q_value, 
                              endpoint = endpoint, 
                              knots = knots_to_use)
            
            setProgress(value = 0.9, detail = "Finalizing results")
          }, error = function(e) {
            warning("Error in rarefaction: ", conditionMessage(e))
            return(NULL)
          })
        }
      } else {
        # Standard processing for smaller datasets
        setProgress(value = 0.5, detail = "Calculating rarefaction curves")
        
        tryCatch({
          out <- iNEXT::iNEXT(otu_mat, q = q_value, 
                            endpoint = endpoint, 
                            knots = knots)
          
          setProgress(value = 0.9, detail = "Finalizing results")
        }, error = function(e) {
          warning("Error in rarefaction: ", conditionMessage(e))
          return(NULL)
        })
      }
      
      # Cache the result if successful
      if (!is.null(out)) {
        saveRDS(out, cache_file)
      }
      
      # Clean up and return
      gc() # Force garbage collection
      return(out)
    })
  })
  
  # Render rarefaction plot
  output$rarefactionPlot <- renderPlotly({
    # Get rarefaction data
    out <- rarefaction_data()
    if (is.null(out)) {
      return(plot_ly() %>% 
             layout(title = "Rarefaction analysis failed",
                   annotations = list(
                     x = 0.5, y = 0.5, 
                     text = "Could not calculate rarefaction curves. Try different parameters.", 
                     showarrow = FALSE
                   )))
    }
    
    # Convert iNEXT plot to ggplot and then to plotly
    tryCatch({
      # Get the basic plot from iNEXT
      inext_plot <- iNEXT::ggiNEXT(out)
      
      # Convert to ggplot and improve aesthetics
      p <- inext_plot +
        theme_minimal() +
        labs(title = paste(input$rarefactionMethod, "Rarefaction Curves"),
             x = "Sample Size (Number of Reads)",
             y = "Diversity") +
        theme(legend.position = "right")
      
      # Convert to plotly
      ggplotly(p)
    }, error = function(e) {
      # If plotting fails, return basic error plot
      warning("Error creating rarefaction plot: ", conditionMessage(e))
      plot_ly() %>% 
        layout(title = "Error creating rarefaction plot",
              annotations = list(
                x = 0.5, y = 0.5, 
                text = paste("Error:", conditionMessage(e)), 
                showarrow = FALSE
              ))
    })
  })
  
  # Download handler for rarefaction plot
  output$downloadRarefaction <- downloadHandler(
    filename = function() {
      paste0("rarefaction_", input$rarefactionMethod, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      # Get rarefaction data and create plot
      out <- rarefaction_data()
      if (is.null(out)) {
        # Create a message plot if data is NULL
        p <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "No rarefaction data available") +
          theme_void()
      } else {
        # Create the plot
        p <- iNEXT::ggiNEXT(out) +
          theme_minimal() +
          labs(title = paste(input$rarefactionMethod, "Rarefaction Curves"),
               x = "Sample Size (Number of Reads)",
               y = "Diversity") +
          theme(legend.position = "right")
      }
      
      # Save plot
      ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
    }
  )
  
  # Taxonomic Tree tab functionality optimized for large datasets
  # Function to prepare taxonomy data in hierarchical format for d3Tree
  prepare_tree_data <- reactive({
    # Create a cache key based on parameters
    threshold <- input$treePruneThreshold
    min_reads <- input$minReads
    prevalence <- input$prevalence
    
    cache_key <- paste0("tax_tree_", threshold, "_", min_reads, "_", prevalence)
    cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
    
    # Check if we have a cache hit
    if (file.exists(cache_file)) {
      # Check if cache is newer than source phyloseq
      ps_path <- "results/phyloseq_object.rds"
      if (file.exists(ps_path) && file.info(cache_file)$mtime > file.info(ps_path)$mtime) {
        # Use cached result
        return(readRDS(cache_file))
      }
    }
    
    # Not in cache, need to compute
    ps <- filtered_ps()
    if (is.null(ps)) return(NULL)
    
    # Check if taxonomy table exists
    has_tax_table <- tryCatch({
      !is.null(phyloseq::tax_table(ps))
    }, error = function(e) {
      FALSE
    })
    
    if (!has_tax_table) return(NULL)
    
    # Detect if this is a large dataset
    large_dataset <- is_large_phyloseq()
    
    # Process with progress tracking for large datasets
    if (large_dataset) {
      withProgress(message = 'Building taxonomic tree...', value = 0, {
        tree_data <- build_taxonomic_tree(ps, threshold/100, setProgress)
        
        # Cache result if successful
        if (!is.null(tree_data)) {
          saveRDS(tree_data, cache_file)
        }
        
        return(tree_data)
      })
    } else {
      # Standard processing for smaller datasets
      tree_data <- build_taxonomic_tree(ps, threshold/100)
      
      # Cache result if successful
      if (!is.null(tree_data)) {
        saveRDS(tree_data, cache_file)
      }
      
      return(tree_data)
    }
  })
  
  # Helper function to build the taxonomy tree with optimization
  build_taxonomic_tree <- function(ps, threshold = 0.001, progress_callback = NULL) {
    # Update progress
    if (!is.null(progress_callback)) {
      progress_callback(0.1, "Creating relative abundance table")
    }
    
    # Create relative abundance (more efficiently for large datasets)
    ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
    
    # Get taxonomy data
    if (!is.null(progress_callback)) {
      progress_callback(0.2, "Extracting taxonomy data")
    }
    
    tax <- as.data.frame(tax_table(ps_rel))
    
    # Get taxa abundance (mean across samples)
    if (!is.null(progress_callback)) {
      progress_callback(0.3, "Calculating mean abundances")
    }
    
    taxa_sums_val <- phyloseq::taxa_sums(ps_rel) / phyloseq::nsamples(ps_rel)
    
    # Optimize for memory with large datasets
    is_large <- tryCatch({
      ntaxa(ps) > 5000 || nsamples(ps) > 100
    }, error = function(e) FALSE)
    
    # Only keep taxa above threshold
    if (!is.null(progress_callback)) {
      progress_callback(0.4, "Filtering by abundance threshold")
    }
    
    if (!is.null(threshold) && threshold > 0) {
      taxa_to_keep <- names(taxa_sums_val)[taxa_sums_val >= threshold]
      
      # Check if we have too many taxa after filtering
      if (length(taxa_to_keep) > 1000 && !is.null(progress_callback)) {
        progress_callback(0.45, paste("Large number of taxa:", length(taxa_to_keep), "- sub-sampling for visualization"))
        
        # For very large datasets, sample the taxa for better visualization
        if (length(taxa_to_keep) > 5000) {
          set.seed(123) # For reproducibility
          taxa_to_keep <- sample(taxa_to_keep, 5000)
        }
      }
      
      tax <- tax[rownames(tax) %in% taxa_to_keep, ]
      taxa_sums_val <- taxa_sums_val[taxa_to_keep]
    }
    
    # Use data.table for more efficient data handling with large datasets
    if (is_large && requireNamespace("data.table", quietly = TRUE)) {
      if (!is.null(progress_callback)) {
        progress_callback(0.5, "Using optimized data handling for large taxonomy")
      }
      
      # Convert to data.table for faster operations
      tax_dt <- as.data.table(tax)
      tax_dt[, ASV_ID := rownames(tax)]
      
      # Add abundance column
      tax_dt[, abundance := taxa_sums_val[ASV_ID]]
      
      # Create hierarchical paths more efficiently
      tax_dt[, Kingdom := as.character(Kingdom)]
      tax_dt[, Phylum := as.character(Phylum)]
      tax_dt[, Class := as.character(Class)]
      tax_dt[, Order := as.character(Order)]
      tax_dt[, Family := as.character(Family)]
      tax_dt[, Genus := as.character(Genus)]
      tax_dt[, Species := as.character(Species)]
      
      # Replace NAs with "Unknown"
      for (col in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
        tax_dt[is.na(get(col)), (col) := "Unknown"]
      }
      
      # Create hierarchical paths
      if (!is.null(progress_callback)) {
        progress_callback(0.6, "Creating taxonomy hierarchy")
      }
      
      tax_dt[, Phylum_path := paste(Kingdom, Phylum, sep = "|")]
      tax_dt[, Class_path := paste(Kingdom, Phylum, Class, sep = "|")]
      tax_dt[, Order_path := paste(Kingdom, Phylum, Class, Order, sep = "|")]
      tax_dt[, Family_path := paste(Kingdom, Phylum, Class, Order, Family, sep = "|")]
      tax_dt[, Genus_path := paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = "|")]
      tax_dt[, Species_path := paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "|")]
      
      # Melt to long format more efficiently
      if (!is.null(progress_callback)) {
        progress_callback(0.7, "Converting to long format")
      }
      
      # Create a list of all taxon paths with their abundances
      kingdom_dt <- unique(tax_dt[, .(Taxon = Kingdom, Level = "Kingdom", abundance = abundance)])
      phylum_dt <- unique(tax_dt[, .(Taxon = Phylum_path, Level = "Phylum", abundance = abundance)])
      class_dt <- unique(tax_dt[, .(Taxon = Class_path, Level = "Class", abundance = abundance)])
      order_dt <- unique(tax_dt[, .(Taxon = Order_path, Level = "Order", abundance = abundance)])
      family_dt <- unique(tax_dt[, .(Taxon = Family_path, Level = "Family", abundance = abundance)])
      genus_dt <- unique(tax_dt[, .(Taxon = Genus_path, Level = "Genus", abundance = abundance)])
      species_dt <- unique(tax_dt[, .(Taxon = Species_path, Level = "Species", abundance = abundance)])
      
      # Combine all levels
      tax_long <- rbindlist(list(kingdom_dt, phylum_dt, class_dt, order_dt, family_dt, genus_dt, species_dt))
      
      # Remove duplicates
      tax_long <- unique(tax_long)
      
    } else {
      # Standard processing for smaller datasets
      if (!is.null(progress_callback)) {
        progress_callback(0.5, "Creating taxonomy hierarchy")
      }
      
      # Combine taxonomy levels to create hierarchy
      tax_hierarchy <- data.frame(
        Kingdom = tax$Kingdom,
        Phylum = paste(tax$Kingdom, tax$Phylum, sep = "|"),
        Class = paste(tax$Kingdom, tax$Phylum, tax$Class, sep = "|"),
        Order = paste(tax$Kingdom, tax$Phylum, tax$Class, tax$Order, sep = "|"),
        Family = paste(tax$Kingdom, tax$Phylum, tax$Class, tax$Order, tax$Family, sep = "|"),
        Genus = paste(tax$Kingdom, tax$Phylum, tax$Class, tax$Order, tax$Family, tax$Genus, sep = "|"),
        Species = paste(tax$Kingdom, tax$Phylum, tax$Class, tax$Order, tax$Family, tax$Genus, tax$Species, sep = "|"),
        stringsAsFactors = FALSE
      )
      
      # Replace NAs with "Unknown"
      tax_hierarchy[] <- lapply(tax_hierarchy, function(x) gsub("NA", "Unknown", x))
      
      # Add abundance values as attributes
      tax_hierarchy$abundance <- taxa_sums_val[rownames(tax)]
      
      # Convert to long format for tree creation
      if (!is.null(progress_callback)) {
        progress_callback(0.7, "Converting to long format")
      }
      
      tax_long <- reshape2::melt(tax_hierarchy, 
                                id.vars = "abundance", 
                                variable.name = "Level", 
                                value.name = "Taxon")
      
      # Remove duplicated taxa (from different ASVs)
      tax_long <- tax_long[!duplicated(tax_long$Taxon), ]
    }
    
    # Function to create nested list for d3Tree with optimization for large datasets
    if (!is.null(progress_callback)) {
      progress_callback(0.8, "Building tree structure")
    }
    
    build_tree <- function(tax_data) {
      # Use a more efficient approach for large datasets
      if (nrow(tax_data) > 10000) {
        # Create a hash table for fast lookups
        taxon_lists <- new.env(hash = TRUE)
        
        # Get unique taxa
        unique_taxa <- unique(tax_data$Taxon)
        
        # First pass: create nodes
        for (taxon in unique_taxa) {
          # Extract the taxon parts
          taxon_parts <- strsplit(taxon, "\\|")[[1]]
          taxon_name <- taxon_parts[length(taxon_parts)]
          
          # Get abundance
          abund_rows <- which(tax_data$Taxon == taxon)
          abund <- if (length(abund_rows) > 0) tax_data$abundance[abund_rows[1]] else 0
          
          # Create list for this taxon
          taxon_lists[[taxon]] <- list(
            name = taxon_name,
            abundance = round(abund * 100, 2),  # Convert to percentage
            children = list()
          )
        }
        
        # Second pass: link parents and children
        for (taxon in unique_taxa) {
          # Skip root taxa
          if (!grepl("\\|", taxon)) next
          
          # Get parent
          parent_parts <- strsplit(taxon, "\\|")[[1]]
          parent <- paste(parent_parts[1:(length(parent_parts) - 1)], collapse = "|")
          
          # Add this taxon to parent's children if parent exists
          if (exists(parent, envir = taxon_lists, inherits = FALSE)) {
            child_node <- taxon_lists[[taxon]]
            parent_node <- taxon_lists[[parent]]
            parent_node$children <- c(parent_node$children, list(child_node))
            taxon_lists[[parent]] <- parent_node
          }
        }
        
        # Get root nodes
        root_nodes <- list()
        for (taxon in unique_taxa) {
          if (!grepl("\\|", taxon)) {
            root_nodes <- c(root_nodes, list(taxon_lists[[taxon]]))
          }
        }
        
      } else {
        # Standard approach for smaller datasets
        # Create empty lists for each unique taxon
        unique_taxa <- unique(tax_data$Taxon)
        taxon_lists <- setNames(vector("list", length(unique_taxa)), unique_taxa)
        
        # For each taxon, collect its children
        for (taxon in unique_taxa) {
          # Extract the taxon parts
          taxon_parts <- strsplit(taxon, "\\|")[[1]]
          taxon_name <- taxon_parts[length(taxon_parts)]
          
          # Get abundance for the taxon
          abund <- tax_data$abundance[tax_data$Taxon == taxon]
          if (length(abund) == 0) abund <- 0
          
          # Create list for this taxon
          taxon_lists[[taxon]] <- list(
            name = taxon_name,
            abundance = round(abund * 100, 2),  # Convert to percentage
            children = list()
          )
        }
        
        # Link parents and children
        for (taxon in unique_taxa) {
          # Skip the root taxa (Kingdom level)
          if (!grepl("\\|", taxon)) next
          
          # Get parent taxon
          parent_parts <- strsplit(taxon, "\\|")[[1]]
          parent <- paste(parent_parts[1:(length(parent_parts) - 1)], collapse = "|")
          
          # Add this taxon to its parent's children
          if (parent %in% names(taxon_lists)) {
            taxon_lists[[parent]]$children <- c(
              taxon_lists[[parent]]$children,
              list(taxon_lists[[taxon]])
            )
          }
        }
        
        # Return only the root (Kingdom) nodes
        root_nodes <- taxon_lists[!grepl("\\|", names(taxon_lists))]
      }
      
      # Force garbage collection
      gc()
      
      # Return tree structure
      return(list(name = "Taxonomy", children = root_nodes))
    }
    
    # Build tree structure
    tree_data <- build_tree(tax_long)
    
    if (!is.null(progress_callback)) {
      progress_callback(1.0, "Tree construction complete")
    }
    
    return(tree_data)
  }
  
  # Generate taxonomic tree on button click
  observeEvent(input$generateTree, {
    # Create the tree
    tree_data <- prepare_tree_data()
    
    if (is.null(tree_data)) {
      # Display message if no data
      output$taxTree <- renderUI({
        div(
          style = "text-align: center; padding-top: 100px;",
          h3("No taxonomy data available or all taxa below threshold.")
        )
      })
    } else {
      # Render the tree
      output$taxTree <- renderD3tree({
        d3tree(
          data = tree_data,
          width = "100%",
          height = 700,
          rootname = "Taxonomy",
          tooltip = TRUE,
          tooltipHtml = htmlwidgets::JS('function(d) { 
            return d.name + "<br>Abundance: " + d.abundance + "%";
          }')
        )
      })
    }
  })
  
  # Download handler for taxonomic tree
  output$downloadTree <- downloadHandler(
    filename = function() {
      paste0("taxonomy_tree_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
    },
    content = function(file) {
      tree_data <- prepare_tree_data()
      if (is.null(tree_data)) {
        # Create a basic HTML file with message
        writeLines('<html><body><h1>No taxonomy data available</h1></body></html>', file)
      } else {
        # Create the tree widget
        tree_widget <- d3tree(
          data = tree_data,
          width = "100%",
          height = 800,
          rootname = "Taxonomy",
          tooltip = TRUE,
          tooltipHtml = htmlwidgets::JS('function(d) { 
            return d.name + "<br>Abundance: " + d.abundance + "%";
          }')
        )
        
        # Save the widget as standalone HTML
        htmlwidgets::saveWidget(tree_widget, file, selfcontained = TRUE)
      }
    }
  )
  
  # Export tab functionality
  # Download handlers for tables
  output$downloadAsvTable <- downloadHandler(
    filename = function() {
      ext <- switch(input$exportAsvFormat,
                   "CSV" = "csv",
                   "TSV" = "tsv",
                   "Excel" = "xlsx")
      paste0("asv_table_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
    },
    content = function(file) {
      ps <- filtered_ps()
      if (is.null(ps)) return(NULL)
      
      # Extract ASV table
      asv_table <- as.data.frame(t(otu_table(ps)))
      
      # Try to add taxonomy if available
      tryCatch({
        if (!is.null(tax_table(ps))) {
          tax <- as.data.frame(tax_table(ps))
          
          # For each taxonomy level, add it to the ASV table if it exists
          for (level in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
            if (level %in% colnames(tax)) {
              # Match by row names
              asv_table[[level]] <- tax[rownames(asv_table), level]
            }
          }
        }
      }, error = function(e) {
        warning("Could not add taxonomy to ASV table:", conditionMessage(e))
      })
      
      # Write to file in the selected format
      if (input$exportAsvFormat == "CSV") {
        write.csv(asv_table, file, row.names = TRUE)
      } else if (input$exportAsvFormat == "TSV") {
        write.table(asv_table, file, sep = "\t", row.names = TRUE)
      } else if (input$exportAsvFormat == "Excel") {
        if (requireNamespace("openxlsx", quietly = TRUE)) {
          openxlsx::write.xlsx(asv_table, file, rowNames = TRUE)
        } else {
          # Fallback to CSV
          write.csv(asv_table, sub("\\.xlsx$", ".csv", file), row.names = TRUE)
        }
      }
    }
  )
  
  output$downloadTaxTable <- downloadHandler(
    filename = function() {
      ext <- switch(input$exportTaxFormat,
                   "CSV" = "csv",
                   "TSV" = "tsv",
                   "Excel" = "xlsx")
      paste0("taxonomy_table_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
    },
    content = function(file) {
      ps <- filtered_ps()
      if (is.null(ps)) return(NULL)
      
      # Check if taxonomy table exists
      if (is.null(tax_table(ps))) {
        # Create a dummy file with message
        write.csv(data.frame(Message = "No taxonomy data available"), file)
        return()
      }
      
      # Extract taxonomy table
      tax_table_df <- as.data.frame(tax_table(ps))
      
      # Add mean abundance information
      ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
      tax_table_df$MeanRelativeAbundance <- taxa_sums(ps_rel) / nsamples(ps_rel)
      
      # Write to file in the selected format
      if (input$exportTaxFormat == "CSV") {
        write.csv(tax_table_df, file, row.names = TRUE)
      } else if (input$exportTaxFormat == "TSV") {
        write.table(tax_table_df, file, sep = "\t", row.names = TRUE)
      } else if (input$exportTaxFormat == "Excel") {
        if (requireNamespace("openxlsx", quietly = TRUE)) {
          openxlsx::write.xlsx(tax_table_df, file, rowNames = TRUE)
        } else {
          # Fallback to CSV
          write.csv(tax_table_df, sub("\\.xlsx$", ".csv", file), row.names = TRUE)
        }
      }
    }
  )
  
  output$downloadMetadata <- downloadHandler(
    filename = function() {
      paste0("sample_metadata_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      ps <- filtered_ps()
      if (is.null(ps)) return(NULL)
      
      # Check if sample data exists
      has_sample_data <- tryCatch({
        !is.null(sample_data(ps))
      }, error = function(e) {
        FALSE
      })
      
      if (has_sample_data) {
        # Extract sample data
        metadata <- as.data.frame(sample_data(ps))
        # Add read counts
        metadata$TotalReads <- sample_sums(ps)
      } else {
        # Create basic sample data with read counts
        metadata <- data.frame(
          Sample = sample_names(ps),
          TotalReads = sample_sums(ps)
        )
      }
      
      # Check if we have unique SampleID column
      if (!"SampleID" %in% colnames(metadata)) {
        metadata$SampleID <- rownames(metadata) 
      }
      
      # Add standardized sample name for reference
      if (requireNamespace("stringr", quietly = TRUE)) {
        metadata$StandardizedName <- stringr::str_replace_all(
          metadata$SampleID, 
          pattern = "_R[12].*|\\.fastq.*|\\.fq.*", 
          replacement = ""
        )
      }
      
      # Write to CSV
      write.csv(metadata, file, row.names = FALSE)
    }
  )
  
  output$downloadReadTracking <- downloadHandler(
    filename = function() {
      paste0("read_tracking_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      track <- tracking_data()
      if (is.null(track)) {
        # Create a dummy file with message
        write.csv(data.frame(Message = "No read tracking data available"), file)
      } else {
        # Write tracking data to CSV
        write.csv(track, file, row.names = FALSE)
      }
    }
  )
  
  # Download handlers for plots
  output$downloadAlphaPlot <- downloadHandler(
    filename = function() {
      ext <- switch(input$exportPlotFormat,
                  "PNG" = "png",
                  "PDF" = "pdf",
                  "SVG" = "svg",
                  "HTML" = "html")
      paste0("alpha_diversity_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
    },
    content = function(file) {
      ps <- filtered_ps()
      if (is.null(ps) || is.null(input$alphaMethod)) {
        # Create a dummy plot
        p <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "No alpha diversity data available") +
          theme_void()
      } else {
        # Calculate alpha diversity
        alpha_div <- estimate_richness(ps, measures = input$alphaMethod)
        
        # Create data frame for plotting
        alpha_df <- data.frame(
          Sample = sample_names(ps),
          Diversity = alpha_div[[input$alphaMethod]]
        )
        
        # Create the plot
        p <- ggplot(alpha_df, aes(x = reorder(Sample, -Diversity), y = Diversity)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Sample", y = paste(input$alphaMethod, "Diversity"),
               title = paste(input$alphaMethod, "Alpha Diversity"))
        
        # Apply log scale if selected
        if (input$logScale) {
          p <- p + scale_y_log10()
        }
      }
      
      # Save in the selected format
      if (input$exportPlotFormat == "HTML") {
        # Convert to plotly and save as HTML
        p_ly <- ggplotly(p)
        htmlwidgets::saveWidget(p_ly, file, selfcontained = TRUE)
      } else {
        # Save as static image
        ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
      }
    }
  )
  
  output$downloadBetaPlot <- downloadHandler(
    filename = function() {
      ext <- switch(input$exportPlotFormat,
                  "PNG" = "png",
                  "PDF" = "pdf",
                  "SVG" = "svg",
                  "HTML" = "html")
      paste0("beta_diversity_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
    },
    content = function(file) {
      ps <- filtered_ps()
      if (is.null(ps) || is.null(input$betaMethod) || is.null(input$ordMethod)) {
        # Create a dummy plot
        p <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "No beta diversity data available") +
          theme_void()
      } else {
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
          p <- ggplot() + 
            annotate("text", x = 0.5, y = 0.5, 
                    label = "UniFrac requires a phylogenetic tree, which is not available") +
            theme_void()
        } else {
          # Calculate distance matrix and perform ordination
          dist_matrix <- phyloseq::distance(ps_rel, method = dist_method)
          
          ord_method <- switch(input$ordMethod,
                              "PCoA" = "PCoA",
                              "NMDS" = "NMDS",
                              "t-SNE" = "tsne",
                              "UMAP" = "umap")
          
          ord <- ordinate(ps_rel, method = ord_method, distance = dist_matrix)
          
          # Get ordination coordinates for plotting
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
          
          # Add sample names for plotting
          ord_data$Sample <- rownames(ord_data)
          
          # Create the plot
          p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, label = Sample)) +
            geom_point(size = 3, alpha = 0.7, color = "steelblue") +
            geom_text(nudge_y = 0.02, check_overlap = TRUE) +
            theme_minimal() +
            labs(x = axis_labels[1], y = axis_labels[2],
                title = paste(input$betaMethod, "Distance -", input$ordMethod, "Ordination"))
        }
      }
      
      # Save in the selected format
      if (input$exportPlotFormat == "HTML") {
        # Convert to plotly and save as HTML
        p_ly <- ggplotly(p)
        htmlwidgets::saveWidget(p_ly, file, selfcontained = TRUE)
      } else {
        # Save as static image
        ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
      }
    }
  )
  
  output$downloadTaxPlot <- downloadHandler(
    filename = function() {
      ext <- switch(input$exportPlotFormat,
                  "PNG" = "png",
                  "PDF" = "pdf",
                  "SVG" = "svg",
                  "HTML" = "html")
      paste0("taxonomy_barplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
    },
    content = function(file) {
      ps <- filtered_ps()
      if (is.null(ps) || is.null(input$taxLevel)) {
        # Create a dummy plot
        p <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "No taxonomy data available") +
          theme_void()
      } else {
        # Check for valid tax_table
        has_tax_table <- tryCatch({
          tt <- phyloseq::tax_table(ps)
          !is.null(tt) && input$taxLevel %in% colnames(tt)
        }, error = function(e) FALSE)
        
        if (!has_tax_table) {
          p <- ggplot() + 
            annotate("text", x = 0.5, y = 0.5, 
                    label = paste("No", input$taxLevel, "data available")) +
            theme_void()
        } else {
          # Agglomerate at selected taxonomic level
          ps_glom <- phyloseq::tax_glom(ps, taxrank = input$taxLevel)
          
          # Transform to relative abundance
          ps_rel <- phyloseq::transform_sample_counts(ps_glom, function(x) x / sum(x))
          
          # Create data frame for plotting
          tax_data <- phyloseq::psmelt(ps_rel)
          
          # Handle NA values in taxonomy
          tax_data[[input$taxLevel]] <- as.character(tax_data[[input$taxLevel]])
          tax_data[[input$taxLevel]][is.na(tax_data[[input$taxLevel]])] <- "Unknown"
          
          # Get top taxa
          n_top <- min(input$topTaxa, length(unique(tax_data$OTU)))
          taxa_sums_val <- phyloseq::taxa_sums(ps_rel)
          top_taxa <- names(sort(taxa_sums_val, decreasing = TRUE)[1:n_top])
          
          # Group remaining as Other if selected
          if (input$otherCat) {
            tax_data[[input$taxLevel]] <- ifelse(tax_data$OTU %in% top_taxa, 
                                              tax_data[[input$taxLevel]], "Other")
          } else {
            tax_data <- tax_data[tax_data$OTU %in% top_taxa, ]
          }
          
          # Create bar plot
          p <- ggplot(tax_data, aes(x = Sample, y = Abundance, fill = !!rlang::sym(input$taxLevel))) +
            geom_bar(stat = "identity") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = "Sample", y = "Relative Abundance", fill = input$taxLevel,
                title = paste("Taxonomy Composition at", input$taxLevel, "Level")) +
            scale_fill_viridis_d()
        }
      }
      
      # Save in the selected format
      if (input$exportPlotFormat == "HTML") {
        # Convert to plotly and save as HTML
        p_ly <- ggplotly(p)
        htmlwidgets::saveWidget(p_ly, file, selfcontained = TRUE)
      } else {
        # Save as static image
        ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
      }
    }
  )
  
  # Download handlers for R objects
  output$downloadPhyloseq <- downloadHandler(
    filename = function() {
      paste0("phyloseq_object_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      ps <- filtered_ps()
      if (is.null(ps)) {
        # Create a dummy phyloseq object
        dummy_ps <- phyloseq(otu_table(matrix(0, nrow=1, ncol=1), taxa_are_rows=TRUE))
        saveRDS(dummy_ps, file)
      } else {
        # Save the filtered phyloseq object
        saveRDS(ps, file)
      }
    }
  )
  
  # ===== MULTI-RUN ANALYSIS FUNCTIONALITY =====
  
  # Reactive value to store detected runs
  detected_runs <- reactiveVal(NULL)
  run_assignments <- reactiveVal(NULL)
  
  # Reactive to detect run identifiers
  detect_run_identifiers <- function() {
    ps <- phyloseq_obj()
    if (is.null(ps)) return(NULL)
    
    # Get sample names
    sample_ids <- sample_names(ps)
    run_ids <- NULL
    
    # Detection method based on user selection
    if (input$runIdentifier == "auto") {
      # Auto-detect based on sample name pattern
      pattern <- input$runPattern
      
      # Extract run information from sample names
      run_matches <- regmatches(sample_ids, regexpr(pattern, sample_ids, ignore.case = TRUE))
      
      # Count unique matches
      unique_runs <- unique(run_matches[run_matches != ""])
      
      if (length(unique_runs) >= 2) {
        # Create run mapping
        run_ids <- rep(NA, length(sample_ids))
        names(run_ids) <- sample_ids
        
        for (i in seq_along(sample_ids)) {
          match <- regmatches(sample_ids[i], regexpr(pattern, sample_ids[i], ignore.case = TRUE))
          if (length(match) > 0 && match != "") {
            run_ids[i] <- match
          }
        }
      }
    } else if (input$runIdentifier == "metadata") {
      # Use a metadata column
      if (input$runColumn != "auto") {
        # Check if we have sample data
        has_sample_data <- !is.null(sample_data(ps))
        
        if (has_sample_data) {
          # Extract sample data
          sample_df <- as.data.frame(sample_data(ps))
          
          # Use the specified column
          if (input$runColumn %in% colnames(sample_df)) {
            run_col_values <- sample_df[[input$runColumn]]
            names(run_col_values) <- rownames(sample_df)
            
            # Check if we have at least 2 different runs
            if (length(unique(run_col_values)) >= 2) {
              run_ids <- run_col_values
            }
          }
        }
      }
    }
    
    return(run_ids)
  }
  
  # Update run ID when the detect button is clicked
  observeEvent(input$detectRuns, {
    run_ids <- detect_run_identifiers()
    detected_runs(run_ids)
    
    # If we have run IDs, create a mapping of samples to runs
    if (!is.null(run_ids) && length(run_ids) > 0) {
      # Create run assignments for samples
      assignments <- data.frame(
        Sample = names(run_ids),
        Run = as.character(run_ids),
        stringsAsFactors = FALSE
      )
      
      # Remove rows with NA runs
      assignments <- assignments[!is.na(assignments$Run), ]
      
      if (nrow(assignments) > 0) {
        run_assignments(assignments)
      }
    }
  })
  
  # Display detected runs
  output$detectedRuns <- renderPrint({
    runs <- detected_runs()
    assignments <- run_assignments()
    
    if (is.null(runs) || length(runs) == 0) {
      cat("No runs detected. Please adjust the run identification method and try again.\n")
    } else {
      unique_runs <- unique(runs[!is.na(runs)])
      
      cat("Detected", length(unique_runs), "runs:\n")
      
      # Count samples per run
      run_counts <- table(runs)
      for (run in unique_runs) {
        cat("- ", run, ": ", run_counts[run], " samples\n", sep = "")
      }
      
      # Report samples without run assignment
      n_unassigned <- sum(is.na(runs))
      if (n_unassigned > 0) {
        cat("\n", n_unassigned, " samples could not be assigned to any run\n", sep = "")
      }
    }
  })
  
  # Update metadata column selection based on available sample data
  observe({
    ps <- phyloseq_obj()
    if (is.null(ps)) return()
    
    # Check if sample data exists
    has_sample_data <- tryCatch({
      !is.null(sample_data(ps))
    }, error = function(e) FALSE)
    
    if (has_sample_data) {
      # Extract sample data
      sample_df <- as.data.frame(sample_data(ps))
      
      # Get column names for potential run identifiers
      meta_cols <- colnames(sample_df)
      
      # Filter to columns that might contain run information
      run_cols <- meta_cols[grepl("run|batch|seq|plate", meta_cols, ignore.case = TRUE)]
      
      if (length(run_cols) > 0) {
        # Create choices with auto-detect first
        choices <- c("Auto-detect" = "auto", setNames(run_cols, run_cols))
        
        # Update the run column selector
        updateSelectInput(session, "runColumn", choices = choices)
        
        # Also update the metadata columns for coloring
        updateSelectInput(session, "runColorMetadata", choices = c("Auto-detect" = "auto", setNames(meta_cols, meta_cols)))
        
        # Update interaction factor options
        updateSelectInput(session, "interactionFactor", choices = c("Auto-detect" = "auto", setNames(meta_cols, meta_cols)))
      }
    }
  })
  
  # Generate checkboxes for run selection dynamically
  output$runSelectionCheckboxes <- renderUI({
    runs <- detected_runs()
    
    if (!is.null(runs)) {
      unique_runs <- unique(runs[!is.na(runs)])
      
      checkboxGroupInput("selectedRuns", "Select runs to include:",
                        choices = setNames(unique_runs, unique_runs),
                        selected = unique_runs)
    } else {
      helpText("No runs detected. Click 'Detect Runs' first.")
    }
  })
  
  # Run Quality Metrics plot
  output$runQualityMetrics <- renderPlotly({
    runs <- detected_runs()
    assignments <- run_assignments()
    
    if (is.null(runs) || is.null(assignments) || nrow(assignments) == 0) {
      # Return empty plot with message
      return(plot_ly() %>% 
              layout(title = "No run data available",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = "Click 'Detect Runs' to identify sequencing runs", 
                      showarrow = FALSE
                    )))
    }
    
    # Get phyloseq object
    ps <- filtered_ps()
    if (is.null(ps)) {
      return(plot_ly() %>% 
              layout(title = "No data available",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = "Phyloseq object not loaded", 
                      showarrow = FALSE
                    )))
    }
    
    # Calculate quality metrics per run
    tryCatch({
      # Get read counts for each sample
      read_counts <- sample_sums(ps)
      
      # Get alpha diversity for each sample
      alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon"))
      
      # Combine with run assignments
      metrics_df <- data.frame(
        Sample = sample_names(ps),
        ReadCount = read_counts,
        Observed = alpha_div$Observed,
        Shannon = alpha_div$Shannon,
        stringsAsFactors = FALSE
      )
      
      # Merge with run assignments
      metrics_df <- merge(metrics_df, assignments, by = "Sample", all.x = TRUE)
      
      # Replace NA runs with "Unassigned"
      metrics_df$Run[is.na(metrics_df$Run)] <- "Unassigned"
      
      # Create long format for plotting
      metrics_long <- reshape2::melt(metrics_df, 
                                   id.vars = c("Sample", "Run"),
                                   measure.vars = c("ReadCount", "Observed", "Shannon"),
                                   variable.name = "Metric", 
                                   value.name = "Value")
      
      # Create plot
      p <- ggplot(metrics_long, aes(x = Run, y = Value, fill = Run)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6, aes(text = Sample)) +
        facet_wrap(~Metric, scales = "free_y") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Quality Metrics by Run") +
        scale_fill_viridis_d()
      
      ggplotly(p, tooltip = c("text", "y"))
    }, error = function(e) {
      # Return error plot
      return(plot_ly() %>% 
              layout(title = "Error calculating run metrics",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = paste("Error:", conditionMessage(e)), 
                      showarrow = FALSE
                    )))
    })
  })
  
  # Run Beta Diversity Comparison
  output$runBetaDiversity <- renderPlotly({
    runs <- detected_runs()
    assignments <- run_assignments()
    
    if (is.null(runs) || is.null(assignments) || nrow(assignments) == 0) {
      # Return empty plot with message
      return(plot_ly() %>% 
              layout(title = "No run data available",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = "Click 'Detect Runs' to identify sequencing runs", 
                      showarrow = FALSE
                    )))
    }
    
    # Get phyloseq object
    ps <- filtered_ps()
    if (is.null(ps)) {
      return(plot_ly() %>% 
              layout(title = "No data available",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = "Phyloseq object not loaded", 
                      showarrow = FALSE
                    )))
    }
    
    # Calculate ordination
    tryCatch({
      # Transform to relative abundance
      ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
      
      # Calculate distance matrix
      dist_method <- input$runDistMethod
      dist_matrix <- phyloseq::distance(ps_rel, method = dist_method)
      
      # Perform ordination
      ord_method <- input$runCompareMethod
      
      if (ord_method == "pcoa") {
        ord <- ordinate(ps_rel, method = "PCoA", distance = dist_matrix)
        
        # Get ordination coordinates
        ord_data <- data.frame(ord$vectors[,1:2])
        colnames(ord_data) <- c("Axis1", "Axis2")
        var_explained <- round(ord$values$Relative_eig[1:2] * 100, 1)
        axis_labels <- paste0("Axis ", 1:2, " (", var_explained, "%)")
      } else if (ord_method == "nmds") {
        ord <- ordinate(ps_rel, method = "NMDS", distance = dist_matrix)
        ord_data <- data.frame(scores(ord)[,1:2])
        colnames(ord_data) <- c("Axis1", "Axis2")
        axis_labels <- c("NMDS1", "NMDS2")
      } else if (ord_method == "tsne") {
        # Try to use Rtsne package
        if (!requireNamespace("Rtsne", quietly = TRUE)) {
          return(plot_ly() %>% 
                  layout(title = "Rtsne package not available",
                        annotations = list(
                          x = 0.5, y = 0.5, 
                          text = "Please install Rtsne package: install.packages('Rtsne')", 
                          showarrow = FALSE
                        )))
        }
        
        # Convert distance matrix to Rtsne input format
        dist_mat <- as.matrix(dist_matrix)
        
        # Run t-SNE
        tsne_result <- Rtsne::Rtsne(dist_mat, is_distance = TRUE, perplexity = min(30, nrow(dist_mat) - 1))
        
        # Create data frame with the results
        ord_data <- data.frame(
          Axis1 = tsne_result$Y[,1],
          Axis2 = tsne_result$Y[,2]
        )
        rownames(ord_data) <- rownames(dist_mat)
        axis_labels <- c("t-SNE 1", "t-SNE 2")
      }
      
      # Add sample names
      ord_data$Sample <- rownames(ord_data)
      
      # Merge with run assignments
      ord_data <- merge(ord_data, assignments, by = "Sample", all.x = TRUE)
      
      # Replace NA runs with "Unassigned"
      ord_data$Run[is.na(ord_data$Run)] <- "Unassigned"
      
      # Create hover text
      ord_data$hover_text <- paste("Sample:", ord_data$Sample, "<br>Run:", ord_data$Run)
      
      # Add color based on selection
      if (input$runColorBy == "run") {
        # Color by run
        p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, color = Run, text = hover_text)) +
          geom_point(size = 3, alpha = 0.7) +
          theme_minimal() +
          labs(x = axis_labels[1], y = axis_labels[2],
               title = paste(toupper(ord_method), "Ordination of", dist_method, "Distances by Run")) +
          scale_color_viridis_d() +
          stat_ellipse(aes(color = Run), type = "t")
      } else {
        # Color by metadata
        meta_column <- input$runColorMetadata
        
        if (meta_column != "auto") {
          # Get sample data
          has_sample_data <- !is.null(sample_data(ps))
          
          if (has_sample_data) {
            # Extract sample data and merge
            sample_df <- as.data.frame(sample_data(ps))
            if (meta_column %in% colnames(sample_df)) {
              sample_df$Sample <- rownames(sample_df)
              ord_data <- merge(ord_data, sample_df[, c("Sample", meta_column)], by = "Sample", all.x = TRUE)
              
              # Create better hover text
              ord_data$hover_text <- paste("Sample:", ord_data$Sample, 
                                         "<br>Run:", ord_data$Run,
                                         "<br>", meta_column, ":", ord_data[[meta_column]])
              
              # Plot with metadata coloring
              p <- ggplot(ord_data, aes_string(x = "Axis1", y = "Axis2", color = meta_column, text = "hover_text")) +
                geom_point(size = 3, alpha = 0.7) +
                theme_minimal() +
                labs(x = axis_labels[1], y = axis_labels[2],
                     title = paste(toupper(ord_method), "Ordination of", dist_method, "Distances by", meta_column)) +
                scale_color_viridis_d()
              
              # Add run ellipses to still show run clustering
              p <- p + stat_ellipse(aes(color = Run, group = Run), type = "t", linetype = "dashed", alpha = 0.5)
            } else {
              # Fallback to run coloring
              p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, color = Run, text = hover_text)) +
                geom_point(size = 3, alpha = 0.7) +
                theme_minimal() +
                labs(x = axis_labels[1], y = axis_labels[2],
                     title = paste(toupper(ord_method), "Ordination of", dist_method, "Distances by Run")) +
                scale_color_viridis_d() +
                stat_ellipse(aes(color = Run), type = "t")
            }
          } else {
            # Fallback to run coloring
            p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, color = Run, text = hover_text)) +
              geom_point(size = 3, alpha = 0.7) +
              theme_minimal() +
              labs(x = axis_labels[1], y = axis_labels[2],
                   title = paste(toupper(ord_method), "Ordination of", dist_method, "Distances by Run")) +
              scale_color_viridis_d() +
              stat_ellipse(aes(color = Run), type = "t")
          }
        } else {
          # Auto-detect a good metadata column for coloring
          has_sample_data <- !is.null(sample_data(ps))
          
          if (has_sample_data) {
            # Extract sample data
            sample_df <- as.data.frame(sample_data(ps))
            
            # Look for good categorical columns
            potential_cols <- c("Group", "Treatment", "Condition", "SampleType", "SampleGroup")
            
            color_col <- NULL
            for (col in potential_cols) {
              if (col %in% colnames(sample_df)) {
                n_unique <- length(unique(sample_df[[col]]))
                if (n_unique > 1 && n_unique < 10) {
                  color_col <- col
                  break
                }
              }
            }
            
            if (!is.null(color_col)) {
              # Merge metadata
              sample_df$Sample <- rownames(sample_df)
              ord_data <- merge(ord_data, sample_df[, c("Sample", color_col)], by = "Sample", all.x = TRUE)
              
              # Create better hover text
              ord_data$hover_text <- paste("Sample:", ord_data$Sample, 
                                         "<br>Run:", ord_data$Run,
                                         "<br>", color_col, ":", ord_data[[color_col]])
              
              # Plot with metadata coloring
              p <- ggplot(ord_data, aes_string(x = "Axis1", y = "Axis2", color = color_col, text = "hover_text")) +
                geom_point(size = 3, alpha = 0.7) +
                theme_minimal() +
                labs(x = axis_labels[1], y = axis_labels[2],
                     title = paste(toupper(ord_method), "Ordination of", dist_method, "Distances by", color_col)) +
                scale_color_viridis_d()
              
              # Add run ellipses to still show run clustering
              p <- p + stat_ellipse(aes(color = Run, group = Run), type = "t", linetype = "dashed", alpha = 0.5)
            } else {
              # Fallback to run coloring
              p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, color = Run, text = hover_text)) +
                geom_point(size = 3, alpha = 0.7) +
                theme_minimal() +
                labs(x = axis_labels[1], y = axis_labels[2],
                     title = paste(toupper(ord_method), "Ordination of", dist_method, "Distances by Run")) +
                scale_color_viridis_d() +
                stat_ellipse(aes(color = Run), type = "t")
            }
          } else {
            # Fallback to run coloring
            p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, color = Run, text = hover_text)) +
              geom_point(size = 3, alpha = 0.7) +
              theme_minimal() +
              labs(x = axis_labels[1], y = axis_labels[2],
                   title = paste(toupper(ord_method), "Ordination of", dist_method, "Distances by Run")) +
              scale_color_viridis_d() +
              stat_ellipse(aes(color = Run), type = "t")
          }
        }
      }
      
      ggplotly(p, tooltip = "text")
    }, error = function(e) {
      # Return error plot
      return(plot_ly() %>% 
              layout(title = "Error calculating ordination",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = paste("Error:", conditionMessage(e)), 
                      showarrow = FALSE
                    )))
    })
  })
  
  # ASV Overlap Plot
  output$asvOverlapPlot <- renderPlot({
    runs <- detected_runs()
    assignments <- run_assignments()
    
    if (is.null(runs) || is.null(assignments) || nrow(assignments) == 0) {
      # Return empty plot with message
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                       label = "No run data available. Click 'Detect Runs' to identify runs.") +
               theme_void())
    }
    
    # Get phyloseq object
    ps <- filtered_ps()
    if (is.null(ps)) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                       label = "Phyloseq object not loaded") +
               theme_void())
    }
    
    # Filter runs if requested
    if (input$filterRunsForOverlap && !is.null(input$selectedRuns)) {
      # Filter to only include selected runs
      assignments <- assignments[assignments$Run %in% input$selectedRuns, ]
      
      if (nrow(assignments) == 0) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, 
                         label = "No samples match the selected runs") +
                 theme_void())
      }
    }
    
    # Generate the overlap visualization
    tryCatch({
      # Get OTU table
      otu <- otu_table(ps)
      
      if (taxa_are_rows(otu)) {
        otu_mat <- as.matrix(otu)
      } else {
        otu_mat <- t(as.matrix(otu))
      }
      
      # Create presence/absence matrix
      if (input$showPrevalentASVsOnly) {
        # Filter by prevalence
        threshold <- input$asvPrevalenceThreshold / 100
        prevalence <- rowSums(otu_mat > 0) / ncol(otu_mat)
        otu_mat <- otu_mat[prevalence >= threshold, ]
      }
      
      # Create binary presence/absence matrix
      binary_mat <- otu_mat > 0
      
      # Create a mapping of samples to runs
      sample_to_run <- assignments$Run
      names(sample_to_run) <- assignments$Sample
      
      # Get list of sample names in the matrix
      sample_names <- colnames(binary_mat)
      
      # Get run for each sample, defaulting to "Unassigned" for missing values
      sample_runs <- sample_to_run[sample_names]
      sample_runs[is.na(sample_runs)] <- "Unassigned"
      
      if (input$asvOverlapMethod == "upset") {
        # Check if UpSetR package is available
        if (!requireNamespace("UpSetR", quietly = TRUE)) {
          return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                           label = "UpSetR package not available. Please install: install.packages('UpSetR')") +
                   theme_void())
        }
        
        # Calculate run-specific ASV presence
        run_asv_presence <- list()
        unique_runs <- unique(sample_runs)
        
        for (run in unique_runs) {
          if (run != "Unassigned") {
            run_samples <- sample_names[sample_runs == run]
            if (length(run_samples) > 0) {
              run_asv_presence[[run]] <- which(rowSums(binary_mat[, run_samples, drop = FALSE]) > 0)
            }
          }
        }
        
        # Only include runs with detected ASVs
        run_asv_presence <- run_asv_presence[sapply(run_asv_presence, length) > 0]
        
        # Create UpSet plot
        if (length(run_asv_presence) >= 2) {
          UpSetR::upset(UpSetR::fromList(run_asv_presence), 
                      nsets = length(run_asv_presence),
                      order.by = "freq",
                      main.bar.color = "steelblue",
                      sets.bar.color = "darkgreen",
                      sets.x.label = "ASVs per Run",
                      mainbar.y.label = "Shared ASVs",
                      text.scale = 1.5)
        } else {
          return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                           label = "At least 2 runs with ASVs are required for UpSet plot") +
                   theme_void())
        }
        
      } else if (input$asvOverlapMethod == "venn") {
        # Check if we have a reasonable number of runs for Venn diagram
        unique_runs <- unique(sample_runs)
        unique_runs <- unique_runs[unique_runs != "Unassigned"]
        
        if (length(unique_runs) < 2) {
          return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                           label = "At least 2 runs are required for Venn diagram") +
                   theme_void())
        }
        
        if (length(unique_runs) > 5) {
          return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                           label = "Venn diagrams only support up to 5 sets. Use UpSet plot for more runs.") +
                   theme_void())
        }
        
        # Check if VennDiagram package is available
        if (!requireNamespace("VennDiagram", quietly = TRUE)) {
          return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                           label = "VennDiagram package not available. Please install: install.packages('VennDiagram')") +
                   theme_void())
        }
        
        # Calculate run-specific ASV presence
        run_asv_presence <- list()
        
        for (run in unique_runs) {
          run_samples <- sample_names[sample_runs == run]
          if (length(run_samples) > 0) {
            run_asv_presence[[run]] <- which(rowSums(binary_mat[, run_samples, drop = FALSE]) > 0)
          }
        }
        
        # Get a nice color palette
        run_colors <- viridis::viridis(length(unique_runs))
        
        # Generate Venn diagram
        venn_plot <- VennDiagram::venn.diagram(
          x = run_asv_presence,
          filename = NULL,
          fill = run_colors,
          alpha = 0.5,
          cex = 2,
          cat.cex = 1.5,
          cat.col = run_colors,
          main = "ASV Overlap Between Runs",
          main.cex = 2
        )
        
        # Convert to a grob (graphics object)
        grid::grid.newpage()
        grid::grid.draw(venn_plot)
        
      } else if (input$asvOverlapMethod == "heatmap") {
        # Calculate ASV presence by run
        unique_runs <- unique(sample_runs)
        
        # Calculate the presence/absence of each ASV in each run
        asv_run_matrix <- matrix(0, nrow = nrow(binary_mat), ncol = length(unique_runs))
        rownames(asv_run_matrix) <- rownames(binary_mat)
        colnames(asv_run_matrix) <- unique_runs
        
        for (i in seq_along(unique_runs)) {
          run <- unique_runs[i]
          run_samples <- sample_names[sample_runs == run]
          if (length(run_samples) > 0) {
            asv_run_matrix[, i] <- rowSums(binary_mat[, run_samples, drop = FALSE]) > 0
          }
        }
        
        # Convert to a data frame for ggplot
        asv_run_df <- reshape2::melt(asv_run_matrix, 
                                   varnames = c("ASV", "Run"), 
                                   value.name = "Present")
        
        # Add taxonomy information if available
        has_taxonomy <- !is.null(tax_table(ps))
        if (has_taxonomy) {
          tax_df <- as.data.frame(tax_table(ps))
          
          # Add taxonomy columns to the data frame
          asv_run_df$Phylum <- tax_df[asv_run_df$ASV, "Phylum"]
          asv_run_df$Class <- tax_df[asv_run_df$ASV, "Class"]
          asv_run_df$Order <- tax_df[asv_run_df$ASV, "Order"]
          asv_run_df$Family <- tax_df[asv_run_df$ASV, "Family"]
          asv_run_df$Genus <- tax_df[asv_run_df$ASV, "Genus"]
          
          # Handle NAs
          for (col in c("Phylum", "Class", "Order", "Family", "Genus")) {
            asv_run_df[[col]][is.na(asv_run_df[[col]])] <- "Unknown"
          }
          
          # Create taxonomy label
          asv_run_df$TaxLabel <- paste(asv_run_df$Phylum, asv_run_df$Family, asv_run_df$Genus, sep = "; ")
        } else {
          asv_run_df$TaxLabel <- asv_run_df$ASV
        }
        
        # Create a more informative ASV label with taxonomy
        asv_run_df$ASVLabel <- paste0(asv_run_df$ASV, " (", asv_run_df$TaxLabel, ")")
        
        # Create the heatmap
        p <- ggplot(asv_run_df, aes(x = Run, y = ASVLabel, fill = Present)) +
          geom_tile() +
          scale_fill_manual(values = c("white", "steelblue"), 
                          labels = c("Absent", "Present")) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Run", y = "ASV", fill = "Presence",
               title = "ASV Presence/Absence Across Runs")
        
        p
      }
    }, error = function(e) {
      # Return error plot
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                       label = paste("Error:", conditionMessage(e))) +
               theme_void())
    })
  })
  
  # Taxonomy by Run Plot
  output$taxonomyRunPlot <- renderPlotly({
    runs <- detected_runs()
    assignments <- run_assignments()
    
    if (is.null(runs) || is.null(assignments) || nrow(assignments) == 0) {
      # Return empty plot with message
      return(plot_ly() %>% 
              layout(title = "No run data available",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = "Click 'Detect Runs' to identify sequencing runs", 
                      showarrow = FALSE
                    )))
    }
    
    # Get phyloseq object
    ps <- filtered_ps()
    if (is.null(ps)) {
      return(plot_ly() %>% 
              layout(title = "No data available",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = "Phyloseq object not loaded", 
                      showarrow = FALSE
                    )))
    }
    
    # Check if taxonomy table exists
    has_tax_table <- tryCatch({
      tt <- phyloseq::tax_table(ps)
      !is.null(tt) && ncol(tt) > 0 && nrow(tt) > 0
    }, error = function(e) FALSE)
    
    if (!has_tax_table || !input$taxRunLevel %in% colnames(phyloseq::tax_table(ps))) {
      return(plot_ly() %>% 
              layout(title = paste("No", input$taxRunLevel, "data available"),
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = paste("Taxonomy table does not contain", input$taxRunLevel), 
                      showarrow = FALSE
                    )))
    }
    
    # Generate the taxonomy by run visualization
    tryCatch({
      # Get run for each sample
      sample_to_run <- assignments$Run
      names(sample_to_run) <- assignments$Sample
      
      # Create a new sample_data with run information
      sample_df <- as.data.frame(sample_data(ps))
      sample_df$Run <- sample_to_run[rownames(sample_df)]
      sample_df$Run[is.na(sample_df$Run)] <- "Unassigned"
      
      # Update the phyloseq object with run information
      sample_data(ps) <- sample_data(sample_df)
      
      # Agglomerate at the selected taxonomic level
      ps_glom <- tax_glom(ps, taxrank = input$taxRunLevel)
      
      # Transform to relative abundance
      ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
      
      # Get the top N taxa
      top_n <- input$taxRunTopN
      top_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE)[1:top_n])
      
      # Subset to top taxa
      ps_rel_top <- prune_taxa(top_taxa, ps_rel)
      
      if (input$taxRunVizType == "bar") {
        # Melt for plotting
        tax_data <- psmelt(ps_rel_top)
        
        # Aggregate by run and taxa
        tax_by_run <- tax_data %>%
          dplyr::group_by(Run, .data[[input$taxRunLevel]]) %>%
          dplyr::summarize(MeanAbundance = mean(Abundance), .groups = "drop")
        
        # Handle NA values
        tax_by_run[[input$taxRunLevel]] <- as.character(tax_by_run[[input$taxRunLevel]])
        tax_by_run[[input$taxRunLevel]][is.na(tax_by_run[[input$taxRunLevel]])] <- "Unknown"
        
        # Create stacked bar plot
        p <- ggplot(tax_by_run, aes_string(x = "Run", y = "MeanAbundance", fill = input$taxRunLevel)) +
          geom_bar(stat = "identity") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Run", y = "Mean Relative Abundance", fill = input$taxRunLevel,
               title = paste("Mean Taxonomic Composition at", input$taxRunLevel, "Level by Run")) +
          scale_fill_viridis_d()
        
        ggplotly(p)
        
      } else if (input$taxRunVizType == "heatmap") {
        # Melt for analysis
        tax_data <- psmelt(ps_rel_top)
        
        # Aggregate by run and taxa
        tax_by_run <- tax_data %>%
          dplyr::group_by(Run, .data[[input$taxRunLevel]]) %>%
          dplyr::summarize(MeanAbundance = mean(Abundance), .groups = "drop")
        
        # Handle NA values
        tax_by_run[[input$taxRunLevel]] <- as.character(tax_by_run[[input$taxRunLevel]])
        tax_by_run[[input$taxRunLevel]][is.na(tax_by_run[[input$taxRunLevel]])] <- "Unknown"
        
        # Reshape to wide format for heatmap
        tax_wide <- reshape2::dcast(tax_by_run, Run ~ .data[[input$taxRunLevel]], value.var = "MeanAbundance")
        row.names(tax_wide) <- tax_wide$Run
        tax_wide$Run <- NULL
        
        # Handle clustering if requested
        if (input$taxRunClusterMethod != "none") {
          # Cluster the runs
          run_dist <- dist(tax_wide)
          run_hclust <- hclust(run_dist, method = input$taxRunClusterMethod)
          run_order <- row.names(tax_wide)[run_hclust$order]
          
          # Cluster the taxa
          taxa_dist <- dist(t(tax_wide))
          taxa_hclust <- hclust(taxa_dist, method = input$taxRunClusterMethod)
          taxa_order <- colnames(tax_wide)[taxa_hclust$order]
          
          # Create heatmap with ordered rows and columns
          tax_wide_ordered <- tax_wide[run_order, taxa_order]
          
          # Convert to long format for ggplot
          tax_long <- reshape2::melt(as.matrix(tax_wide_ordered), 
                                   varnames = c("Run", input$taxRunLevel), 
                                   value.name = "Abundance")
          
          # Make sure Run and taxa are ordered factors for ggplot
          tax_long$Run <- factor(tax_long$Run, levels = run_order)
          tax_long[[input$taxRunLevel]] <- factor(tax_long[[input$taxRunLevel]], levels = taxa_order)
        } else {
          # Convert to long format without clustering
          tax_long <- reshape2::melt(as.matrix(tax_wide),
                                   varnames = c("Run", input$taxRunLevel),
                                   value.name = "Abundance")
        }
        
        # Create heatmap
        p <- ggplot(tax_long, aes_string(x = input$taxRunLevel, y = "Run", fill = "Abundance")) +
          geom_tile() +
          scale_fill_viridis_c() +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = input$taxRunLevel, y = "Run", fill = "Abundance",
               title = paste("Heatmap of", input$taxRunLevel, "Abundance by Run"))
        
        ggplotly(p)
        
      } else if (input$taxRunVizType == "pcoa") {
        # Melt for analysis
        tax_data <- psmelt(ps_rel_top)
        
        # Aggregate by run and taxa
        tax_by_run <- tax_data %>%
          dplyr::group_by(Run, .data[[input$taxRunLevel]]) %>%
          dplyr::summarize(MeanAbundance = mean(Abundance), .groups = "drop")
        
        # Handle NA values
        tax_by_run[[input$taxRunLevel]] <- as.character(tax_by_run[[input$taxRunLevel]])
        tax_by_run[[input$taxRunLevel]][is.na(tax_by_run[[input$taxRunLevel]])] <- "Unknown"
        
        # Reshape to wide format for ordination
        tax_wide <- reshape2::dcast(tax_by_run, Run ~ .data[[input$taxRunLevel]], value.var = "MeanAbundance")
        row.names(tax_wide) <- tax_wide$Run
        tax_wide$Run <- NULL
        
        # Replace NAs with zeros
        tax_wide[is.na(tax_wide)] <- 0
        
        # Calculate distance matrix
        dist_mat <- vegdist(tax_wide, method = "bray")
        
        # Perform PCoA ordination
        pcoa_result <- ape::pcoa(dist_mat)
        
        # Create data frame for plotting
        pcoa_df <- data.frame(
          Run = rownames(tax_wide),
          Axis1 = pcoa_result$vectors[,1],
          Axis2 = pcoa_result$vectors[,2]
        )
        
        # Calculate variance explained
        var_explained <- round(pcoa_result$values$Relative_eig[1:2] * 100, 1)
        
        # Create PCoA plot
        p <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, label = Run, color = Run)) +
          geom_point(size = 4) +
          geom_text(hjust = 0, vjust = 0, nudge_x = 0.01, nudge_y = 0.01) +
          theme_minimal() +
          labs(x = paste0("PCoA 1 (", var_explained[1], "%)"), 
               y = paste0("PCoA 2 (", var_explained[2], "%)"),
               title = paste("PCoA of", input$taxRunLevel, "Composition by Run")) +
          scale_color_viridis_d()
        
        ggplotly(p)
      }
    }, error = function(e) {
      # Return error plot
      return(plot_ly() %>% 
              layout(title = "Error in taxonomy visualization",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = paste("Error:", conditionMessage(e)), 
                      showarrow = FALSE
                    )))
    })
  })
  
  # Run Batch Effect Analysis
  # Store batch effect test results
  batch_effect_results <- reactiveVal(NULL)
  
  # Run batch effect analysis when button is clicked
  observeEvent(input$runBatchTest, {
    runs <- detected_runs()
    assignments <- run_assignments()
    
    if (is.null(runs) || is.null(assignments) || nrow(assignments) == 0) {
      batch_effect_results(data.frame(
        Test = input$batchEffectMethod,
        Result = "No run data available. Click 'Detect Runs' to identify runs.",
        p_value = NA,
        statistic = NA
      ))
      return()
    }
    
    # Get phyloseq object
    ps <- filtered_ps()
    if (is.null(ps)) {
      batch_effect_results(data.frame(
        Test = input$batchEffectMethod,
        Result = "Phyloseq object not loaded",
        p_value = NA,
        statistic = NA
      ))
      return()
    }
    
    # Run the selected batch effect analysis
    tryCatch({
      # Get run for each sample
      sample_to_run <- assignments$Run
      names(sample_to_run) <- assignments$Sample
      
      # Create a new sample_data with run information
      sample_df <- as.data.frame(sample_data(ps))
      sample_df$Run <- sample_to_run[rownames(sample_df)]
      sample_df$Run[is.na(sample_df$Run)] <- "Unassigned"
      
      # Update the phyloseq object with run information
      sample_data(ps) <- sample_data(sample_df)
      
      # Transform to relative abundance
      ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
      
      # Calculate distance matrix
      dist_method <- input$batchDistMethod
      dist_matrix <- phyloseq::distance(ps_rel, method = dist_method)
      
      # Run the appropriate batch effect test
      if (input$batchEffectMethod == "betadisper") {
        # Check if we have more than one run group
        run_groups <- unique(sample_df$Run)
        if (length(run_groups) < 2) {
          batch_effect_results(data.frame(
            Test = "Beta Dispersion",
            Result = "At least 2 runs required for beta dispersion test",
            p_value = NA,
            statistic = NA
          ))
          return()
        }
        
        # Run beta dispersion test (test of homogeneity of multivariate dispersions)
        bd <- vegan::betadisper(dist_matrix, sample_df$Run)
        anova_result <- anova(bd)
        perm_result <- vegan::permutest(bd, permutations = 999)
        
        # Store results
        batch_effect_results(data.frame(
          Test = "Beta Dispersion",
          Result = ifelse(perm_result$tab$`Pr(>F)`[1] < 0.05, 
                         "Significant dispersion differences between runs (p < 0.05)", 
                         "No significant dispersion differences between runs (p ≥ 0.05)"),
          p_value = perm_result$tab$`Pr(>F)`[1],
          statistic = perm_result$tab$F[1]
        ))
        
        # Store BD object for plotting
        bd_object <- bd
        
      } else if (input$batchEffectMethod == "adonis") {
        # Run PERMANOVA test
        adonis_result <- vegan::adonis2(dist_matrix ~ Run, data = sample_df)
        
        # Store results
        batch_effect_results(data.frame(
          Test = "PERMANOVA (adonis2)",
          Result = ifelse(adonis_result$`Pr(>F)`[1] < 0.05, 
                         "Significant differences in community composition between runs (p < 0.05)", 
                         "No significant differences in community composition between runs (p ≥ 0.05)"),
          p_value = adonis_result$`Pr(>F)`[1],
          statistic = adonis_result$F[1]
        ))
        
      } else if (input$batchEffectMethod == "interaction") {
        # Run vs. Sample Type interaction test
        interaction_factor <- input$interactionFactor
        
        if (interaction_factor == "auto") {
          # Try to find a good factor automatically
          potential_cols <- c("Group", "Treatment", "Condition", "SampleType", "SampleGroup")
          
          # Find the first column that works
          for (col in potential_cols) {
            if (col %in% colnames(sample_df)) {
              unique_vals <- unique(sample_df[[col]])
              if (!is.numeric(unique_vals) && length(unique_vals) > 1 && length(unique_vals) < 10) {
                interaction_factor <- col
                break
              }
            }
          }
        }
        
        if (interaction_factor != "auto" && interaction_factor %in% colnames(sample_df)) {
          # Run PERMANOVA with interaction term
          adonis_formula <- as.formula(paste("dist_matrix ~", interaction_factor, "* Run"))
          adonis_result <- vegan::adonis2(adonis_formula, data = sample_df)
          
          # Extract interaction term
          interaction_row <- which(rownames(adonis_result) == paste0(interaction_factor, ":Run"))
          
          if (length(interaction_row) > 0) {
            # Store results focusing on the interaction term
            batch_effect_results(data.frame(
              Test = paste("Interaction:", interaction_factor, "x Run"),
              Result = ifelse(adonis_result$`Pr(>F)`[interaction_row] < 0.05, 
                             paste0("Significant interaction between ", interaction_factor, " and Run (p < 0.05)\n",
                                  "This indicates that run effects differ across ", interaction_factor, " groups."), 
                             paste0("No significant interaction between ", interaction_factor, " and Run (p ≥ 0.05)\n",
                                  "This suggests batch effects are consistent across ", interaction_factor, " groups.")),
              p_value = adonis_result$`Pr(>F)`[interaction_row],
              statistic = adonis_result$F[interaction_row]
            ))
          } else {
            # If interaction term wasn't found, store the main effect of Run
            run_row <- which(rownames(adonis_result) == "Run")
            
            batch_effect_results(data.frame(
              Test = paste("Main effect of Run with", interaction_factor),
              Result = ifelse(adonis_result$`Pr(>F)`[run_row] < 0.05, 
                             "Significant effect of Run after accounting for sample type (p < 0.05)", 
                             "No significant effect of Run after accounting for sample type (p ≥ 0.05)"),
              p_value = adonis_result$`Pr(>F)`[run_row],
              statistic = adonis_result$F[run_row]
            ))
          }
        } else {
          # No suitable factor found
          batch_effect_results(data.frame(
            Test = "Interaction test",
            Result = "No suitable sample type factor found for interaction test",
            p_value = NA,
            statistic = NA
          ))
        }
      } else if (input$batchEffectMethod == "rda") {
        # Redundancy Analysis (RDA)
        # Check if vegan package is available
        if (!requireNamespace("vegan", quietly = TRUE)) {
          batch_effect_results(data.frame(
            Test = "RDA",
            Result = "Vegan package is required but not available",
            p_value = NA,
            statistic = NA
          ))
          return()
        }
        
        # Convert distance matrix to an ordination object
        # We'll need the OTU table for RDA
        otu_table_matrix <- as.matrix(otu_table(ps_rel))
        if (!taxa_are_rows(ps_rel)) {
          otu_table_matrix <- t(otu_table_matrix)
        }
        
        # Run RDA with Run as constraining variable
        rda_result <- vegan::rda(otu_table_matrix ~ Run, data = sample_df)
        
        # Calculate significance with permutation test
        rda_test <- vegan::anova.cca(rda_result, permutations = 999)
        
        # Extract proportion of variance explained by Run factor
        variance_explained <- vegan::RsquareAdj(rda_result)$r.squared
        
        # Store results
        batch_effect_results(data.frame(
          Test = "RDA (Redundancy Analysis)",
          Result = ifelse(rda_test$`Pr(>F)`[1] < 0.05, 
                         paste0("Significant batch effect detected (p < 0.05)\n",
                              "Run explains ", round(variance_explained * 100, 1), "% of community variation."), 
                         paste0("No significant batch effect detected (p ≥ 0.05)\n",
                              "Run explains only ", round(variance_explained * 100, 1), "% of community variation.")),
          p_value = rda_test$`Pr(>F)`[1],
          statistic = rda_test$F[1]
        ))
        
      } else if (input$batchEffectMethod == "mantel") {
        # Mantel test to compare run distance matrix with community distance
        # Create a run distance matrix (0 for same run, 1 for different)
        run_dist <- matrix(0, nrow = nrow(sample_df), ncol = nrow(sample_df))
        rownames(run_dist) <- rownames(sample_df)
        colnames(run_dist) <- rownames(sample_df)
        
        # Fill run distance matrix
        for (i in 1:(nrow(sample_df)-1)) {
          for (j in (i+1):nrow(sample_df)) {
            if (sample_df$Run[i] != sample_df$Run[j]) {
              run_dist[i, j] <- run_dist[j, i] <- 1
            }
          }
        }
        
        # Convert to distance object
        run_dist_obj <- as.dist(run_dist)
        
        # Run Mantel test
        mantel_result <- vegan::mantel(dist_matrix, run_dist_obj, method = "spearman", permutations = 999)
        
        # Store results
        batch_effect_results(data.frame(
          Test = "Mantel Test",
          Result = ifelse(mantel_result$signif < 0.05, 
                         paste0("Significant correlation between run membership and community composition (p < 0.05)\n",
                              "Mantel r = ", round(mantel_result$statistic, 3), " indicates batch effects are present."), 
                         paste0("No significant correlation between run membership and community composition (p ≥ 0.05)\n",
                              "Mantel r = ", round(mantel_result$statistic, 3), " suggests minimal batch effects.")),
          p_value = mantel_result$signif,
          statistic = mantel_result$statistic
        ))
        
      } else if (input$batchEffectMethod == "anosim") {
        # Analysis of Similarities (ANOSIM)
        anosim_result <- vegan::anosim(dist_matrix, sample_df$Run, permutations = 999)
        
        # Store results
        batch_effect_results(data.frame(
          Test = "ANOSIM",
          Result = ifelse(anosim_result$signif < 0.05, 
                         paste0("Significant differences between runs (p < 0.05)\n",
                              "ANOSIM R = ", round(anosim_result$statistic, 3), " indicates strong batch effects."), 
                         paste0("No significant differences between runs (p ≥ 0.05)\n",
                              "ANOSIM R = ", round(anosim_result$statistic, 3), " suggests minimal batch effects.")),
          p_value = anosim_result$signif,
          statistic = anosim_result$statistic
        ))
      }
      
    }, error = function(e) {
      # Store error result
      batch_effect_results(data.frame(
        Test = input$batchEffectMethod,
        Result = paste("Error:", conditionMessage(e)),
        p_value = NA,
        statistic = NA
      ))
    })
  })
  
  # Display batch effect test results
  output$batchTestResults <- renderPrint({
    results <- batch_effect_results()
    
    if (is.null(results)) {
      cat("No analysis run yet. Click 'Run Analysis' to perform batch effect testing.\n")
    } else {
      cat("BATCH EFFECT ANALYSIS:", results$Test, "\n\n")
      cat("RESULT:", results$Result, "\n\n")
      
      if (!is.na(results$p_value)) {
        cat("p-value:", format(results$p_value, digits = 4), "\n")
      }
      
      if (!is.na(results$statistic)) {
        if (results$Test == "Beta Dispersion") {
          cat("F-statistic:", format(results$statistic, digits = 4), "\n")
        } else if (grepl("PERMANOVA", results$Test)) {
          cat("pseudo F-statistic:", format(results$statistic, digits = 4), "\n")
        }
      }
    }
  })
  
  # Batch effect visualization
  output$batchEffectPlot <- renderPlotly({
    results <- batch_effect_results()
    
    if (is.null(results)) {
      # Return empty plot with message
      return(plot_ly() %>% 
              layout(title = "No analysis run yet",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = "Click 'Run Analysis' to perform batch effect testing", 
                      showarrow = FALSE
                    )))
    }
    
    # Visualization depends on the test used
    if (results$Test == "Beta Dispersion") {
      # Get phyloseq object
      ps <- filtered_ps()
      if (is.null(ps)) return(NULL)
      
      # Get run assignments
      runs <- detected_runs()
      assignments <- run_assignments()
      
      if (is.null(runs) || is.null(assignments)) return(NULL)
      
      # Create a PCoA plot with ellipses showing dispersion
      tryCatch({
        # Get run for each sample
        sample_to_run <- assignments$Run
        names(sample_to_run) <- assignments$Sample
        
        # Create a new sample_data with run information
        sample_df <- as.data.frame(sample_data(ps))
        sample_df$Run <- sample_to_run[rownames(sample_df)]
        sample_df$Run[is.na(sample_df$Run)] <- "Unassigned"
        
        # Update the phyloseq object with run information
        sample_data(ps) <- sample_data(sample_df)
        
        # Transform to relative abundance
        ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
        
        # Calculate distance matrix
        dist_method <- input$batchDistMethod
        dist_matrix <- phyloseq::distance(ps_rel, method = dist_method)
        
        # Run beta dispersion analysis
        bd <- vegan::betadisper(dist_matrix, sample_df$Run)
        
        # Extract PCoA scores
        bd_data <- data.frame(
          Sample = rownames(bd$vectors),
          PC1 = bd$vectors[,1],
          PC2 = bd$vectors[,2],
          Run = sample_df$Run
        )
        
        # Get centroids
        centroids <- data.frame(
          Run = rownames(bd$centroids),
          PC1 = bd$centroids[,1],
          PC2 = bd$centroids[,2]
        )
        
        # Create PCoA plot with dispersion visualization
        p <- ggplot(bd_data, aes(x = PC1, y = PC2, color = Run, text = Sample)) +
          geom_point(size = 3, alpha = 0.7) +
          stat_ellipse(aes(group = Run, color = Run), type = "t") +
          geom_point(data = centroids, aes(x = PC1, y = PC2, color = Run), 
                    size = 5, shape = 8) +
          theme_minimal() +
          labs(x = paste0("PCoA 1"), y = paste0("PCoA 2"),
               title = "Beta Dispersion Analysis - Multivariate Homogeneity of Group Dispersions") +
          scale_color_viridis_d()
        
        ggplotly(p, tooltip = c("text", "color"))
        
      }, error = function(e) {
        # Return error plot
        return(plot_ly() %>% 
                layout(title = "Error in beta dispersion plot",
                      annotations = list(
                        x = 0.5, y = 0.5, 
                        text = paste("Error:", conditionMessage(e)), 
                        showarrow = FALSE
                      )))
      })
    } else if (grepl("PERMANOVA|adonis", results$Test) || grepl("Main effect", results$Test)) {
      # For PERMANOVA, show NMDS plot
      ps <- filtered_ps()
      if (is.null(ps)) return(NULL)
      
      # Get run assignments
      runs <- detected_runs()
      assignments <- run_assignments()
      
      if (is.null(runs) || is.null(assignments)) return(NULL)
      
      # Create a NMDS/PCoA plot
      tryCatch({
        # Get run for each sample
        sample_to_run <- assignments$Run
        names(sample_to_run) <- assignments$Sample
        
        # Create a new sample_data with run information
        sample_df <- as.data.frame(sample_data(ps))
        sample_df$Run <- sample_to_run[rownames(sample_df)]
        sample_df$Run[is.na(sample_df$Run)] <- "Unassigned"
        
        # Update the phyloseq object with run information
        sample_data(ps) <- sample_data(sample_df)
        
        # Transform to relative abundance
        ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
        
        # Calculate distance matrix
        dist_method <- input$batchDistMethod
        dist_matrix <- phyloseq::distance(ps_rel, method = dist_method)
        
        # Perform PCoA ordination (usually more stable than NMDS)
        ord <- ordinate(ps_rel, method = "PCoA", distance = dist_matrix)
        
        # Create data frame for plotting
        ord_data <- data.frame(
          Sample = rownames(ord$vectors),
          Axis1 = ord$vectors[,1],
          Axis2 = ord$vectors[,2],
          Run = sample_df$Run
        )
        
        # Get variance explained
        var_explained <- round(ord$values$Relative_eig[1:2] * 100, 1)
        
        # Create PCoA plot
        p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, color = Run, text = Sample)) +
          geom_point(size = 3, alpha = 0.7) +
          stat_ellipse(aes(group = Run), type = "t") +
          theme_minimal() +
          labs(x = paste0("PCoA 1 (", var_explained[1], "%)"), 
               y = paste0("PCoA 2 (", var_explained[2], "%)"),
               title = paste("PCoA showing batch differences (PERMANOVA:", 
                            format(results$p_value, digits = 3), ")")) +
          scale_color_viridis_d()
        
        ggplotly(p, tooltip = c("text", "color"))
        
      }, error = function(e) {
        # Return error plot
        return(plot_ly() %>% 
                layout(title = "Error in PERMANOVA plot",
                      annotations = list(
                        x = 0.5, y = 0.5, 
                        text = paste("Error:", conditionMessage(e)), 
                        showarrow = FALSE
                      )))
      })
    } else if (grepl("Interaction", results$Test)) {
      # For interaction test, create a grouped PCoA plot
      ps <- filtered_ps()
      if (is.null(ps)) return(NULL)
      
      # Get run assignments
      runs <- detected_runs()
      assignments <- run_assignments()
      
      if (is.null(runs) || is.null(assignments)) return(NULL)
      
      # Extract the interaction factor
      interaction_factor <- gsub("Interaction: ([^ ]+) x Run.*", "\\1", results$Test)
      
      # Create plot showing interaction
      tryCatch({
        # Get run for each sample
        sample_to_run <- assignments$Run
        names(sample_to_run) <- assignments$Sample
        
        # Create a new sample_data with run information
        sample_df <- as.data.frame(sample_data(ps))
        sample_df$Run <- sample_to_run[rownames(sample_df)]
        sample_df$Run[is.na(sample_df$Run)] <- "Unassigned"
        
        # Check if we have the interaction factor
        if (interaction_factor != "Run" && interaction_factor %in% colnames(sample_df)) {
          # Update the phyloseq object with run information
          sample_data(ps) <- sample_data(sample_df)
          
          # Transform to relative abundance
          ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
          
          # Calculate distance matrix
          dist_method <- input$batchDistMethod
          dist_matrix <- phyloseq::distance(ps_rel, method = dist_method)
          
          # Perform PCoA ordination
          ord <- ordinate(ps_rel, method = "PCoA", distance = dist_matrix)
          
          # Create data frame for plotting
          ord_data <- data.frame(
            Sample = rownames(ord$vectors),
            Axis1 = ord$vectors[,1],
            Axis2 = ord$vectors[,2],
            Run = sample_df$Run
          )
          
          # Add the interaction factor
          ord_data[[interaction_factor]] <- sample_df[ord_data$Sample, interaction_factor]
          
          # Get variance explained
          var_explained <- round(ord$values$Relative_eig[1:2] * 100, 1)
          
          # Create a faceted plot to show the interaction
          p <- ggplot(ord_data, aes(x = Axis1, y = Axis2, color = Run, 
                                    text = paste("Sample:", Sample, 
                                               "<br>Run:", Run, 
                                               "<br>", interaction_factor, ":", ord_data[[interaction_factor]]))) +
            geom_point(size = 3, alpha = 0.7) +
            facet_wrap(as.formula(paste("~", interaction_factor))) +
            theme_minimal() +
            labs(x = paste0("PCoA 1 (", var_explained[1], "%)"), 
                 y = paste0("PCoA 2 (", var_explained[2], "%)"),
                 title = paste("PCoA by", interaction_factor, "and Run (interaction p =", 
                              format(results$p_value, digits = 3), ")")) +
            scale_color_viridis_d()
          
          ggplotly(p, tooltip = "text")
        } else {
          # Fallback to showing a message
          return(plot_ly() %>% 
                  layout(title = paste("Could not find", interaction_factor, "in sample data"),
                        annotations = list(
                          x = 0.5, y = 0.5, 
                          text = "Interaction factor not available for visualization", 
                          showarrow = FALSE
                        )))
        }
      }, error = function(e) {
        # Return error plot
        return(plot_ly() %>% 
                layout(title = "Error in interaction plot",
                      annotations = list(
                        x = 0.5, y = 0.5, 
                        text = paste("Error:", conditionMessage(e)), 
                        showarrow = FALSE
                      )))
      })
    } else {
      # Default visualization for other tests
      return(plot_ly() %>% 
              layout(title = results$Test,
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = results$Result, 
                      showarrow = FALSE
                    )))
    }
  })
  
  # Normalization-related functionality
  normalized_ps <- reactiveVal(NULL)
  
  # Compute normalization when requested
  observeEvent(input$computeNormalization, {
    ps <- filtered_ps()
    if (is.null(ps)) {
      return()
    }
    
    # Apply the selected normalization method
    tryCatch({
      # Start with relative abundance transformation as baseline
      ps_norm <- transform_sample_counts(ps, function(x) x / sum(x))
      method_desc <- "Relative Abundance"
      
      if (input$runNormMethod == "css") {
        # CSS normalization (requires metagenomeSeq)
        if (requireNamespace("metagenomeSeq", quietly = TRUE)) {
          # Convert to metagenomeSeq format
          otu_mat <- as.matrix(otu_table(ps))
          if (!taxa_are_rows(otu_table(ps))) {
            otu_mat <- t(otu_mat)
          }
          
          # Create metagenomeSeq object
          mgs_obj <- metagenomeSeq::newMRexperiment(otu_mat)
          
          # Apply CSS normalization
          mgs_norm <- metagenomeSeq::cumNorm(mgs_obj)
          
          # Get normalized counts
          css_counts <- metagenomeSeq::MRcounts(mgs_norm, normalized = TRUE)
          
          # Create new otu_table
          otu_norm <- otu_table(css_counts, taxa_are_rows = TRUE)
          
          # Replace the otu_table in the phyloseq object
          ps_norm <- phyloseq(otu_norm, 
                             tax_table(ps), 
                             sample_data(ps), 
                             phy_tree(ps))
          
          method_desc <- "CSS (Cumulative Sum Scaling)"
        } else {
          method_desc <- "CSS (failed - metagenomeSeq package not available)"
        }
      } else if (input$runNormMethod == "rle") {
        # RLE normalization (requires DESeq2)
        if (requireNamespace("DESeq2", quietly = TRUE)) {
          # Convert to DESeq2 format - needs to handle zeros
          # Get OTU table
          otu_mat <- as.matrix(otu_table(ps))
          if (!taxa_are_rows(otu_table(ps))) {
            otu_mat <- t(otu_mat)
          }
          
          # Add pseudocount to handle zeros
          otu_mat_pseudo <- otu_mat + 1
          
          # Create DESeq2 object
          sample_df <- as.data.frame(sample_data(ps))
          
          # Need a design formula - use a column from sample data or create one
          design_col <- "group"
          if (!"group" %in% colnames(sample_df)) {
            # Create a dummy grouping variable
            sample_df$group <- factor(rep("A", nrow(sample_df)))
          }
          
          # Create DESeq2 object
          dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = otu_mat_pseudo,
            colData = sample_df,
            design = as.formula(paste("~", design_col))
          )
          
          # Estimate size factors
          dds <- DESeq2::estimateSizeFactors(dds)
          
          # Get normalized counts
          norm_counts <- DESeq2::counts(dds, normalized = TRUE)
          
          # Create new otu_table
          otu_norm <- otu_table(norm_counts, taxa_are_rows = TRUE)
          
          # Replace the otu_table in the phyloseq object
          ps_norm <- phyloseq(otu_norm, 
                             tax_table(ps), 
                             sample_data(ps), 
                             phy_tree(ps))
          
          method_desc <- "RLE (Relative Log Expression from DESeq2)"
        } else {
          method_desc <- "RLE (failed - DESeq2 package not available)"
        }
      } else if (input$runNormMethod == "tmm") {
        # TMM normalization (requires edgeR)
        if (requireNamespace("edgeR", quietly = TRUE)) {
          # Get OTU table
          otu_mat <- as.matrix(otu_table(ps))
          if (!taxa_are_rows(otu_table(ps))) {
            otu_mat <- t(otu_mat)
          }
          
          # Create DGEList object
          dge <- edgeR::DGEList(counts = otu_mat)
          
          # Calculate normalization factors
          dge <- edgeR::calcNormFactors(dge, method = "TMM")
          
          # Get normalized counts
          norm_counts <- edgeR::cpm(dge)
          
          # Create new otu_table
          otu_norm <- otu_table(norm_counts, taxa_are_rows = TRUE)
          
          # Replace the otu_table in the phyloseq object
          ps_norm <- phyloseq(otu_norm, 
                             tax_table(ps), 
                             sample_data(ps), 
                             phy_tree(ps))
          
          method_desc <- "TMM (Trimmed Mean of M-values from edgeR)"
        } else {
          method_desc <- "TMM (failed - edgeR package not available)"
        }
      } else if (input$runNormMethod == "clr") {
        # CLR transformation (requires compositions package)
        if (requireNamespace("compositions", quietly = TRUE)) {
          # Get OTU table
          otu_mat <- as.matrix(otu_table(ps))
          if (!taxa_are_rows(otu_table(ps))) {
            otu_mat <- t(otu_mat)
          }
          
          # Add small pseudocount to handle zeros
          otu_mat_pseudo <- otu_mat + 0.5
          
          # Apply CLR transformation
          clr_counts <- t(compositions::clr(t(otu_mat_pseudo)))
          
          # Create new otu_table
          otu_norm <- otu_table(clr_counts, taxa_are_rows = TRUE)
          
          # Replace the otu_table in the phyloseq object
          ps_norm <- phyloseq(otu_norm, 
                             tax_table(ps), 
                             sample_data(ps), 
                             phy_tree(ps))
          
          method_desc <- "CLR (Centered Log-Ratio Transformation)"
        } else {
          method_desc <- "CLR (failed - compositions package not available)"
        }
      } else if (input$runNormMethod == "vst") {
        # Variance Stabilizing Transformation (requires DESeq2)
        if (requireNamespace("DESeq2", quietly = TRUE)) {
          # Get OTU table
          otu_mat <- as.matrix(otu_table(ps))
          if (!taxa_are_rows(otu_table(ps))) {
            otu_mat <- t(otu_mat)
          }
          
          # Round to integers (VST requires count data)
          otu_mat <- round(otu_mat)
          
          # Create DESeq2 object
          sample_df <- as.data.frame(sample_data(ps))
          
          # Need a design formula - use a column from sample data or create one
          design_col <- "group"
          if (!"group" %in% colnames(sample_df)) {
            # Create a dummy grouping variable
            sample_df$group <- factor(rep("A", nrow(sample_df)))
          }
          
          # Create DESeq2 object
          dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = otu_mat,
            colData = sample_df,
            design = as.formula(paste("~", design_col))
          )
          
          # Apply variance stabilizing transformation
          vst_data <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
          vst_counts <- DESeq2::assay(vst_data)
          
          # Create new otu_table
          otu_norm <- otu_table(vst_counts, taxa_are_rows = TRUE)
          
          # Replace the otu_table in the phyloseq object
          ps_norm <- phyloseq(otu_norm, 
                             tax_table(ps), 
                             sample_data(ps), 
                             phy_tree(ps))
          
          method_desc <- "VST (Variance Stabilizing Transformation)"
        } else {
          method_desc <- "VST (failed - DESeq2 package not available)"
        }
      } else if (input$runNormMethod == "alr") {
        # ALR transformation (requires compositions package)
        if (requireNamespace("compositions", quietly = TRUE)) {
          # Get OTU table
          otu_mat <- as.matrix(otu_table(ps))
          if (!taxa_are_rows(otu_table(ps))) {
            otu_mat <- t(otu_mat)
          }
          
          # Add pseudocount to handle zeros
          otu_mat_pseudo <- otu_mat + 0.5
          
          # Apply ALR transformation - using the last taxon as denominator
          alr_counts <- t(compositions::alr(t(otu_mat_pseudo)))
          
          # Create new otu_table
          otu_norm <- otu_table(alr_counts, taxa_are_rows = TRUE)
          
          # Replace the otu_table in the phyloseq object
          ps_norm <- phyloseq(otu_norm, 
                             tax_table(ps), 
                             sample_data(ps), 
                             phy_tree(ps))
          
          method_desc <- "ALR (Additive Log-Ratio Transformation)"
        } else {
          method_desc <- "ALR (failed - compositions package not available)"
        }
      }
      
      # Store the normalized phyloseq object
      normalized_ps(list(ps = ps_norm, method = method_desc))
      
    }, error = function(e) {
      # Return error information
      normalized_ps(list(
        ps = NULL, 
        method = paste("Error:", conditionMessage(e))
      ))
    })
  })
  
  # Display normalization summary
  output$normalizationSummary <- renderPrint({
    norm <- normalized_ps()
    
    if (is.null(norm)) {
      if (input$runNormMethod == "none" || input$runNormMethod == "relative") {
        cat("Using", ifelse(input$runNormMethod == "none", "raw counts", "relative abundance"), "for analysis.\n")
        cat("This is applied automatically to all analyses.\n")
      } else {
        cat("Normalization not yet computed. Click 'Compute Normalization' to apply the selected method.\n")
      }
    } else {
      cat("Normalization method:", norm$method, "\n\n")
      
      if (!is.null(norm$ps)) {
        # Print summary statistics for the normalized data
        cat("Summary of normalized abundance values:\n")
        otu <- otu_table(norm$ps)
        otu_mat <- as.matrix(otu)
        
        cat("Min:", format(min(otu_mat), digits = 4), "\n")
        cat("Max:", format(max(otu_mat), digits = 4), "\n")
        cat("Mean:", format(mean(otu_mat), digits = 4), "\n")
        cat("Median:", format(median(otu_mat), digits = 4), "\n")
        
        # Calculate coefficient of variation for raw vs. normalized data
        ps_orig <- filtered_ps_orig()
        if (!is.null(ps_orig)) {
          otu_raw <- as.matrix(otu_table(ps_orig))
          # Only compare if matrices have same dimensions
          if (dim(otu_raw)[1] == dim(otu_mat)[1] && dim(otu_raw)[2] == dim(otu_mat)[2]) {
            # Calculate CV for each sample
            sample_cv_raw <- apply(otu_raw, 2, function(x) sd(x) / mean(x) * 100)
            sample_cv_norm <- apply(otu_mat, 2, function(x) sd(x) / mean(x) * 100)
            
            cat("\nCoefficient of variation (CV):\n")
            cat("- Raw data mean CV: ", format(mean(sample_cv_raw, na.rm = TRUE), digits = 4), "%\n")
            cat("- Normalized data mean CV: ", format(mean(sample_cv_norm, na.rm = TRUE), digits = 4), "%\n")
            
            if (mean(sample_cv_norm, na.rm = TRUE) < mean(sample_cv_raw, na.rm = TRUE)) {
              cat("✓ Normalization reduced variability.\n")
            } else {
              cat("Note: This normalization method may not reduce variability for this dataset.\n")
            }
          }
        }
        
        if (input$applyRunNormalization) {
          cat("\nThis normalization is now being applied to all analyses in the dashboard.\n")
        } else {
          cat("\nTo use this normalization in other analyses, check 'Apply normalization to all analyses'.\n")
        }
      }
    }
  })
  
  # Normalization comparison plot
  output$normalizationComparisonPlot <- renderPlotly({
    norm <- normalized_ps()
    ps <- filtered_ps()
    
    if (is.null(ps)) {
      return(plot_ly() %>% 
              layout(title = "No data available",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = "Phyloseq object not loaded", 
                      showarrow = FALSE
                    )))
    }
    
    # Get run assignments for coloring
    runs <- detected_runs()
    assignments <- run_assignments()
    
    # Create comparison visualization
    tryCatch({
      # Create a comparison plot showing raw, relative, and normalized abundances
      # Start with raw counts
      raw_sums <- sample_sums(ps)
      raw_df <- data.frame(
        Sample = names(raw_sums),
        Value = raw_sums,
        Type = "Raw counts"
      )
      
      # Add relative abundance
      rel_ps <- transform_sample_counts(ps, function(x) x / sum(x) * mean(raw_sums))
      rel_sums <- sample_sums(rel_ps)
      rel_df <- data.frame(
        Sample = names(rel_sums),
        Value = rel_sums,
        Type = "Relative abundance (scaled)"
      )
      
      # Create combined data frame
      plot_df <- rbind(raw_df, rel_df)
      
      # Add normalized counts if available
      if (!is.null(norm) && !is.null(norm$ps)) {
        # Get normalized counts and scale to same range
        norm_sums <- sample_sums(norm$ps)
        
        # Handle zero/negative values that might result from transformations like CLR
        if (any(norm_sums <= 0)) {
          # For transformations that can yield negative values (like CLR)
          # We'll shift values to make them positive for visualization
          if (min(norm_sums) <= 0) {
            norm_sums <- norm_sums + abs(min(norm_sums)) + 1
          }
        }
        
        scaling_factor <- mean(raw_sums) / mean(norm_sums)
        norm_sums_scaled <- norm_sums * scaling_factor
        
        norm_df <- data.frame(
          Sample = names(norm_sums),
          Value = norm_sums_scaled,
          Type = paste(norm$method, "(scaled)")
        )
        
        # Add to plot data
        plot_df <- rbind(plot_df, norm_df)
      }
      
      # Calculate coefficient of variation for each normalization method
      # This shows how the methods affect the variability across samples
      cv_data <- plot_df %>%
        dplyr::group_by(Type) %>%
        dplyr::summarize(CoeffVar = sd(Value) / mean(Value) * 100) %>%
        dplyr::arrange(CoeffVar)
      
      # Add run information if available
      if (!is.null(assignments) && nrow(assignments) > 0) {
        # Create sample to run mapping
        sample_to_run <- assignments$Run
        names(sample_to_run) <- assignments$Sample
        
        # Add run column to plotting data
        plot_df$Run <- sample_to_run[plot_df$Sample]
        plot_df$Run[is.na(plot_df$Run)] <- "Unassigned"
        
        # Create comparison box plots by run
        p <- ggplot(plot_df, aes(x = Run, y = Value, fill = Type, text = Sample)) +
          geom_boxplot(position = position_dodge(width = 0.8)) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Run", y = "Scaled Read Count", fill = "Normalization",
               title = paste0("Comparison of Normalization Methods Across Runs\n",
                             "Coefficient of Variation: ", 
                             paste(cv_data$Type, sprintf("(%.1f%%)", cv_data$CoeffVar), collapse = ", "))) +
          scale_y_continuous(labels = scales::comma)
      } else {
        # Create simple boxplot comparison without runs
        p <- ggplot(plot_df, aes(x = Type, y = Value, fill = Type, text = Sample)) +
          geom_boxplot() +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Normalization Method", y = "Scaled Read Count", fill = "Normalization",
               title = paste0("Comparison of Normalization Methods\n",
                             "Coefficient of Variation: ", 
                             paste(cv_data$Type, sprintf("(%.1f%%)", cv_data$CoeffVar), collapse = ", "))) +
          scale_y_continuous(labels = scales::comma)
      }
      
      ggplotly(p, tooltip = c("text", "y", "fill"))
      
    }, error = function(e) {
      # Return error plot
      return(plot_ly() %>% 
              layout(title = "Error in normalization comparison",
                    annotations = list(
                      x = 0.5, y = 0.5, 
                      text = paste("Error:", conditionMessage(e)), 
                      showarrow = FALSE
                    )))
    })
  })
  
  # Override filtered_ps to use normalized version if requested
  filtered_ps_orig <- filtered_ps
  filtered_ps <- reactive({
    # Check if normalization should be applied
    if (input$applyRunNormalization && input$runNormMethod != "none") {
      norm <- normalized_ps()
      
      if (!is.null(norm) && !is.null(norm$ps)) {
        # Return the normalized version
        return(norm$ps)
      }
      
      # If normalization is requested but not computed, apply relative abundance
      if (input$runNormMethod == "relative") {
        ps <- filtered_ps_orig()
        if (!is.null(ps)) {
          return(transform_sample_counts(ps, function(x) x / sum(x)))
        }
      }
    }
    
    # Return the original filtered phyloseq
    return(filtered_ps_orig())
  })
  
  output$downloadBiom <- downloadHandler(
    filename = function() {
      paste0("biom_export_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".biom")
    },
    content = function(file) {
      ps <- filtered_ps()
      if (is.null(ps)) {
        # Create a text file with message
        writeLines("No data available to export to BIOM format", file)
      } else {
        # Try to convert to BIOM format
        tryCatch({
          if (requireNamespace("biomformat", quietly = TRUE)) {
            # Convert phyloseq to BIOM
            biom_object <- biomformat::make_biom(
              data = as(otu_table(ps), "matrix"),
              sample_metadata = as(sample_data(ps), "data.frame"),
              observation_metadata = as(tax_table(ps), "matrix")
            )
            
            # Write BIOM file
            biomformat::write_biom(biom_object, file)
          } else {
            # Create a text file with message
            writeLines("Package 'biomformat' is required but not installed.\nInstall with: install.packages('biomformat')", file)
          }
        }, error = function(e) {
          writeLines(paste("Error exporting to BIOM format:", conditionMessage(e)), file)
        })
      }
    }
  )
  
  # Download handler for comprehensive report
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste0("dada2_analysis_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
    },
    content = function(file) {
      # Create a temporary report file
      tempReport <- tempfile(fileext = ".Rmd")
      
      # Generate report content based on selected sections
      report_content <- c(
        "---",
        "title: \"DADA2 Microbiome Analysis Report\"",
        "date: \"`r format(Sys.time(), '%Y-%m-%d %H:%M:%S')`\"",
        "output: ",
        "  html_document:",
        "    toc: true",
        "    toc_float: true",
        "    theme: cosmo",
        "---",
        "",
        "```{r setup, include=FALSE}",
        "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.width = 10, fig.height = 7)",
        "library(phyloseq)",
        "library(ggplot2)",
        "library(dplyr)",
        "library(vegan)",
        "library(viridis)",
        "```",
        "",
        "```{r load_data}",
        "# Load the data from the current session",
        "ps <- readRDS(\"results/phyloseq_object.rds\")",
        "track <- NULL",
        "try({",
        "  track <- readRDS(\"results/read_tracking.rds\")",
        "}, silent = TRUE)",
        "```",
        ""
      )
      
      # Add selected sections
      if ("Overview" %in% input$reportSections) {
        report_content <- c(report_content,
          "# Overview",
          "",
          "```{r overview}",
          "# Summary statistics",
          "cat(\"**Number of samples:** \", phyloseq::nsamples(ps), \"\\n\\n\")",
          "cat(\"**Number of ASVs:** \", phyloseq::ntaxa(ps), \"\\n\\n\")",
          "cat(\"**Total reads:** \", format(sum(phyloseq::sample_sums(ps)), big.mark = \",\"), \"\\n\\n\")",
          "cat(\"**Mean reads per sample:** \", format(round(mean(phyloseq::sample_sums(ps))), big.mark = \",\"), \"\\n\\n\")",
          "cat(\"**Median reads per sample:** \", format(median(phyloseq::sample_sums(ps)), big.mark = \",\"), \"\\n\\n\")",
          "```",
          "",
          "## Sample Read Counts",
          "",
          "```{r sample_reads}",
          "# Create data frame with sample sums",
          "df <- data.frame(",
          "  Sample = sample_names(ps),",
          "  Reads = sample_sums(ps)",
          ")",
          "",
          "# Create plot",
          "ggplot(df, aes(x = reorder(Sample, -Reads), y = Reads)) +",
          "  geom_bar(stat = \"identity\", fill = \"steelblue\") +",
          "  theme_minimal() +",
          "  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +",
          "  labs(x = \"Sample\", y = \"Read Count\") +",
          "  scale_y_continuous(labels = scales::comma)",
          "```",
          "",
          "## Read Processing Stats",
          "",
          "```{r read_tracking, eval=!is.null(track)}",
          "# Check if tracking data is available",
          "if (!is.null(track) && is.data.frame(track)) {",
          "  # Identify numeric columns",
          "  numeric_cols <- sapply(track, is.numeric)",
          "  id_vars <- c(\"Sample\")",
          "  measure_vars <- setdiff(names(track)[numeric_cols], id_vars)",
          "  ",
          "  if (length(measure_vars) > 0) {",
          "    # Reshape data for plotting",
          "    track_long <- reshape2::melt(track, id.vars = \"Sample\", ",
          "                               measure.vars = measure_vars,",
          "                               variable.name = \"Step\", value.name = \"Reads\")",
          "    ",
          "    # Create plot",
          "    ggplot(track_long, aes(x = Step, y = Reads, group = Sample, color = Sample)) +",
          "      geom_line() +",
          "      geom_point() +",
          "      theme_minimal() +",
          "      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +",
          "      labs(y = \"Read Count\") +",
          "      scale_y_continuous(labels = scales::comma)",
          "  } else {",
          "    cat(\"No numeric columns found in read tracking data.\")",
          "  }",
          "} else {",
          "  cat(\"No read tracking data available.\")",
          "}",
          "```",
          "",
          "## Taxonomic Overview",
          "",
          "```{r taxa_overview}",
          "# Check if taxonomy table exists",
          "has_tax <- !is.null(tax_table(ps)) && \"Phylum\" %in% colnames(tax_table(ps))",
          "",
          "if (has_tax) {",
          "  # Agglomerate at phylum level",
          "  ps_phylum <- phyloseq::tax_glom(ps, taxrank = \"Phylum\")",
          "  ",
          "  # Transform to relative abundance",
          "  ps_phylum_rel <- phyloseq::transform_sample_counts(ps_phylum, function(x) x / sum(x))",
          "  ",
          "  # Create data frame for plotting",
          "  phylum_data <- phyloseq::psmelt(ps_phylum_rel)",
          "  ",
          "  # Handle NA values in Phylum",
          "  phylum_data$Phylum <- as.character(phylum_data$Phylum)",
          "  phylum_data$Phylum[is.na(phylum_data$Phylum)] <- \"Unknown\"",
          "  ",
          "  # Get top phyla (up to 10)",
          "  taxa_sums_val <- phyloseq::taxa_sums(ps_phylum_rel)",
          "  n_phyla <- min(10, length(unique(phylum_data$Phylum)))",
          "  top_phyla <- names(sort(taxa_sums_val, decreasing = TRUE)[1:n_phyla])",
          "  ",
          "  # Group remaining as Other",
          "  phylum_data$Phylum <- ifelse(phylum_data$OTU %in% top_phyla, ",
          "                              phylum_data$Phylum, \"Other\")",
          "  ",
          "  # Create plot",
          "  ggplot(phylum_data, aes(x = Sample, y = Abundance, fill = Phylum)) +",
          "    geom_bar(stat = \"identity\") +",
          "    theme_minimal() +",
          "    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +",
          "    labs(x = \"Sample\", y = \"Relative Abundance\") +",
          "    scale_fill_viridis_d()",
          "} else {",
          "  cat(\"No taxonomy data available at Phylum level.\")",
          "}",
          "```",
          ""
        )
      }
      
      if ("Quality Control" %in% input$reportSections) {
        report_content <- c(report_content,
          "# Quality Control",
          "",
          "## Read Count Distribution",
          "",
          "```{r read_dist}",
          "# Create data frame with read counts",
          "read_counts <- data.frame(",
          "  Reads = sample_sums(ps)",
          ")",
          "",
          "# Create histogram",
          "ggplot(read_counts, aes(x = Reads)) +",
          "  geom_histogram(bins = 30, fill = \"steelblue\", color = \"black\") +",
          "  theme_minimal() +",
          "  labs(x = \"Read Count\", y = \"Number of Samples\") +",
          "  scale_x_continuous(labels = scales::comma)",
          "```",
          "",
          "## ASV Length Distribution",
          "",
          "```{r asv_lengths}",
          "# Get ASV sequences",
          "asv_seqs <- colnames(otu_table(ps))",
          "",
          "# Calculate sequence length distribution",
          "seq_lengths <- nchar(asv_seqs)",
          "",
          "# Create data frame for plotting",
          "length_df <- data.frame(table(seq_lengths))",
          "length_df$seq_lengths <- as.numeric(as.character(length_df$seq_lengths))",
          "",
          "# Create plot",
          "ggplot(length_df, aes(x = seq_lengths, y = Freq)) +",
          "  geom_bar(stat = \"identity\", fill = \"steelblue\") +",
          "  theme_minimal() +",
          "  labs(x = \"Sequence Length\", y = \"Count\") +",
          "  scale_x_continuous(breaks = sort(unique(seq_lengths)))",
          "```",
          "",
          "## Sequence Length Statistics",
          "",
          "```{r seq_stats}",
          "cat(\"**Total ASVs:** \", length(asv_seqs), \"\\n\\n\")",
          "cat(\"**Min length:** \", min(seq_lengths), \"\\n\\n\")",
          "cat(\"**Max length:** \", max(seq_lengths), \"\\n\\n\")",
          "cat(\"**Mean length:** \", round(mean(seq_lengths), 1), \"\\n\\n\")",
          "cat(\"**Median length:** \", median(seq_lengths), \"\\n\\n\")",
          "```",
          ""
        )
      }
      
      if ("Alpha Diversity" %in% input$reportSections) {
        report_content <- c(report_content,
          "# Alpha Diversity",
          "",
          "Alpha diversity measures the diversity within individual samples.",
          "",
          "## Diversity Metrics",
          "",
          "```{r alpha_div}",
          "# Calculate alpha diversity",
          "alpha_div <- estimate_richness(ps, measures = c(\"Observed\", \"Shannon\", \"Simpson\"))",
          "",
          "# Prepare data for plotting",
          "alpha_long <- reshape2::melt(alpha_div, var.name = \"Metric\", value.name = \"Value\")",
          "alpha_long$Sample <- rep(sample_names(ps), each = 3)",
          "",
          "# Create plots for each metric",
          "ggplot(alpha_long, aes(x = reorder(Sample, -Value), y = Value)) +",
          "  geom_bar(stat = \"identity\", fill = \"steelblue\") +",
          "  facet_wrap(~Metric, scales = \"free_y\") +",
          "  theme_minimal() +",
          "  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +",
          "  labs(x = \"Sample\", y = \"Diversity\")",
          "```",
          "",
          "## Diversity Statistics",
          "",
          "```{r alpha_stats}",
          "# Calculate summary statistics for each metric",
          "alpha_summary <- sapply(alpha_div, function(x) {",
          "  c(min = min(x), max = max(x), mean = mean(x), median = median(x))",
          "})",
          "",
          "# Print results",
          "knitr::kable(alpha_summary, digits = 2)",
          "```",
          ""
        )
      }
      
      if ("Beta Diversity" %in% input$reportSections) {
        report_content <- c(report_content,
          "# Beta Diversity",
          "",
          "Beta diversity measures the difference in microbial composition between samples.",
          "",
          "## Ordination",
          "",
          "```{r beta_div}",
          "# Transform to relative abundance",
          "ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))",
          "",
          "# Calculate Bray-Curtis distance",
          "dist_bc <- phyloseq::distance(ps_rel, method = \"bray\")",
          "",
          "# Perform PCoA ordination",
          "ord_pcoa <- ordinate(ps_rel, method = \"PCoA\", distance = dist_bc)",
          "",
          "# Get ordination coordinates",
          "ord_data <- data.frame(ord_pcoa$vectors[,1:2])",
          "colnames(ord_data) <- c(\"Axis1\", \"Axis2\")",
          "var_explained <- round(ord_pcoa$values$Relative_eig[1:2] * 100, 1)",
          "axis_labels <- paste0(\"Axis \", 1:2, \" (\", var_explained, \"%)\")",
          "",
          "# Add sample names",
          "ord_data$Sample <- rownames(ord_data)",
          "",
          "# Create PCoA plot",
          "ggplot(ord_data, aes(x = Axis1, y = Axis2, label = Sample)) +",
          "  geom_point(size = 3, alpha = 0.7, color = \"steelblue\") +",
          "  geom_text(nudge_y = 0.02, check_overlap = TRUE) +",
          "  theme_minimal() +",
          "  labs(x = axis_labels[1], y = axis_labels[2],",
          "       title = \"PCoA of Bray-Curtis Distances\")",
          "```",
          "",
          "## Sample Distance Heatmap",
          "",
          "```{r dist_heatmap}",
          "# Convert distance matrix to regular matrix",
          "dist_mat <- as.matrix(dist_bc)",
          "",
          "# Create heatmap",
          "pheatmap::pheatmap(dist_mat,",
          "                  main = \"Bray-Curtis Distance Heatmap\",",
          "                  display_numbers = FALSE,",
          "                  fontsize = 8,",
          "                  color = viridis::viridis(100))",
          "```",
          ""
        )
      }
      
      if ("Taxonomy" %in% input$reportSections) {
        report_content <- c(report_content,
          "# Taxonomy",
          "",
          "This section shows the taxonomic composition of samples at different levels.",
          "",
          "## Phylum Composition",
          "",
          "```{r tax_phylum}",
          "# Check if taxonomy table exists",
          "has_tax <- !is.null(tax_table(ps)) && \"Phylum\" %in% colnames(tax_table(ps))",
          "",
          "if (has_tax) {",
          "  # Agglomerate at phylum level",
          "  ps_phylum <- phyloseq::tax_glom(ps, taxrank = \"Phylum\")",
          "  ",
          "  # Transform to relative abundance",
          "  ps_phylum_rel <- phyloseq::transform_sample_counts(ps_phylum, function(x) x / sum(x))",
          "  ",
          "  # Create data frame for plotting",
          "  phylum_data <- phyloseq::psmelt(ps_phylum_rel)",
          "  ",
          "  # Handle NA values in Phylum",
          "  phylum_data$Phylum <- as.character(phylum_data$Phylum)",
          "  phylum_data$Phylum[is.na(phylum_data$Phylum)] <- \"Unknown\"",
          "  ",
          "  # Get top phyla (up to 10)",
          "  taxa_sums_val <- phyloseq::taxa_sums(ps_phylum_rel)",
          "  n_phyla <- min(10, length(unique(phylum_data$Phylum)))",
          "  top_phyla <- names(sort(taxa_sums_val, decreasing = TRUE)[1:n_phyla])",
          "  ",
          "  # Group remaining as Other",
          "  phylum_data$Phylum <- ifelse(phylum_data$OTU %in% top_phyla, ",
          "                              phylum_data$Phylum, \"Other\")",
          "  ",
          "  # Create plot",
          "  ggplot(phylum_data, aes(x = Sample, y = Abundance, fill = Phylum)) +",
          "    geom_bar(stat = \"identity\") +",
          "    theme_minimal() +",
          "    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +",
          "    labs(x = \"Sample\", y = \"Relative Abundance\") +",
          "    scale_fill_viridis_d()",
          "} else {",
          "  cat(\"No taxonomy data available at Phylum level.\")",
          "}",
          "```",
          "",
          "## Genus Composition",
          "",
          "```{r tax_genus}",
          "# Check if taxonomy table exists",
          "has_tax <- !is.null(tax_table(ps)) && \"Genus\" %in% colnames(tax_table(ps))",
          "",
          "if (has_tax) {",
          "  # Agglomerate at genus level",
          "  ps_genus <- phyloseq::tax_glom(ps, taxrank = \"Genus\")",
          "  ",
          "  # Transform to relative abundance",
          "  ps_genus_rel <- phyloseq::transform_sample_counts(ps_genus, function(x) x / sum(x))",
          "  ",
          "  # Create data frame for plotting",
          "  genus_data <- phyloseq::psmelt(ps_genus_rel)",
          "  ",
          "  # Handle NA values in Genus",
          "  genus_data$Genus <- as.character(genus_data$Genus)",
          "  genus_data$Genus[is.na(genus_data$Genus)] <- \"Unknown\"",
          "  ",
          "  # Get top genera (up to 20)",
          "  taxa_sums_val <- phyloseq::taxa_sums(ps_genus_rel)",
          "  n_genera <- min(20, length(unique(genus_data$Genus)))",
          "  top_genera <- names(sort(taxa_sums_val, decreasing = TRUE)[1:n_genera])",
          "  ",
          "  # Group remaining as Other",
          "  genus_data$Genus <- ifelse(genus_data$OTU %in% top_genera, ",
          "                            genus_data$Genus, \"Other\")",
          "  ",
          "  # Create plot",
          "  ggplot(genus_data, aes(x = Sample, y = Abundance, fill = Genus)) +",
          "    geom_bar(stat = \"identity\") +",
          "    theme_minimal() +",
          "    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +",
          "    labs(x = \"Sample\", y = \"Relative Abundance\") +",
          "    scale_fill_viridis_d(option = \"turbo\")",
          "} else {",
          "  cat(\"No taxonomy data available at Genus level.\")",
          "}",
          "```",
          "",
          "## Top Taxa",
          "",
          "```{r top_taxa}",
          "# Check if taxonomy table exists",
          "has_tax <- !is.null(tax_table(ps))",
          "",
          "if (has_tax) {",
          "  # Transform to relative abundance",
          "  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))",
          "  ",
          "  # Get mean abundance of each taxon",
          "  mean_abund <- taxa_sums(ps_rel) / nsamples(ps_rel)",
          "  ",
          "  # Get top 20 taxa",
          "  top20 <- names(sort(mean_abund, decreasing = TRUE)[1:20])",
          "  ",
          "  # Get taxonomy for these taxa",
          "  top_tax <- as.data.frame(tax_table(ps)[top20, ])",
          "  top_tax$Abundance <- mean_abund[rownames(top_tax)]",
          "  top_tax$Abundance <- round(top_tax$Abundance * 100, 2)",
          "  ",
          "  # Remove rownames for display",
          "  rownames(top_tax) <- NULL",
          "  ",
          "  # Display table",
          "  knitr::kable(top_tax, caption = \"Top 20 taxa by abundance, showing mean relative abundance (%)\")",
          "} else {",
          "  cat(\"No taxonomy data available.\")",
          "}",
          "```",
          ""
        )
      }
      
      if ("ASV Table" %in% input$reportSections) {
        report_content <- c(report_content,
          "# ASV Table",
          "",
          "This section provides summary information about the ASV table.",
          "",
          "## ASV Table Summary",
          "",
          "```{r asv_summary}",
          "# ASV table dimensions",
          "otu <- otu_table(ps)",
          "if (taxa_are_rows(otu)) {",
          "  n_taxa <- nrow(otu)",
          "  n_samples <- ncol(otu)",
          "} else {",
          "  n_taxa <- ncol(otu)",
          "  n_samples <- nrow(otu)",
          "}",
          "",
          "cat(\"**Number of ASVs:** \", n_taxa, \"\\n\\n\")",
          "cat(\"**Number of samples:** \", n_samples, \"\\n\\n\")",
          "cat(\"**Sparsity (% of zero values):** \", round(sum(otu == 0) / (n_taxa * n_samples) * 100, 1), \"%\\n\\n\")",
          "```",
          "",
          "## Sample Summary",
          "",
          "```{r sample_summary}",
          "# Create summary data frame for samples",
          "sample_df <- data.frame(",
          "  Sample = sample_names(ps),",
          "  ReadCount = sample_sums(ps),",
          "  ASVCount = colSums(otu > 0)",
          ")",
          "",
          "# Order by read count",
          "sample_df <- sample_df[order(sample_df$ReadCount, decreasing = TRUE), ]",
          "",
          "# Display sample summary",
          "knitr::kable(head(sample_df, 20), caption = \"Top 20 samples by read count\")",
          "```",
          "",
          "## ASV Summary",
          "",
          "```{r asv_top}",
          "# Create summary data frame for ASVs",
          "asv_df <- data.frame(",
          "  ASV = taxa_names(ps),",
          "  Abundance = taxa_sums(ps),",
          "  RelativeAbundance = taxa_sums(ps) / sum(taxa_sums(ps)) * 100,",
          "  Prevalence = rowSums(otu > 0) / n_samples * 100",
          ")",
          "",
          "# Add taxonomy if available",
          "if (!is.null(tax_table(ps))) {",
          "  tax <- as.data.frame(tax_table(ps))",
          "  for (col in colnames(tax)) {",
          "    asv_df[[col]] <- tax[asv_df$ASV, col]",
          "  }",
          "}",
          "",
          "# Order by abundance",
          "asv_df <- asv_df[order(asv_df$Abundance, decreasing = TRUE), ]",
          "",
          "# Format columns",
          "asv_df$RelativeAbundance <- round(asv_df$RelativeAbundance, 4)",
          "asv_df$Prevalence <- round(asv_df$Prevalence, 1)",
          "",
          "# Display ASV summary",
          "knitr::kable(head(asv_df, 20), caption = \"Top 20 ASVs by abundance\")",
          "```",
          ""
        )
      }
      
      # Add a conclusion section
      report_content <- c(report_content,
        "# Conclusion",
        "",
        "This report was automatically generated from the DADA2 Results Dashboard on `r format(Sys.time(), '%Y-%m-%d %H:%M:%S')`.",
        "",
        "For more detailed analysis, please use the interactive dashboard and R functions.",
        ""
      )
      
      # Write the report to the temporary file
      writeLines(report_content, tempReport)
      
      # Render the report
      rmarkdown::render(tempReport, output_file = file, quiet = TRUE)
    }
  )
}

# Display system information for performance tuning
print_system_info <- function() {
  cat("System Information for Performance Tuning:\n")
  cat("R Version:", R.version.string, "\n")
  cat("Available Cores:", parallel::detectCores(), "\n")
  cat("Memory Limit:", format(utils::memory.limit(), big.mark=","), "MB\n")
  
  # Check available memory
  if (requireNamespace("pryr", quietly = TRUE)) {
    cat("Free Memory:", format(pryr::mem_used(), big.mark=","), "bytes\n")
  }
  
  # Check if we're running in RStudio
  is_rstudio <- Sys.getenv("RSTUDIO") == "1"
  cat("Running in RStudio:", ifelse(is_rstudio, "Yes", "No"), "\n")
  
  # Check if parallel processing is available
  cat("Parallel Processing:", ifelse(parallel_cores > 1, "Enabled", "Disabled"), "\n")
  if (parallel_cores > 1) {
    cat("Using", parallel_cores, "cores for parallel operations\n")
  }
  
  # Check cache directory
  cat("Cache Directory:", cache_dir, "\n")
  cat("Cache Enabled:", ifelse(dir.exists(cache_dir), "Yes", "No"), "\n")
  
  # Print optional packages availability
  optional_pkgs <- c("data.table", "future", "promises", "memoise")
  for (pkg in optional_pkgs) {
    cat(pkg, "package:", ifelse(requireNamespace(pkg, quietly = TRUE), "Available", "Not available"), "\n")
  }
}

# Option to start in low-memory mode
low_memory_mode <- FALSE
if (exists("args") && "low-memory" %in% args) {
  low_memory_mode <- TRUE
  cat("Starting in low-memory mode\n")
}

# Create a simplified UI version to run first while loading data
ui_simple <- dashboardPage(
  dashboardHeader(title = "DADA2 Results (Loading...)"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "overview",
              fluidRow(
                box(width = 12, title = "DADA2 Dashboard - Loading Data", status = "primary",
                    div(
                      id = "loading-content",
                      tags$h3("Initializing dashboard..."),
                      tags$p("Please wait while we prepare the dashboard."),
                      tags$p("This may take a few moments for large datasets."),
                      tags$hr(),
                      tags$div(class = "progress progress-striped active",
                               tags$div(class = "progress-bar progress-bar-info", 
                                        style = "width: 100%",
                                        tags$span("Loading data...")
                               )
                      )
                    ))
              ),
              fluidRow(
                box(width = 12, title = "System Information", status = "info",
                    verbatimTextOutput("systemInfo"))
              )
      )
    )
  )
)

# Print system information before starting app
print_system_info()

# Run the application with memory optimization
app <- shinyApp(ui = ui, server = server)

# Run with custom options for better performance
options(shiny.maxRequestSize = 100 * 1024^2)  # Allow up to 100MB file uploads
options(future.globals.maxSize = 500 * 1024^2)  # Allow large data in future
options(shiny.reactlog = FALSE)  # Disable reactlog for better performance

# Run app with appropriate settings
if (low_memory_mode) {
  # Low memory mode settings
  runApp(app, port = 4321, launch.browser = TRUE, 
         quiet = TRUE, display.mode = "normal")
} else {
  # Standard mode with better performance
  runApp(app, port = 4321, launch.browser = TRUE, 
         quiet = TRUE, display.mode = "normal")
}