#!/usr/bin/env Rscript

# DADA2 Workflow Launcher App
# A Shiny application to configure and launch the optimized DADA2 workflow

# Function to set CRAN mirror automatically with smart detection
auto_set_cran_mirror <- function() {
  # Check if a CRAN mirror is already set
  repos <- getOption("repos")
  
  # If CRAN is already set to a working URL, use it
  if ("CRAN" %in% names(repos) && repos["CRAN"] != "@CRAN@" && repos["CRAN"] != "") {
    cat("Using existing CRAN mirror:", repos["CRAN"], "\n")
    return(repos["CRAN"])
  }
  
  # Primary mirrors to try (in order of preference)
  primary_mirrors <- c(
    "https://cloud.r-project.org",           # Global CDN
    "https://cran.rstudio.com"               # RStudio mirror
  )
  
  # Regional mirrors based on timezones
  timezone <- Sys.timezone()
  regional_mirrors <- list(
    US = c(
      "https://mirrors.nics.utk.edu/cran",    # USA (Tennessee)
      "https://cran.case.edu"                 # USA (Cleveland)
    ),
    Europe = c(
      "https://cran.wu.ac.at",                # Austria
      "https://cran.stat.unipd.it"            # Italy
    ),
    Asia = c(
      "https://mirror.las.iastate.edu/CRAN",  # China
      "https://cran.asia"                     # Singapore
    ),
    Australia = c(
      "https://cran.ms.unimelb.edu.au"        # Australia
    )
  )
  
  # Add regional mirrors based on timezone
  if (grepl("America|US", timezone)) {
    mirrors_to_try <- c(primary_mirrors, regional_mirrors$US)
  } else if (grepl("Europe|GB|DE|FR|IT|ES", timezone)) {
    mirrors_to_try <- c(primary_mirrors, regional_mirrors$Europe)
  } else if (grepl("Asia|Japan|China|Hong|Singapore|Korea", timezone)) {
    mirrors_to_try <- c(primary_mirrors, regional_mirrors$Asia) 
  } else if (grepl("Australia", timezone)) {
    mirrors_to_try <- c(primary_mirrors, regional_mirrors$Australia)
  } else {
    # Default to primary mirrors if region can't be determined
    mirrors_to_try <- primary_mirrors
  }
  
  # Try each mirror until one works
  for (mirror in mirrors_to_try) {
    tryCatch({
      # Try to download a small file from the mirror to test connectivity
      test_url <- paste0(mirror, "/PACKAGES")
      
      # Set a reasonable timeout
      timeout_option <- options(timeout = 5)
      on.exit(options(timeout_option), add = TRUE)
      
      # Try to connect to the mirror
      con <- url(test_url, "rb")
      close(con)
      
      # If we get here, the mirror is working
      cat("Setting CRAN mirror to", mirror, "\n")
      options(repos = c(CRAN = mirror))
      return(mirror)
    }, error = function(e) {
      # This mirror failed, try the next one
      cat("Mirror", mirror, "not available:", conditionMessage(e), "\n")
    })
  }
  
  # If all else fails, use the cloud.r-project.org mirror
  cat("No working mirrors found. Defaulting to cloud.r-project.org...\n")
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  return("https://cloud.r-project.org")
}

# Call the function to automatically set a CRAN mirror
auto_set_cran_mirror()

# Check for required packages and install if missing
required_packages <- c("shiny", "shinydashboard", "shinyFiles", "shinyjs", 
                      "rmarkdown", "DT", "future")

# Define special handling for Bioconductor packages
bioc_packages <- c("dada2")

# Install and load required packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    tryCatch({
      install.packages(pkg)
    }, error = function(e) {
      cat("Error installing package", pkg, ":", conditionMessage(e), "\n")
      cat("Trying to find alternative mirror...\n")
      # Re-attempt mirror selection and try again
      mirror <- auto_set_cran_mirror()
      install.packages(pkg, repos = mirror)
    })
  }
  library(pkg, character.only = TRUE)
}

# Check if we need BiocManager
if (any(!sapply(bioc_packages, requireNamespace, quietly = TRUE))) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    tryCatch({
      install.packages("BiocManager")
    }, error = function(e) {
      cat("Error installing BiocManager:", conditionMessage(e), "\n")
      cat("Trying to find alternative mirror...\n")
      # Re-attempt mirror selection and try again
      mirror <- auto_set_cran_mirror()
      install.packages("BiocManager", repos = mirror)
    })
  }
  
  # Only check for dada2 - don't try to install automatically as it's large
  if (!requireNamespace("dada2", quietly = TRUE)) {
    cat("NOTE: The dada2 package is not installed. It will be required for workflow execution.\n")
    cat("It can be installed using: BiocManager::install('dada2')\n")
  } else {
    library(dada2)
  }
}

# UI Definition
ui <- dashboardPage(
  skin = "blue",
  
  # Header
  dashboardHeader(
    title = "DADA2 Workflow Launcher",
    titleWidth = 300
  ),
  
  # Sidebar
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Workflow Setup", tabName = "setup", icon = icon("sliders")),
      menuItem("Advanced Settings", tabName = "advanced", icon = icon("cogs")),
      menuItem("Run Workflow", tabName = "run", icon = icon("play-circle")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  # Body
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f8f9fa;
        }
        .box {
          box-shadow: 0 1px 3px rgba(0,0,0,0.12), 0 1px 2px rgba(0,0,0,0.24);
        }
        .btn-primary {
          background-color: #3c8dbc;
          border-color: #367fa9;
        }
        .btn-primary:hover {
          background-color: #367fa9;
          border-color: #204d74;
        }
        .shiny-notification {
          position: fixed;
          top: 70px;
          right: 10px;
        }
      "))
    ),
    
    tabItems(
      # Setup Tab
      tabItem(tabName = "setup",
        fluidRow(
          box(
            title = "Input Data Configuration", width = 12, status = "primary",
            solidHeader = TRUE,
            fluidRow(
              column(6,
                radioButtons("runMode", "Processing Mode:",
                  choices = list(
                    "Single Run" = "single",
                    "Multi-Run (Multiple Sequencing Runs)" = "multi"
                  ),
                  selected = "single"
                ),
                conditionalPanel(
                  condition = "input.runMode == 'single'",
                  shinyDirButton("dataDir", "Select Data Directory", "Choose directory containing FASTQ files")
                ),
                conditionalPanel(
                  condition = "input.runMode == 'multi'",
                  shinyDirButton("runDir", "Select Runs Directory", "Choose directory containing run subdirectories")
                )
              ),
              column(6,
                textOutput("selectedDirText"),
                tags$br(),
                checkboxInput("bigData", "Enable Big Data Mode (optimized for memory usage)", FALSE),
                numericInput("numCores", "Number of CPU Cores to Use", 
                             value = max(1, parallel::detectCores() - 1),
                             min = 1, max = parallel::detectCores())
              )
            )
          ),
          
          box(
            title = "Output Configuration", width = 12, status = "primary",
            solidHeader = TRUE,
            fluidRow(
              column(6,
                shinyDirButton("outputDir", "Select Output Directory", "Choose where to save results"),
                textOutput("selectedOutDirText"),
                tags$br(),
                textInput("outputFileBase", "Output File Base Name:", "dada2_workflow_report")
              ),
              column(6,
                checkboxInput("generateReport", "Generate Report", TRUE),
                conditionalPanel(
                  condition = "input.generateReport",
                  selectInput("reportFormat", "Report Format:",
                    choices = list(
                      "HTML" = "html_document",
                      "PDF" = "pdf_document"
                    ),
                    selected = "html_document"
                  )
                ),
                checkboxInput("launchDashboard", "Launch Dashboard After Completion", TRUE)
              )
            )
          )
        )
      ),
      
      # Advanced Settings Tab
      tabItem(tabName = "advanced",
        fluidRow(
          box(
            title = "Taxonomy Settings", width = 12, status = "primary",
            solidHeader = TRUE,
            fluidRow(
              column(6,
                selectInput("taxonomyMethods", "Taxonomy Assignment Methods:",
                  choices = list(
                    "SILVA DADA2 (RDP classifier)" = "SILVA_DADA2",
                    "IDTAXA with SILVA" = "IDTAXA_SILVA",
                    "RDP Database" = "RDP_DADA2",
                    "All Available Methods" = "ALL"
                  ),
                  selected = "SILVA_DADA2",
                  multiple = TRUE
                ),
                numericInput("taxonomyConfThresh", "Taxonomy Confidence Threshold (%)", 
                             value = 80, min = 50, max = 100)
              ),
              column(6,
                checkboxInput("phylogeneticTree", "Build Phylogenetic Tree", TRUE),
                conditionalPanel(
                  condition = "input.phylogeneticTree",
                  selectInput("treeMethod", "Tree Construction Method:",
                    choices = list(
                      "Neighbor-Joining (Fast)" = "NJ",
                      "Maximum Likelihood (Accurate)" = "ML"
                    ),
                    selected = "NJ"
                  ),
                  sliderInput("bootstrapReps", "Bootstrap Replicates:", 
                              min = 0, max = 100, value = 0, step = 10)
                )
              )
            )
          ),
          
          box(
            title = "Filtering & Processing Settings", width = 12, status = "primary",
            solidHeader = TRUE,
            fluidRow(
              column(6,
                checkboxInput("autoOptimizeTrunc", "Auto-Optimize Truncation Lengths", TRUE),
                conditionalPanel(
                  condition = "!input.autoOptimizeTrunc",
                  numericInput("truncLenF", "Forward Read Truncation Length:", 0),
                  numericInput("truncLenR", "Reverse Read Truncation Length:", 0)
                ),
                checkboxInput("autoOptimizeMaxEE", "Auto-Optimize Expected Error Thresholds", TRUE),
                conditionalPanel(
                  condition = "!input.autoOptimizeMaxEE",
                  numericInput("maxEEF", "Forward MaxEE:", 2, min = 0, max = 10, step = 0.5),
                  numericInput("maxEER", "Reverse MaxEE:", 2, min = 0, max = 10, step = 0.5)
                )
              ),
              column(6,
                checkboxInput("performRarefaction", "Perform Rarefaction Analysis", TRUE),
                checkboxInput("enableCheckpointing", "Enable Checkpointing System", TRUE),
                selectInput("poolingStrategy", "Reads Pooling Strategy:",
                  choices = list(
                    "Pseudo-pooling (Recommended)" = "pseudo",
                    "No Pooling" = "FALSE",
                    "True Pooling" = "TRUE"
                  ),
                  selected = "pseudo"
                ),
                conditionalPanel(
                  condition = "input.runMode == 'multi'",
                  checkboxInput("batchEffectAnalysis", "Perform Batch Effect Analysis", TRUE)
                )
              )
            )
          )
        )
      ),
      
      # Run Tab
      tabItem(tabName = "run",
        fluidRow(
          box(
            title = "Workflow Summary", width = 12, status = "primary",
            solidHeader = TRUE,
            verbatimTextOutput("workflowSummary"),
            tags$br(),
            actionButton("runWorkflow", "Run DADA2 Workflow", 
                        icon = icon("play"), 
                        class = "btn-lg btn-primary"),
            tags$br(), tags$br(),
            div(id = "progressBox", style = "display: none",
              progressBar(id = "progressBar", value = 0, 
                        title = "Running workflow...",
                        display_pct = TRUE),
              verbatimTextOutput("workflowLog")
            )
          )
        )
      ),
      
      # Help Tab
      tabItem(tabName = "help",
        fluidRow(
          box(
            title = "Help & Documentation", width = 12, status = "primary",
            solidHeader = TRUE,
            tags$h4("About the DADA2 Optimized Workflow"),
            tags$p("This application provides a user-friendly interface to run the optimized DADA2 workflow for processing 16S rRNA gene amplicon data. The workflow identifies exact amplicon sequence variants (ASVs) with higher resolution than traditional OTU-based methods."),
            tags$h4("Key Features of the Optimized Workflow"),
            tags$ul(
              tags$li(tags$b("Automatic parameter optimization"), "- Determines optimal filtering parameters based on your data"),
              tags$li(tags$b("Multi-run support"), "- Process and integrate data from multiple sequencing runs"),
              tags$li(tags$b("Memory optimization"), "- Efficiently handles large datasets with batched processing"),
              tags$li(tags$b("Checkpoint system"), "- Resume pipeline from intermediate steps if interrupted"),
              tags$li(tags$b("Taxonomy confidence"), "- Bootstrap confidence values for taxonomic assignments"),
              tags$li(tags$b("Comprehensive reporting"), "- Generate detailed reports and visualizations")
            ),
            tags$h4("Input Requirements"),
            tags$p("Your input data should be paired-end Illumina sequencing reads in FASTQ format, preferably gzipped. Files should follow naming conventions with _R1 and _R2 designations for forward and reverse reads."),
            tags$h4("Directory Structure"),
            tags$p("For single-run mode, all FASTQ files should be in one directory. For multi-run mode, each run should be in a separate subdirectory."),
            tags$h4("Need More Help?"),
            tags$p("For more detailed documentation, please refer to:", 
                  tags$a("DADA2 Tutorial", href = "https://benjjneb.github.io/dada2/tutorial.html", target = "_blank"), 
                  "and", 
                  tags$a("DADA2 Big Data Processing", href = "https://benjjneb.github.io/dada2/bigdata.html", target = "_blank"))
          )
        )
      )
    )
  )
)

# Server Definition
server <- function(input, output, session) {
  # Directory selection handlers
  volumes <- c(Home = "~", "R Installation" = R.home(), getVolumes()())
  
  shinyDirChoose(input, "dataDir", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "runDir", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "outputDir", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  # Reactive values to store selected directories
  selectedDataDir <- reactiveVal("")
  selectedRunDir <- reactiveVal("")
  selectedOutputDir <- reactiveVal("")
  
  # Update selected data directory text
  observe({
    if (input$runMode == "single" && !is.null(input$dataDir)) {
      path <- parseDirPath(volumes, input$dataDir)
      if (length(path) > 0) {
        selectedDataDir(path)
        output$selectedDirText <- renderText({
          paste("Selected data directory:", path)
        })
      }
    } else if (input$runMode == "multi" && !is.null(input$runDir)) {
      path <- parseDirPath(volumes, input$runDir)
      if (length(path) > 0) {
        selectedRunDir(path)
        output$selectedDirText <- renderText({
          paste("Selected runs directory:", path)
        })
      }
    }
  })
  
  # Update selected output directory text
  observe({
    if (!is.null(input$outputDir)) {
      path <- parseDirPath(volumes, input$outputDir)
      if (length(path) > 0) {
        selectedOutputDir(path)
        output$selectedOutDirText <- renderText({
          paste("Selected output directory:", path)
        })
      }
    }
  })
  
  # Generate workflow summary
  output$workflowSummary <- renderText({
    req(input$runMode)
    
    # Required directories based on run mode
    if (input$runMode == "single") {
      req(selectedDataDir())
    } else {
      req(selectedRunDir())
    }
    
    summary_text <- paste0(
      "DADA2 Workflow Configuration Summary:\n\n",
      "Processing Mode: ", ifelse(input$runMode == "single", "Single Run", "Multi-Run"), "\n",
      ifelse(input$runMode == "single", 
            paste0("Data Directory: ", selectedDataDir(), "\n"),
            paste0("Runs Directory: ", selectedRunDir(), "\n")),
      "Output Directory: ", ifelse(selectedOutputDir() != "", selectedOutputDir(), "Default"), "\n",
      "CPU Cores: ", input$numCores, "\n",
      "Big Data Mode: ", ifelse(input$bigData, "Enabled", "Disabled"), "\n",
      "Generate Report: ", ifelse(input$generateReport, paste0("Yes (", input$reportFormat, ")"), "No"), "\n",
      "Launch Dashboard: ", ifelse(input$launchDashboard, "Yes", "No"), "\n\n",
      
      "Advanced Settings:\n",
      "Taxonomy Methods: ", paste(input$taxonomyMethods, collapse = ", "), "\n",
      "Build Phylogenetic Tree: ", ifelse(input$phylogeneticTree, 
                                        paste0("Yes (", input$treeMethod, ", ", 
                                              ifelse(input$bootstrapReps > 0, 
                                                    paste0(input$bootstrapReps, " bootstrap replicates)"), 
                                                    "no bootstrapping)")), 
                                        "No"), "\n",
      "Auto-Optimize Parameters: ", 
      ifelse(input$autoOptimizeTrunc && input$autoOptimizeMaxEE, "All Parameters", 
            ifelse(!input$autoOptimizeTrunc && !input$autoOptimizeMaxEE, "None", 
                  ifelse(input$autoOptimizeTrunc, "Truncation Only", "MaxEE Only"))), "\n",
      "Pooling Strategy: ", input$poolingStrategy, "\n",
      "Rarefaction Analysis: ", ifelse(input$performRarefaction, "Yes", "No"), "\n",
      "Checkpointing: ", ifelse(input$enableCheckpointing, "Enabled", "Disabled"), "\n",
      ifelse(input$runMode == "multi" && input$batchEffectAnalysis, 
            "Batch Effect Analysis: Enabled\n", "")
    )
    
    return(summary_text)
  })
  
  # Run workflow button handler
  observeEvent(input$runWorkflow, {
    # Required directories based on run mode
    if (input$runMode == "single") {
      if (selectedDataDir() == "") {
        showNotification("Please select a data directory first.", type = "error")
        return()
      }
      data_path <- selectedDataDir()
    } else {
      if (selectedRunDir() == "") {
        showNotification("Please select a runs directory first.", type = "error")
        return()
      }
      data_path <- selectedRunDir()
    }
    
    # Show progress box
    shinyjs::show("progressBox")
    
    # Construct the rendering parameters
    params <- list(
      multi_run = input$runMode == "multi",
      run_dir = ifelse(input$runMode == "multi", data_path, NULL),
      big_data = input$bigData,
      generate_report = input$generateReport,
      output_format = input$reportFormat,
      output_dir = ifelse(selectedOutputDir() != "", selectedOutputDir(), "reports"),
      output_file = input$outputFileBase
    )
    
    # Add parameter overrides based on advanced settings
    if (!input$autoOptimizeTrunc) {
      params$truncLen_forward <- input$truncLenF
      params$truncLen_reverse <- input$truncLenR
    }
    
    if (!input$autoOptimizeMaxEE) {
      params$maxEE_forward <- input$maxEEF
      params$maxEE_reverse <- input$maxEER
    }
    
    # Set taxonomy methods
    if ("ALL" %in% input$taxonomyMethods) {
      params$taxonomy_methods <- list(
        SILVA_DADA2 = TRUE,
        IDTAXA_SILVA = TRUE,
        RDP_DADA2 = TRUE
      )
    } else {
      params$taxonomy_methods <- list(
        SILVA_DADA2 = "SILVA_DADA2" %in% input$taxonomyMethods,
        IDTAXA_SILVA = "IDTAXA_SILVA" %in% input$taxonomyMethods,
        RDP_DADA2 = "RDP_DADA2" %in% input$taxonomyMethods
      )
    }
    
    # Set phylogenetic tree parameters
    if (input$phylogeneticTree) {
      params$build_tree <- TRUE
      params$tree_method <- input$treeMethod
      params$bootstrap_reps <- input$bootstrapReps
    } else {
      params$build_tree <- FALSE
    }
    
    # Set other parameters
    params$pooling_option <- input$poolingStrategy
    params$perform_rarefaction <- input$performRarefaction
    params$enable_checkpointing <- input$enableCheckpointing
    params$batch_effect_analysis <- input$runMode == "multi" && input$batchEffectAnalysis
    params$taxonomy_conf_threshold <- input$taxonomyConfThresh
    
    # Create a temporary file for logging
    log_file <- tempfile(pattern = "dada2_log_", fileext = ".txt")
    
    # Set up a process to run the workflow
    # We need to use callr or system2 to run the workflow in a separate R process
    # and capture the output
    
    # First create a temporary R script that will run the workflow
    script_file <- tempfile(pattern = "dada2_script_", fileext = ".R")
    
    script_content <- paste0(
      "library(rmarkdown)\n\n",
      "# Set number of cores\n",
      "future::plan(future::multisession, workers = ", input$numCores, ")\n\n",
      "# Set parameters\n",
      "params <- list(\n",
      "  multi_run = ", ifelse(params$multi_run, "TRUE", "FALSE"), ",\n",
      "  run_dir = ", ifelse(is.null(params$run_dir), "NULL", paste0("'", params$run_dir, "'")), ",\n",
      "  big_data = ", ifelse(params$big_data, "TRUE", "FALSE"), ",\n",
      "  generate_report = ", ifelse(params$generate_report, "TRUE", "FALSE"), ",\n",
      "  output_format = '", params$output_format, "',\n",
      "  output_dir = '", params$output_dir, "',\n",
      "  output_file = '", params$output_file, "'\n",
      ")\n\n",
      "# Run the workflow\n",
      "tryCatch({\n",
      "  cat('Starting DADA2 workflow execution...\\n')\n",
      "  render('dada2_workflow_optimize.Rmd', params = params, envir = new.env())\n",
      "  cat('DADA2 workflow completed successfully.\\n')\n",
      "  if (", ifelse(input$launchDashboard, "TRUE", "FALSE"), ") {\n",
      "    cat('Launching dashboard...\\n')\n",
      "    source('run_dashboard.R')\n",
      "  }\n",
      "}, error = function(e) {\n",
      "  cat('Error in DADA2 workflow execution:\\n')\n",
      "  cat(conditionMessage(e), '\\n')\n",
      "})\n"
    )
    
    # Write the script to file
    writeLines(script_content, script_file)
    
    # Run the script and capture output
    system2("Rscript", 
           args = c(script_file),
           stdout = log_file, 
           stderr = log_file,
           wait = FALSE)
    
    # Start a timer to update the progress bar and log display
    progress_value <- 0
    
    # Function to update the progress
    update_progress <- function() {
      # Read the log file
      if (file.exists(log_file)) {
        log_content <- readLines(log_file)
        
        # Update the log display
        output$workflowLog <- renderText({
          paste(tail(log_content, 20), collapse = "\n")
        })
        
        # Update progress based on keywords in the log
        if (any(grepl("Starting DADA2 workflow", log_content))) {
          progress_value <<- 5
        }
        if (any(grepl("quality filtering", log_content))) {
          progress_value <<- 20
        }
        if (any(grepl("denoising", log_content))) {
          progress_value <<- 40
        }
        if (any(grepl("merging", log_content))) {
          progress_value <<- 60
        }
        if (any(grepl("chimera removal", log_content))) {
          progress_value <<- 70
        }
        if (any(grepl("taxonomy assignment", log_content))) {
          progress_value <<- 80
        }
        if (any(grepl("phyloseq object", log_content))) {
          progress_value <<- 90
        }
        if (any(grepl("workflow completed successfully", log_content))) {
          progress_value <<- 100
          if (any(grepl("Launching dashboard", log_content))) {
            showNotification("Dashboard is launching in a separate window.", 
                           type = "message", duration = 10)
          }
        }
        
        # Check for errors
        if (any(grepl("Error in DADA2 workflow", log_content))) {
          showNotification("Error in DADA2 workflow. Check the log for details.", 
                         type = "error", duration = NULL)
        }
        
        # Update progress bar
        updateProgressBar(session = session, id = "progressBar", value = progress_value)
        
        # Continue updating if not done
        if (progress_value < 100 && !any(grepl("Error in DADA2 workflow", log_content))) {
          invalidateLater(1000, session)
        } else {
          # Execution finished
          if (progress_value == 100) {
            showNotification("DADA2 workflow completed successfully!", 
                           type = "message", duration = NULL)
          }
        }
      } else {
        # Log file doesn't exist yet, keep checking
        invalidateLater(1000, session)
      }
    }
    
    # Start progress updates
    update_progress()
  })
}

# Run the application
shinyApp(ui = ui, server = server)