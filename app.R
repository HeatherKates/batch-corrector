# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(readxl)
library(sva)
library(dplyr)
library(tibble)
library(shinyjs)
options(shiny.maxRequestSize = 30*1024^3)
# UI
ui <- dashboardPage(
  dashboardHeader(title = "Batch Corrector"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
      menuItem("Configure Data", tabName = "configure", icon = icon("cogs")),
      menuItem("Batch Correction", tabName = "correction", icon = icon("magic")),
      menuItem("Download Results", tabName = "download", icon = icon("download"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    tabItems(
      # Upload tab
      tabItem(tabName = "upload",
              fluidRow(
                # Add demo data option at the top
                box(title = "Try the Demo", status = "success", solidHeader = TRUE, width = 12,
                    p("Want to try the app with example data first?"),
                    actionButton("load_demo", "Load Test Data", 
                                 class = "btn-success btn-lg",
                                 icon = icon("play")),
                    br(), br(),
                    div(id = "demo_info", style = "display: none;",
                        div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; border: 1px solid #c3e6cb;",
                            p(strong("Demo data loaded!"), 
                              "This includes simulated RNA-seq count data with 2000 features, 48 samples across 3 batches and 2 biological conditions. 
                Go to the Configure Data tab to proceed."))
                    ),
                    hr()
                )
              ),
              fluidRow(
                box(title = "Upload Feature Data", status = "primary", solidHeader = TRUE, width = 6,
                    fileInput("count_file", "Choose Feature (counts,intensities, etc.) Data File (.csv or .xlsx)",
                              accept = c(".csv", ".xlsx")),
                    conditionalPanel(
                      condition = "output.count_uploaded",
                      h4("Count Data Preview:"),
                      DT::dataTableOutput("count_preview")
                    )
                ),
                
                box(title = "Upload Sample Metadata", status = "primary", solidHeader = TRUE, width = 6,
                    fileInput("sample_file", "Choose Sample Metadata File (.csv or .xlsx)",
                              accept = c(".csv", ".xlsx")),
                    conditionalPanel(
                      condition = "output.sample_uploaded",
                      h4("Sample Data Preview:"),
                      DT::dataTableOutput("sample_preview")
                    )
                )
              )
      ),
      
      # Configure tab
      tabItem(tabName = "configure",
              fluidRow(
                box(title = "Configure Count Data", status = "warning", solidHeader = TRUE, width = 6,
                    conditionalPanel(
                      condition = "output.count_uploaded",
                      selectInput("feature_col", "Select Feature Name Column:", choices = NULL),
                      selectInput("drop_cols", "Select Columns to Drop:", choices = NULL, multiple = TRUE),
                      br(),
                      verbatimTextOutput("count_config_summary")
                    )
                ),
                
                box(title = "Configure Sample Data", status = "warning", solidHeader = TRUE, width = 6,
                    conditionalPanel(
                      condition = "output.sample_uploaded",
                      selectInput("batch_col", "Select Batch Variable Column:", choices = NULL),
                      selectInput("bio_col", "Select Biological Variable Column:", choices = NULL),
                      selectInput("sample_name_col", "Select Sample Name Column:", choices = NULL),
                      br(),
                      verbatimTextOutput("sample_config_summary")
                    )
                )
              ),
              
              fluidRow(
                box(title = "Preview Sample Matching", status = "info", solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.ready_to_preview",
                      actionButton("preview_matching", "Preview Sample Matching", 
                                   class = "btn-info btn-lg"),
                      br(), br()
                    ),
                    conditionalPanel(
                      condition = "output.matching_previewed",
                      verbatimTextOutput("matching_summary"),
                      br(),
                      h4("Data before batch correction:"),
                      plotlyOutput("pca_before", height = "400px")
                    )
                )
              )
      ),
      
      # Correction tab
      tabItem(tabName = "correction",
              fluidRow(
                # Add this box after the sample configuration boxes in the Configure tab
                box(title = "Select Batch Correction Method", status = "info", solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.ready_to_preview",
                      radioButtons("batch_method", "Choose batch correction method:",
                                   choices = list(
                                     "ComBat-Seq (recommended for RNA-seq other count-based data)" = "combat_seq",
                                     "ComBat (recommended for metabolomics, intensity data)" = "combat_original"
                                   ),
                                   selected = "combat_seq"),
                      br(),
                      div(style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #007bff;",
                          h5("Method Guide:"),
                          tags$ul(
                            tags$li(strong("ComBat-Seq:"), "Uses negative binomial model, works directly with count data. Best for RNA-seq or other data with count-like properties and overdispersion."),
                            tags$li(strong("ComBat:"), "Uses log transformation + linear model. Best for metabolomics peak intensities, microarray data, or other continuous measurements with log-normal distributions.")
                          ),
                          p(strong("Note:"), "PCA plots will always use log-transformed data for better visualization, regardless of which batch correction method you choose.")
                      )
                    )
                ),
                box(title = "Batch Correction", status = "success", solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.ready_for_correction",
                      actionButton("run_combat", "Run ComBat Batch Correction", 
                                   class = "btn-success btn-lg"),
                      br(), br()
                    ),
                    conditionalPanel(
                      condition = "output.correction_complete",
                      h4("Data after batch correction:"),
                      plotlyOutput("pca_after", height = "400px")
                    )
                )
              )
      ),
      
      # Download tab
      tabItem(tabName = "download",
              fluidRow(
                box(title = "Download Results", status = "info", solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.correction_complete",
                      h4("Download batch-corrected data:"),
                      downloadButton("download_corrected", "Download Batch-Corrected Matrix", 
                                     class = "btn-primary"),
                      br(), br(),
                      downloadButton("download_metadata", "Download Updated Sample Metadata", 
                                     class = "btn-primary"),
                      br(), br(),
                      h4("Usage Notes:"),
                      tags$ul(
                        tags$li("The batch-corrected matrix should be used for visualization and clustering, but NOT for differential expression analysis"),
                        tags$li("For differential expression (DESeq2, limma-voom, edgeR), use the original raw counts with the batch factor from the metadata"),
                        tags$li("Include the batch correction factor as a covariate in your model design")
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
  
  # Reactive values to store data
  values <- reactiveValues(
    count_data = NULL,
    sample_data = NULL,
    matched_data = NULL,
    corrected_data = NULL,
    pca_before = NULL,
    pca_after = NULL,
    matching_previewed = FALSE
  )
  
  # File upload functions
  load_file <- function(file_path) {
    ext <- tools::file_ext(file_path)
    if(ext == "csv") {
      return(read.csv(file_path, stringsAsFactors = FALSE))
    } else if(ext %in% c("xlsx", "xls")) {
      return(read_excel(file_path))
    }
    return(NULL)
  }
  
  # Upload count data
  observeEvent(input$count_file, {
    req(input$count_file)
    values$count_data <- load_file(input$count_file$datapath)
    updateSelectInput(session, "feature_col", choices = colnames(values$count_data))
    updateSelectInput(session, "drop_cols", choices = colnames(values$count_data))
    values$matching_previewed <- FALSE  # Reset preview status
    updateRadioButtons(session, "batch_method", selected = "combat_seq")  # Reset to default
    
  })
  
  # Upload sample data
  observeEvent(input$sample_file, {
    req(input$sample_file)
    values$sample_data <- load_file(input$sample_file$datapath)
    updateSelectInput(session, "batch_col", choices = colnames(values$sample_data))
    updateSelectInput(session, "bio_col", choices = colnames(values$sample_data))
    updateSelectInput(session, "sample_name_col", choices = colnames(values$sample_data))
    values$matching_previewed <- FALSE  # Reset preview status
    updateRadioButtons(session, "batch_method", selected = "combat_seq")  # Reset to default
    
  })
  # Load demo data
  observeEvent(input$load_demo, {
    tryCatch({
      # Load test data files
      if(file.exists("test-data/counts-data.csv") && file.exists("test-data/sample-data.csv")) {
        values$count_data <- read.csv("test-data/counts-data.csv", stringsAsFactors = FALSE)
        values$sample_data <- read.csv("test-data/sample-data.csv", stringsAsFactors = FALSE)
        
        # Update dropdown choices
        updateSelectInput(session, "feature_col", choices = colnames(values$count_data), selected = "gene_id")
        updateSelectInput(session, "drop_cols", choices = colnames(values$count_data), 
                          selected = c("description", "notes", "quality_flag"))
        updateSelectInput(session, "batch_col", choices = colnames(values$sample_data), selected = "batch")
        updateSelectInput(session, "bio_col", choices = colnames(values$sample_data), selected = "biological_var")
        updateSelectInput(session, "sample_name_col", choices = colnames(values$sample_data), selected = "sample_name")
        
        # Reset states
        values$matching_previewed <- FALSE
        updateRadioButtons(session, "batch_method", selected = "combat_seq")
        
        # Show demo info message
        shinyjs::show("demo_info")
        
        showNotification("Demo data loaded successfully! Go to Configure Data tab.", type = "message")
        
      } else {
        showNotification("Test data files not found. Please ensure test-data/counts-data.csv and test-data/sample-data.csv exist.", 
                         type = "error")
      }
      
    }, error = function(e) {
      showNotification(paste("Error loading demo data:", e$message), type = "error")
    })
  })
  
  # Reset preview when configuration changes
  observeEvent({
    list(input$feature_col, input$drop_cols, input$batch_col, input$bio_col, input$sample_name_col)
  }, {
    values$matching_previewed <- FALSE
  })
  
  # Output conditions
  output$count_uploaded <- reactive({ !is.null(values$count_data) })
  output$sample_uploaded <- reactive({ !is.null(values$sample_data) })
  output$ready_to_preview <- reactive({
    !is.null(values$count_data) && !is.null(values$sample_data) &&
      !is.null(input$feature_col) && !is.null(input$batch_col) && 
      !is.null(input$bio_col) && !is.null(input$sample_name_col) &&
      input$feature_col != "" && input$batch_col != "" && 
      input$bio_col != "" && input$sample_name_col != ""
  })
  output$matching_previewed <- reactive({ values$matching_previewed })
  output$ready_for_correction <- reactive({ 
    values$matching_previewed && !is.null(values$matched_data) 
  })
  output$correction_complete <- reactive({ 
    !is.null(values$corrected_data) 
  })
  
  outputOptions(output, "count_uploaded", suspendWhenHidden = FALSE)
  outputOptions(output, "sample_uploaded", suspendWhenHidden = FALSE)
  outputOptions(output, "ready_to_preview", suspendWhenHidden = FALSE)
  outputOptions(output, "matching_previewed", suspendWhenHidden = FALSE)
  outputOptions(output, "ready_for_correction", suspendWhenHidden = FALSE)
  outputOptions(output, "correction_complete", suspendWhenHidden = FALSE)
  
  # Preview tables
  output$count_preview <- DT::renderDataTable({
    req(values$count_data)
    DT::datatable(values$count_data[1:min(10, nrow(values$count_data)), 1:min(10, ncol(values$count_data))], 
                  options = list(scrollX = TRUE))
  })
  
  output$sample_preview <- DT::renderDataTable({
    req(values$sample_data)
    DT::datatable(values$sample_data, options = list(scrollX = TRUE))
  })
  
  # Preview sample matching
  observeEvent(input$preview_matching, {
    req(values$count_data, values$sample_data, 
        input$feature_col, input$batch_col, input$bio_col, input$sample_name_col)
    
    tryCatch({
      # Process count data
      count_processed <- values$count_data
      
      # Remove manually dropped columns
      if(!is.null(input$drop_cols) && length(input$drop_cols) > 0) {
        count_processed <- count_processed[, !colnames(count_processed) %in% input$drop_cols, drop = FALSE]
      }
      
      feature_col <- input$feature_col
      potential_sample_cols <- setdiff(colnames(count_processed), feature_col)
      
      # Process sample data
      sample_processed <- values$sample_data[, c(input$sample_name_col, input$batch_col, input$bio_col)]
      colnames(sample_processed) <- c("sample_name", "batch", "biological_var")
      # Match samples using substring matching
      matched_samples <- c()
      sample_mapping <- c()
      
      for(sample_name in sample_processed$sample_name) {
        # Normalize by converting both - and . to the same character for matching
        normalized_sample <- gsub("[-.]", ".", sample_name)
        normalized_cols <- gsub("[-.]", ".", potential_sample_cols)
        
        # Find matches using normalized versions
        match_indices <- grep(normalized_sample, normalized_cols, fixed = TRUE)
        if(length(match_indices) > 0) {
          # Use the original column name (not the normalized version)
          matched_samples <- c(matched_samples, potential_sample_cols[match_indices[1]])
          sample_mapping <- c(sample_mapping, sample_name)
        }
      }
      
      if(length(matched_samples) > 0) {
        # Create count matrix with only matched samples
        count_matrix <- as.matrix(count_processed[, matched_samples, drop = FALSE])
        rownames(count_matrix) <- count_processed[[feature_col]]
        
        # Subset sample data to matched samples only
        sample_data_matched <- sample_processed[sample_processed$sample_name %in% sample_mapping, ]
        
        # Reorder sample data to match count matrix column order
        sample_order <- match(sample_mapping, sample_data_matched$sample_name)
        sample_data_matched <- sample_data_matched[sample_order, ]
        
        # Store matched data
        values$matched_data <- list(
          count_matrix = count_matrix,
          sample_data = sample_data_matched,
          matched_samples = matched_samples
        )
        
        # Perform PCA before correction
        values$pca_before <- perform_pca(count_matrix, sample_data_matched)
        
        # Set preview status
        values$matching_previewed <- TRUE
        
        showNotification(paste("Successfully matched", length(matched_samples), "samples!"), 
                         type = "message")
        
      } else {
        # No matches found
        values$matched_data <- NULL
        values$pca_before <- NULL
        values$matching_previewed <- FALSE
        showNotification("No sample matches found! Check your sample names.", type = "error")
      }
      
    }, error = function(e) {
      # If there's an error in processing
      values$matched_data <- NULL
      values$pca_before <- NULL
      values$matching_previewed <- FALSE
      showNotification(paste("Error processing data:", e$message), type = "error")
    })
  }) # End
  
  # PCA function - always log transform for visualization
  perform_pca <- function(count_matrix, sample_data) {
    # Always log transform for PCA visualization
    # Add small pseudocount to handle zeros
    count_matrix <- log2(count_matrix + 1)
    
    # Select top 1000 most variable features
    if(nrow(count_matrix) > 1000) {
      vars <- apply(count_matrix, 1, var, na.rm = TRUE)
      top_features <- names(sort(vars, decreasing = TRUE)[1:1000])
      count_matrix <- count_matrix[top_features, ]
    }
    
    # Remove features with zero variance
    vars <- apply(count_matrix, 1, var, na.rm = TRUE)
    count_matrix <- count_matrix[vars > 0 & !is.na(vars), ]
    
    # Transpose for PCA (samples as rows)
    pca_data <- t(count_matrix)
    pca_result <- prcomp(pca_data, scale. = TRUE, center = TRUE)
    
    # Create data frame for plotting
    pca_df <- data.frame(
      PC1 = pca_result$x[, 1],
      PC2 = pca_result$x[, 2],
      sample_name = colnames(count_matrix),
      batch = as.factor(sample_data$batch),
      biological_var = as.factor(sample_data$biological_var)
    )
    
    variance_explained <- summary(pca_result)$importance[2, 1:2] * 100
    
    return(list(pca_df = pca_df, variance_explained = variance_explained))
  }
  
  # Configuration summaries
  output$count_config_summary <- renderText({
    req(values$count_data, input$feature_col)
    
    dropped <- if(is.null(input$drop_cols)) "None" else paste(input$drop_cols, collapse = ", ")
    
    paste0("Feature column: ", input$feature_col, "\n",
           "Dropped columns: ", dropped, "\n",
           "Note: Only columns matching sample metadata will be retained as sample data")
  })
  
  output$sample_config_summary <- renderText({
    req(values$sample_data, input$batch_col, input$bio_col, input$sample_name_col)
    
    paste0("Sample name column: ", input$sample_name_col, "\n",
           "Batch column: ", input$batch_col, "\n",
           "Biological variable column: ", input$bio_col)
  })
  
  output$matching_summary <- renderText({
    req(values$matched_data)
    
    paste0("Successfully matched ", ncol(values$matched_data$count_matrix), " samples\n",
           "Features: ", nrow(values$matched_data$count_matrix), "\n",
           "Batches: ", length(unique(values$matched_data$sample_data$batch)), "\n",
           "Biological groups: ", length(unique(values$matched_data$sample_data$biological_var)))
  })
  
  # PCA plots
  create_pca_plots <- function(pca_result, title) {
    req(pca_result)
    
    pca_df <- pca_result$pca_df
    var_exp <- pca_result$variance_explained
    
    # Plot colored by biological variable
    p1 <- plot_ly(pca_df, x = ~PC1, y = ~PC2, color = ~biological_var,
                  text = ~sample_name, hovertemplate = "%{text}<extra></extra>",
                  type = "scatter", mode = "markers") %>%
      layout(title = paste(title, "- Colored by Biological Variable"),
             xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
             yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)")))
    
    # Plot colored by batch
    p2 <- plot_ly(pca_df, x = ~PC1, y = ~PC2, color = ~batch,
                  text = ~sample_name, hovertemplate = "%{text}<extra></extra>",
                  type = "scatter", mode = "markers") %>%
      layout(title = paste(title, "- Colored by Batch"),
             xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
             yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)")))
    
    return(subplot(p1, p2, nrows = 1, margin = 0.05))
  }
  
  output$pca_before <- renderPlotly({
    req(values$pca_before)
    create_pca_plots(values$pca_before, "Data before batch correction")
  })
  
  # Run batch correction
  observeEvent(input$run_combat, {
    req(values$matched_data, input$batch_method)
    
    method_name <- ifelse(input$batch_method == "combat_seq", "ComBat-Seq", "ComBat")
    showModal(modalDialog(paste("Running", method_name, "batch correction..."), footer = NULL))
    
    tryCatch({
      count_matrix <- values$matched_data$count_matrix
      sample_data <- values$matched_data$sample_data
      
      # Convert factors
      batch_factor <- as.factor(sample_data$batch)
      bio_factor <- as.factor(sample_data$biological_var)
      
      if(input$batch_method == "combat_seq") {
        print("=== ComBat-Seq Debug ===")
        print(paste("Original data range:", paste(range(count_matrix), collapse = " to ")))
        
        batch_numeric <- as.numeric(batch_factor)
        bio_numeric <- as.numeric(bio_factor)
        
        corrected_matrix <- ComBat_seq(counts = count_matrix,
                                       batch = batch_numeric,
                                       group = bio_numeric)
        
        print(paste("Corrected data range:", paste(range(corrected_matrix), collapse = " to ")))
        print(paste("Data identical?", identical(count_matrix, corrected_matrix)))
        
      } else if(input$batch_method == "combat_original") {
        # Metabolomics/microarray: Log transform then use ComBat
        # Check for non-positive values
        if(any(count_matrix <= 0, na.rm = TRUE)) {
          showNotification("Warning: Non-positive values detected. Adding small constant (0.01) before log transformation.", 
                           type = "warning")
          count_matrix[count_matrix <= 0] <- 0.01
        }
        
        # Log2 transform
        log_matrix <- log2(count_matrix)
        
        # Create model matrix for biological variable
        mod <- model.matrix(~bio_factor)
        
        # Use original ComBat
        corrected_log_matrix <- ComBat(dat = log_matrix, 
                                       batch = batch_factor, 
                                       mod = mod, 
                                       par.prior = TRUE)
        
        # Back-transform to original scale
        corrected_matrix <- 2^corrected_log_matrix
      }
      
      values$corrected_data <- corrected_matrix
      
      # Perform PCA on corrected data
      values$pca_after <- perform_pca(corrected_matrix, sample_data)
      
      removeModal()
      showNotification(paste(method_name, "batch correction completed successfully!"), type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error in", method_name, "batch correction:", e$message), type = "error")
    })
  })
  
  output$pca_after <- renderPlotly({
    req(values$pca_after)
    create_pca_plots(values$pca_after, "Data after batch correction")
  })
  
  # Download handlers
  output$download_corrected <- downloadHandler(
    filename = function() {
      paste0("batch_corrected_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$corrected_data, values$matched_data)
      
      # Add feature names back
      corrected_df <- as.data.frame(values$corrected_data)
      corrected_df <- tibble::rownames_to_column(corrected_df, var = "feature")
      
      # Add documentation header
      header <- c(
        "# Batch-corrected data",
        paste("# Method used:", ifelse(input$batch_method == "combat_seq", "ComBat-Seq", "ComBat")),
        "# Generated by Batch Corrector app",
        paste("# Date:", Sys.time()),
        "# IMPORTANT: Use this data for visualization and clustering only",
        "# For statistical analysis, consider using raw data with batch as covariate",
        ""
      )
      
      writeLines(header, file)
      write.csv(corrected_df, file, row.names = FALSE, append = TRUE)
    }
  )
  
  output$download_metadata <- downloadHandler(
    filename = function() {
      paste0("sample_metadata_with_batch_factor_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$matched_data)
      
      metadata_df <- values$matched_data$sample_data
      metadata_df$batch_correction_factor <- metadata_df$batch
      
      # Add documentation header
      header <- c(
        "# Sample metadata with batch correction factor",
        "# Generated by Batch Corrector app",
        paste("# Date:", Sys.time()),
        "# Usage in differential expression:",
        "# DESeq2: design = ~ batch_correction_factor + biological_var",
        "# limma-voom: design = model.matrix(~ batch_correction_factor + biological_var)",
        "# edgeR: design = model.matrix(~ batch_correction_factor + biological_var)",
        ""
      )
      
      writeLines(header, file)
      write.csv(metadata_df, file, row.names = FALSE, append = TRUE)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)