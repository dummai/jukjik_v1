#!/usr/bin/env Rscript
#
# Drug Synergy Analysis using SynergyFinder Plus (R/Bioconductor)
#
# Calculates drug synergy scores using ZIP, Bliss, Loewe, and HSA models.
# Generates publication-quality dose-response plots and synergy heatmaps.
#
# Input files:
# - SynergyFinder_input.csv: Drug combination data
# - SynergyFinder_BlockID.csv: Block ID summary with experiment names
#
# Output:
# - synergy_summary.csv: Combined summary of all synergy scores
# - Dose-response plots (JPEG, 300 DPI)
# - Synergy heatmaps (JPEG, 300 DPI)
# - Detailed per-block synergy data
#

# ========== USER ADJUSTABLE PARAMETERS ==========
outlier_detection <- "non"        # Options: "non", "part", "all"
curve_fit_method <- "LL4"         # LL4 (4-param logistic)
dose_response_plot <- TRUE        # TRUE/FALSE
synergy_models <- c("ZIP", "Bliss", "Loewe", "HSA")
input_data_type <- "inhibition"   # "inhibition" or "viability"
# ================================================

# Install and load required packages
install_if_missing <- function(pkg, bioc = FALSE) {
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
  } else {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  library(pkg, character.only = TRUE)
}

cat("Loading required packages...\n")
install_if_missing("synergyfinder", bioc = TRUE)
install_if_missing("ggplot2")
install_if_missing("openxlsx", bioc = FALSE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  input_path <- args[1]
  block_id_path <- args[2]
  output_dir <- if (length(args) >= 3) args[3] else "/Users/kwan/Documents/dengue/DENV_AI/synergy_results"
} else {
  input_path <- "/Users/kwan/Documents/dengue/DENV_AI/SynergyFinder_input.csv"
  block_id_path <- "/Users/kwan/Documents/dengue/DENV_AI/SynergyFinder_BlockID.csv"
  output_dir <- "/Users/kwan/Documents/dengue/DENV_AI/synergy_results"
}

cat("\n========================================\n")
cat("Drug Synergy Analysis\n")
cat("Using SynergyFinder Plus\n")
cat("========================================\n\n")

cat("Parameters:\n")
cat(paste0("  Outlier detection: ", outlier_detection, "\n"))
cat(paste0("  Curve fit method: ", curve_fit_method, "\n"))
cat(paste0("  Dose response plot: ", dose_response_plot, "\n"))
cat(paste0("  Synergy models: ", paste(synergy_models, collapse = ", "), "\n"))
cat(paste0("  Data type: ", input_data_type, "\n\n"))

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "dose_response"), showWarnings = FALSE)
dir.create(file.path(output_dir, "heatmaps"), showWarnings = FALSE)
dir.create(file.path(output_dir, "data"), showWarnings = FALSE)

# Load input data
cat(paste0("Loading input data from: ", input_path, "\n"))
input_data <- read.csv(input_path, stringsAsFactors = FALSE)
cat(paste0("  Total rows: ", nrow(input_data), "\n"))

# Load block ID data
cat(paste0("Loading block ID data from: ", block_id_path, "\n"))
block_data <- read.csv(block_id_path, stringsAsFactors = FALSE)
cat(paste0("  Total blocks: ", nrow(block_data), "\n\n"))

# Rename columns to match synergyfinder format
colnames(input_data) <- c("block_id", "drug_row", "drug_col", 
                          "conc_r", "conc_c", "response", "conc_r_unit")
input_data$conc_c_unit <- input_data$conc_r_unit

# Summary results
summary_results <- data.frame()

# Process each block
for (i in 1:nrow(block_data)) {
  block_id <- block_data$block_id[i]
  drug1 <- block_data$drug1[i]
  drug2 <- block_data$drug2[i]
  experiment <- block_data$experiment[i]
  
  cat(paste0("\n", paste(rep("=", 60), collapse = ""), "\n"))
  cat(paste0("Processing Block ", block_id, ": ", drug1, " + ", drug2, 
             " (", experiment, ")\n"))
  cat(paste0(paste(rep("=", 60), collapse = ""), "\n"))
  
  # Filter data for this block
  block_input <- input_data[input_data$block_id == block_id, ]
  
  if (nrow(block_input) == 0) {
    cat(paste0("  Warning: No data found for block ", block_id, "\n"))
    next
  }
  
  # Count data types
  n_drug1_single <- sum(block_input$conc_r > 0 & block_input$conc_c == 0)
  n_drug2_single <- sum(block_input$conc_r == 0 & block_input$conc_c > 0)
  n_combinations <- sum(block_input$conc_r > 0 & block_input$conc_c > 0)
  n_control <- sum(block_input$conc_r == 0 & block_input$conc_c == 0)
  
  cat(paste0("  Drug1 single doses: ", n_drug1_single, "\n"))
  cat(paste0("  Drug2 single doses: ", n_drug2_single, "\n"))
  cat(paste0("  Combinations: ", n_combinations, "\n"))
  cat(paste0("  Control: ", n_control, "\n"))
  
  tryCatch({
    # Reshape data for synergyfinder
    cat("\n  Reshaping data...\n")
    data_processed <- ReshapeData(
      data = block_input,
      data_type = input_data_type,
      impute = TRUE,
      noise = FALSE
    )
    
    # Calculate synergy scores
    cat("  Calculating synergy scores...\n")
    data_with_synergy <- CalculateSynergy(
      data = data_processed,
      method = synergy_models,
      correct_baseline = outlier_detection,
      adjusted = TRUE
    )
    
    # Extract synergy scores
    scores <- data_with_synergy$scores
    
    # Calculate summary statistics for each model
    for (model in synergy_models) {
      score_col <- paste0(model, "_synergy")
      if (score_col %in% colnames(scores)) {
        model_scores <- scores[[score_col]]
        model_scores <- model_scores[!is.na(model_scores)]
        
        if (length(model_scores) > 0) {
          mean_syn <- mean(model_scores, na.rm = TRUE)
          std_syn <- sd(model_scores, na.rm = TRUE)
          min_syn <- min(model_scores, na.rm = TRUE)
          max_syn <- max(model_scores, na.rm = TRUE)
          
          summary_results <- rbind(summary_results, data.frame(
            block_id = block_id,
            experiment = experiment,
            drug1 = drug1,
            drug2 = drug2,
            model = model,
            mean_synergy = mean_syn,
            std_synergy = std_syn,
            min_synergy = min_syn,
            max_synergy = max_syn,
            n_combinations = length(model_scores)
          ))
          
          cat(paste0("    ", model, ": Mean = ", round(mean_syn, 2), 
                     ", Std = ", round(std_syn, 2), "\n"))
        }
      }
    }
    
    # Generate dose-response plots
    if (dose_response_plot) {
      cat("\n  Generating dose-response plots...\n")
      
      dose_path <- file.path(output_dir, "dose_response", 
                              paste0("block", block_id, "_dose_response.jpg"))
      
      tryCatch({
        p <- PlotDoseResponse(
          data = data_with_synergy,
          block_ids = block_id,
          save_file = FALSE
        )
        
        ggsave(
          filename = dose_path,
          plot = p,
          device = "jpeg",
          width = 12,
          height = 6,
          units = "in",
          dpi = 300
        )
        cat(paste0("    Saved: ", dose_path, "\n"))
      }, error = function(e) {
        cat(paste0("    Warning: Could not generate dose-response plot: ", e$message, "\n"))
      })
    }
    
    # Generate synergy heatmaps
    cat("\n  Generating synergy heatmaps...\n")
    
    for (model in synergy_models) {
      heatmap_path <- file.path(output_dir, "heatmaps",
                                  paste0("block", block_id, "_", model, ".jpg"))
      
      tryCatch({
        p <- PlotSynergy(
          data = data_with_synergy,
          type = "heatmap",
          method = model,
          block_ids = block_id,
          save_file = FALSE,
          high_value_color = "#2166AC",
          low_value_color = "#B2182B"
        )
        
        ggsave(
          filename = heatmap_path,
          plot = p,
          device = "jpeg",
          width = 8,
          height = 6,
          units = "in",
          dpi = 300
        )
        cat(paste0("    Saved: ", heatmap_path, "\n"))
      }, error = function(e) {
        cat(paste0("    Warning: Could not generate ", model, " heatmap: ", e$message, "\n"))
      })
    }
    
    # Save block synergy data
    block_output_path <- file.path(output_dir, "data", 
                                    paste0("block", block_id, "_synergy.csv"))
    write.csv(scores, block_output_path, row.names = FALSE)
    cat(paste0("\n  Saved block data to: ", block_output_path, "\n"))
    
  }, error = function(e) {
    cat(paste0("  Error processing block ", block_id, ": ", e$message, "\n"))
  })
}

# Save summary results
if (nrow(summary_results) > 0) {
  cat("\n========================================\n")
  cat("Saving Summary Results\n")
  cat("========================================\n")
  
  # Save CSV summary
  summary_path <- file.path(output_dir, "synergy_summary.csv")
  write.csv(summary_results, summary_path, row.names = FALSE)
  cat(paste0("Saved synergy summary to: ", summary_path, "\n"))
  
  # Pivot to wide format
  summary_wide <- reshape(summary_results, 
                            idvar = c("block_id", "experiment", "drug1", "drug2"),
                            timevar = "model",
                            direction = "wide")
  
  # Save Excel summary
  excel_path <- file.path(output_dir, "synergy_summary.xlsx")
  write.xlsx(summary_wide, excel_path, rowNames = FALSE)
  cat(paste0("Saved Excel summary to: ", excel_path, "\n"))
}

cat(paste0("\n", paste(rep("=", 60), collapse = ""), "\n"))
cat("Synergy analysis complete!\n")
cat(paste0("Output directory: ", output_dir, "\n"))
cat(paste0(paste(rep("=", 60), collapse = ""), "\n"))
