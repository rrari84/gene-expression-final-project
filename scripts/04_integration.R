#!/usr/bin/env Rscript
# 04_integration.R
# ---------------------------------------------------------------------
# Integrate DEG results from multiple datasets to identify robust genes
# - Should be run from the project root directory
# - Assumes you have already run:
#    00_setup_environment.R (to install packages)
#    01_preprocess_microarray_from_CEL.R (optional, if using raw CEL files)
#    02_download_data.R (to get GEO datasets)
#    03_deg_analysis.R (to generate DEG results)
#
# Usage (from shell, in project root):
#   Rscript scripts/04_integration.R
#
# This script will combine DEG results from results/deg/ and identify
# genes that are consistently up- or down-regulated across datasets



library(tidyverse)
library(data.table)

# 1. Setup paths

# Input directory with per-dataset DEG results
deg_dir <- "results/deg"

# Output directory for integration results
int_dir <- "results/integration"

# Create output directory if it doesn't exist
if (!dir.exists(int_dir)) {
  dir.create(int_dir, recursive = TRUE)
  message("Created output directory: ", int_dir)
}

# 2. Load all DEG results

message("Loading DEG results from: ", deg_dir)

# Get all DEG result files
deg_files <- list.files(deg_dir, pattern = "_deg_results\\.(csv|rds)$",
                        full.names = TRUE)

if (length(deg_files) == 0) {
  stop("No DEG result files found in ", deg_dir)
}

message("Found ", length(deg_files), " dataset(s) to integrate")

# Load all datasets into a list
deg_list <- list()

for (file in deg_files) {
  # Extract dataset ID from filename
  dataset_id <- basename(file) %>%
    str_remove("_deg_results\\.(csv|rds)$")
  
  message(" Loading: ", dataset_id)
  
  # Load data (handle both CSV and RDS)
  if (str_ends(file, ".csv")) {
    deg_data <- read_csv(file, show_col_types = FALSE)
  } else {
    deg_data <- readRDS(file)
  }
  
  # Add dataset column if not present
  if(!"dataset" %in% colnames(deg_data)) {
    deg_data$dataset <- dataset_id
  }
  
  deg_list[[dataset_id]] <- deg_data
}

# Combine all datasets
all_deg <- bind_rows(deg_list)

message("Total rows across datasets: ", nrow(all_deg))

# 3. Standardize column names (if needed)

# Check for common column name variations
if("padj" %in% colnames(all_deg)) {
  # Already standardized
} else if("adj.P.Val" %in% colnames(all_deg)) {
  all_deg <- all_deg %>% rename(padj = adj.P.Val)
} else if ("FDR" %in% colnames(all_deg)) {
  all_deg <- all_deg %>% rename(padj = FDR)
}

if(!"log2FC" %in% colnames(all_deg) && "logFC" %in% colnames(all_deg)) {
  all_deg <- all_deg %>% rename(log2FC = logFC)
}

# 4. Define significant genes per dataset

# Set thresholds
padj_threshold <-0.05
log2fc_threshold <- 1

message("\nApplying significance thresholds:")
message("  padj < ", padj_threshold)
message("  |log2FC| > ", log2fc_threshold)

# Classify genes as up/down/not significant in each dataset
all_deg <- all_deg %>%
  mutate(
    direction = case_when(
      padj < padj_threshold & log2FC > log2fc_threshold ~ "up",
      padj < padj_threshold & log2FC < -log2fc_threshold ~ "down",
      TRUE ~ "not_sig"
    )
  )

# Summary statistics
sig_summary <- all_deg %>%
  filter(direction != "not_sig") %>%
  group_by(dataset, direction) %>%
  summarise(n_genes = n(), .groups = "drop")

print(sig_summary)

# 5. Count support across datasets

message("\nCounting gene support across datasets...")

# For each gene, count:
# - How many datasets show it as significantly UP
# - How many datasets show it as significantly DOWN
# - Total datasets tested

gene_support <- all_deg %>%
  group_by(gene) %>%
  summarise(
    n_datasets_total = n_distinct(dataset),
    n_up = sum(direction == "up", na.rm = TRUE),
    n_down = sum(direction == "down", na.rm = TRUE),
    n_not_sig = sum(direction == "not_sig", na.rm = TRUE),
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    median_log2FC = median(log2FC, na.rm = TRUE),
    mean_padj = mean(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Determine overall direction (must be consistent)
    robust_direction = case_when(
      n_up > 0 & n_down == 0 ~ "up",
      n_down > 0 & n_up == 0 ~ "down",
      n_up > 0 & n_down > 0 ~ "inconsistent",
      TRUE ~ "not_sig"
    ),
    # Support count (for consistent direction only)
    n_support = case_when(
      robust_direction == "up" ~ n_up,
      robust_direction == "down" ~ n_down,
      TRUE ~ 0
    )
  ) %>%
  arrange(desc(n_support), desc(abs(mean_log2FC)))

# 6. Define robust DEGs

# Set minimum number of datasets for "robust" classification
min_datasets_robust <-2 # Adjust based on how many datasets you have

message("\nDefining robust DEGs:")
message("  Minimum datasets required: ", min_datasets_robust)

# Extract robust UP genes
robust_up <- gene_support %>%
  filter(robust_direction == "up", n_support >= min_datasets_robust) %>%
  select(gene, n_datasets_supporting = n_support, n_datasets_total,
         mean_log2FC, median_log2FC, mean_padj) %>%
  arrange(desc(n_datasets_supporting), desc(mean_log2FC))

message("  Robust UP-regulated genes: ", nrow(robust_up))

# Extract robust DOWN genes
robust_down <- gene_support %>%
  filter(robust_direction == "down", n_support >= min_datasets_robust) %>%
  select(gene, n_datasets_supporting = n_support, n_datasets_total,
         mean_log2FC, median_log2FC, mean_padj) %>%
  arrange(desc(n_datasets_supporting), mean_log2FC)

message("  Robust DOWN-regulated genes: ", nrow(robust_down))

# 7. Optional: Fisher's method for meta analysis p-values

# Fisher's method combines p-values: -2 * sum(log(p_i))
# This follows chi-squared distributionwith 2k degrees of freedom

message("\nComputing meta-analysis p-values using Fisher's method...")

fisher_meta <- all_deg %>%
  filter(direction != "not_sig") %>%
  group_by(gene) %>%
  summarise(
    n_datasets = n(),
    fisher_stat = -2 * sum(log(pvalue), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    df = 2 * n_datasets,
    fisher_pvalue = pchisq(fisher_stat, df = df, lower.tail = FALSE)
  ) %>%
  arrange(fisher_pvalue)

# Add Fisher p-values to robust gene lists
robust_up <- robust_up %>%
  left_join(fisher_meta %>% select(gene, fisher_pvalue), by = "gene")

robust_down <- robust_down %>%
  left_join(fisher_meta %>% select(gene, fisher_pvalue), by = "gene")

# 8. Save results

message("\nSaving results to: ", int_dir)

# Save robust UP genes
write_csv(robust_up, file.path(int_dir, "robust_up_genes.csv"))
message("  Saved: robust_up_genes.csv (", nrow(robust_up), " genes)")

# Save robust DOWN genes
write_csv(robust_down, file.path(int_dir, "robust_down_genes.csv"))
message ("  Saved: robust_down_genes.csv (", nrow(robust_down)," genes)")

# Save full gene support summary
write_csv(gene_support, file.path(int_dir, "gene_support_summary.csv"))
message ("  Saved: gene_support_summary.csv (", nrow(gene_support)," genes)")

# Save combined DEG results
write_csv(all_deg, file.path(int_dir, "all_datasets_deg_combined.csv"))
message ("  Saved: all_datasets_deg_combined.csv")

# 9. Generate summary report
message("\n", paste(rep("=", 70), collapse = ""))
message("INTEGRATION SUMMARY")
message(paste(rep("=", 70), collapse = ""))
message("Total datasets analyzed: ", length(deg_files))
message("Total unique genes tested: ", n_distinct(all_deg$gene))
message("\nRobust DEGs (present in >=", min_datasets_robust, " datasets):")
message("  UP-regulated: ", nrow(robust_up))
message("  DOWN-regulated: ", nrow(robust_down))
message("  Total robust: ", nrow(robust_up) + nrow(robust_down))
message("\nTop 10 most robust UP-regulated genes:")
print(robust_up %>% head(10) %>% select(gene, n_datasets_supporting, mean_log2FC))
message("\nTop 10 most robust DOWN-regulated genes:")
print(robust_down %>% head(10) %>% select(gene, n_datasets_supporting, mean_log2FC))
message(paste(rep("=", 70), collapse = ""))

message("\nIntegration complete. Results saved to ", int_dir)
message("Output files saved to: ", int_dir)
message("\nNext: Run 05_enrichment.R for pathway analysis")


    
    