



read_csv_file <- function(file_path) {
  if (file.exists(file_path)) {
    return(read.csv(file_path, stringsAsFactors = FALSE))
  } else {
    cat("Target file was not found at:", file_path, "\n")
    return(NULL)
  }
}

comparisons <- c("Adenoma_vs_Normal", "Cancer_vs_Adenoma", "Metastasis_vs_Cancer")



# LIMMA Differential Expression Check (Training Set)

cat("\nLIMMA Output Verification (Training Set) \n")

significant_probes_pool <- c()

for (comp in comparisons) {
  cat("\nRunning checks for comparison group:", comp, "\n")
  
  # Load the corresponding output matrices
  full_raw  <- read_csv_file(file.path("../output/LIMMA_Results", paste0("DEG_", comp, ".csv")))
  sig_raw   <- read_csv_file(file.path("../output/LIMMA_Results", paste0("DEG_SIGNIFICANT_", comp, ".csv")))
  mapped_q1 <- read_csv_file(file.path("../output/LIMMA_Results", paste0("DEG_", comp, "_Q1Mapped.csv")))
  
  # Probe-level metrics
  if (!is.null(full_raw) && !is.null(sig_raw)) {
    cat("  Total number of probes evaluated in layout: ", nrow(full_raw), "\n")
    cat("  Number of statistically significant probes: ", nrow(sig_raw), "\n")
    
   
    significant_probes_pool <- c(significant_probes_pool, sig_raw$ProbeID)
  }
  
  # Gene-level metrics after collapsing multiple probes
  if (!is.null(mapped_q1)) {
    sig_genes <- mapped_q1[mapped_q1$adj.P.Val < 0.05 & abs(mapped_q1$logFC) > 1, ]
    
    cat("  Total unique gene symbols mapped:           ", nrow(mapped_q1), "\n")
    cat("  Number of significant unique genes:         ", nrow(sig_genes), "\n")
  }
}

#  Calculate union of unique significant probes across all groups
if (length(significant_probes_pool) > 0) {
  total_unique_probes <- length(unique(significant_probes_pool))
  cat("\nGlobal Union: Total unique significant probes across all 3 comparisons:", total_unique_probes, "\n")
}



#  StepMiner Stability and Filtering Checks

cat("\nSection 2: StepMiner Calculations and Sensitivity Tests \n")

# Process within-stage permutation results
permutation_data <- read_csv_file("../output/StepMiner_Permutation_Stability.csv")
if (!is.null(permutation_data)) {
  cat("Total monotonic probes passing into the permutation step:", nrow(permutation_data), "\n")
  
  mean_transition_stability <- mean(permutation_data$Pct_Same_Transition) * 100
  cat("Mean percentage of permutations retaining the same transition label:", 
      round(mean_transition_stability, 1), "%\n")
}

# Process F-threshold sensitivity metrics
f_sensitivity_table <- read_csv_file("../output/F_Threshold_Sensitivity.csv")
if (!is.null(f_sensitivity_table)) {
  cat("\nProbe count distribution per F-statistic threshold tuning:\n")
  print(f_sensitivity_table)
}

# Process final filtered StepMiner targets (F > 10 and mapped to a gene symbol)
final_stepminer_results <- read_csv_file("../output/StepMiner_Q8_Final_Results.csv")
if (!is.null(final_stepminer_results)) {
  cat("\nFinal StepMiner selections passing standard filters (F > 10 + valid symbol):\n")
  cat("  Total probe count selected:    ", nrow(final_stepminer_results), "\n")
  cat("  Total unique gene symbol count:", length(unique(final_stepminer_results$SYMBOL)), "\n")
  
  cat("\nBreakdown by biological transition and directional step orientation:\n")
  print(table(final_stepminer_results$Transition, final_stepminer_results$Direction))
}



# Random Forest Biomarker Panel Check

cat("\n Section 3: Final Random Forest Feature Set Verification \n")

importance_metrics <- read_csv_file("../output/Stage_Classification_Importance.csv")
if (!is.null(importance_metrics)) {
  total_rf_features <- length(unique(importance_metrics$Gene))
  cat("Total unique features embedded in the finalized Random Forest biomarker panel:", total_rf_features, "\n")
}

cat("\ncompleted successfully.\n")


# Colon-Only Normal Sensitivity & OOB Error

library(GEOquery)
gse <- getGEO(filename="../Data/input/GSE41258_series_matrix.txt.gz", getGPL=FALSE)
raw <- pData(gse)
tcol <- names(raw)[sapply(raw, function(x) any(grepl("^tissue:", x)))][1]
tissue <- sub("^tissue: ", "", raw[[tcol]]); names(tissue) <- rownames(raw)

pheno_train <- readRDS("../output/pheno_train.rds")
pheno_test  <- readRDS("../output/pheno_test.rds")
y_test <- readRDS("../output/y_test.rds")
pred   <- readRDS("../output/predictions.rds")

cat("\nTRAIN, Normal stage by tissue:\n")
print(table(tissue[rownames(pheno_train)[pheno_train$Stage == "Normal Colon"]]))
cat("\nTEST, Normal stage by tissue:\n")
print(table(tissue[rownames(pheno_test)[pheno_test$Stage == "Normal Colon"]]))

ids <- names(pred); tis <- tissue[ids]; is_norm <- y_test == "Normal Colon"
cat("\nTEST normals: tissue x predicted class\n")
print(table(Tissue = tis[is_norm], Predicted = pred[is_norm]))

keep <- !(is_norm & tis != "Normal Colon")
cm2 <- table(Predicted = pred[keep], Actual = y_test[keep])
cat("\nTest samples kept (colon-only normals):", sum(keep), "\n"); print(cm2)
cat("Accuracy with colon-only normals:", round(sum(diag(cm2)) / sum(cm2), 3), "\n")

step <- read.csv("../output/StepMiner_Q8_Final_Results.csv")
p1 <- step$ProbeID[step$Transition == "Normal Colon → Adenoma (Polyp)"]
ex <- exprs(gse); if (max(ex, na.rm = TRUE) > 50) ex <- log2(ex + 1)
m  <- function(lbl) rowMeans(ex[p1, names(tissue)[tissue == lbl], drop = FALSE])
dl <- m("Normal Liver") - m("Normal Colon"); dg <- m("Normal Lung") - m("Normal Colon")
cat("\nNormal->Adenoma probes:", length(p1),
    "| |liver-colon| > 2:", sum(abs(dl) > 2), "| |lung-colon| > 2:", sum(abs(dg) > 2), "\n")



#OOB Error
rf_model <- readRDS("../output/Stage_Classification_Final.rds")

oob_error <- rf_model$err.rate[500, "OOB"] * 100

print(oob_error)

