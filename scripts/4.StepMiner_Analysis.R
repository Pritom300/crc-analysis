

options(stringsAsFactors = FALSE)



expr_clean <- readRDS("../output/expr_clean.rds")
pheno_clean <- readRDS("../output/pheno_clean.rds")



library(dplyr)




fitStep <- function(data_vector) {
  
  if (any(is.na(data_vector)) || var(data_vector) == 0) {
    return(NULL)
  }
  
  n <- length(data_vector)
  sstot <- sum((data_vector - mean(data_vector))^2)
  
  best_sse <- Inf
  best_index <- NA
  
  for (i in 2:(n-2)) {
    
    m1 <- mean(data_vector[1:i])
    m2 <- mean(data_vector[(i+1):n])
    
    sse1 <- sum((data_vector[1:i] - m1)^2)
    sse2 <- sum((data_vector[(i+1):n] - m2)^2)
    
    sse <- sse1 + sse2
    
    if (sse < best_sse) {
      best_sse <- sse
      best_index <- i
    }
  }
  
  if (is.na(best_index)) return(NULL)
  
  m1 <- mean(data_vector[1:best_index])
  m2 <- mean(data_vector[(best_index+1):n])
  
  statistic <- ((sstot - best_sse) / 3) / (best_sse / (n - 4))
  
  return(list(index = best_index, statistic = statistic, m1 = m1, m2 = m2))
}


stage_order <- c(
  "Normal Colon",
  "Adenoma (Polyp)",
  "Primary Colorectal Cancer",
  "Metastasis (Liver/Lung)"
)

pheno_clean$Stage <- factor(pheno_clean$Stage, levels = stage_order)


ordered_samples <- rownames(pheno_clean)[order(pheno_clean$Stage)]

expr_ordered <- expr_clean[, ordered_samples]

cat("Samples ordered:", ncol(expr_ordered), "\n")





# LOADED MULTI-STAGE DEG 

deg1 <- read.csv("../output/LIMMA_Results/DEG_Adenoma_vs_Normal.csv")
deg2 <- read.csv("../output/LIMMA_Results/DEG_Cancer_vs_Adenoma.csv")
deg3 <- read.csv("../output/LIMMA_Results/DEG_Metastasis_vs_Cancer.csv")

all_deg <- unique(c(deg1$ProbeID, deg2$ProbeID, deg3$ProbeID))

cat(" Total DEG candidates:", length(all_deg), "\n")


# KEEPInG ONLY AVAILABLE GENES

common_genes <- intersect(all_deg, rownames(expr_ordered))

expr_deg <- expr_ordered[common_genes, ]

cat("Genes after matching:", nrow(expr_deg), "\n")


#  MONOTONIC FILTER for PROGRESSION SIGNAL

is_monotonic <- function(x) {
  all(diff(x) >= 0) || all(diff(x) <= 0)
}

# Uses stage mean ONLY for filtering (NOT for stepminer)
expr_stage_mean <- sapply(levels(pheno_clean$Stage), function(stage) {
  samples <- rownames(pheno_clean)[pheno_clean$Stage == stage]
  rowMeans(expr_clean[, samples, drop = FALSE])
})

expr_stage_mean <- expr_stage_mean[common_genes, ]

monotonic_genes <- apply(expr_stage_mean, 1, is_monotonic)

expr_final <- expr_deg[monotonic_genes, ]

cat(" Monotonic genes:", nrow(expr_final), "\n")


# RUNNING STEPMINER INSPIRED (SAMPLE LEVEL)

cat("\n Running StepMiner...\n")

step_results <- data.frame()

for (gene in rownames(expr_final)) {
  
  vec <- as.numeric(expr_final[gene, ])
  
  res <- fitStep(vec)
  
  if (is.null(res)) next
  
  direction <- if(res$m1 < res$m2) "Up-Step" else "Down-Step"
  
  step_results <- rbind(step_results, data.frame(
    ProbeID = gene,
    StepIndex = res$index,
    F_Stat = res$statistic,
    Direction = direction
  ))
}
step_results_all <- step_results

step_results <- step_results[order(-step_results$F_Stat), ]

# threshold 
step_results <- step_results[step_results$F_Stat > 10, ]

cat("\n Final StepMiner genes:", nrow(step_results), "\n")
print(head(step_results, 10))


# SAVE

write.csv(step_results, "../output/StepMiner_Q1_Final_Results.csv", row.names = FALSE)


# SAVE WITH SYMBOL 


# Added SYMBOL before saving
library(AnnotationDbi)
library(hgu133plus2.db)

symbols <- mapIds(
  hgu133plus2.db,
  keys = step_results$ProbeID,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

step_results$SYMBOL <- symbols
step_results <- step_results[!is.na(step_results$SYMBOL), ]

write.csv(step_results, "../output/StepMiner_Q1_Final_Results.csv", row.names = FALSE)


cat("\n Completed! \n")



# F-THRESHOLD SENSITIVITY CHECK 




thresholds <- c(5, 10, 15, 20, 25)
threshold_summary <- data.frame()

for (th in thresholds) {
  n_genes <- sum(step_results_all$F_Stat > th, na.rm = TRUE)
  threshold_summary <- rbind(threshold_summary, data.frame(F_Threshold = th, N_Genes = n_genes))
}

cat("\n F-THRESHOLD SENSITIVITY:\n")
print(threshold_summary)
write.csv(threshold_summary, "../output/F_Threshold_Sensitivity.csv", row.names = FALSE)





saveRDS(ordered_samples, "../output/ordered_samples.rds")
saveRDS(expr_final, "../output/expr_final.rds")
saveRDS(step_results, "../output/step_results.rds")

saveRDS(expr_clean, "../output/expr_clean.rds")
saveRDS(pheno_clean, "../output/pheno_clean.rds")


