

library(hgu133plus2.db)
library(AnnotationDbi)
library(randomForest)


expr_train <- readRDS("../output/expr_train.rds")
expr_test <- readRDS("../output/expr_test.rds")

pheno_test <- readRDS("../output/pheno_test.rds")

expr_clean <- readRDS("../output/expr_clean.rds")      
pheno_clean <- readRDS("../output/pheno_clean.rds") 



collapse_to_symbol <- function(expr_mat) {
  if (max(expr_mat, na.rm = TRUE) > 50) {
    expr_mat <- log2(expr_mat + 1)
  }
  probe_to_symbol <- mapIds(
    hgu133plus2.db,
    keys = rownames(expr_mat),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  keep <- !is.na(probe_to_symbol)
  expr_mat <- expr_mat[keep, ]
  rownames(expr_mat) <- probe_to_symbol[keep]
  expr_mat <- aggregate(expr_mat, by = list(Gene = rownames(expr_mat)), FUN = max)
  rownames(expr_mat) <- expr_mat$Gene
  expr_mat$Gene <- NULL
  as.matrix(expr_mat)
}



cat("\n Building symbol-level matrix for TRAIN...\n")
expr_clean <- collapse_to_symbol(expr_train)
cat(" Train matrix:", dim(expr_clean), "\n")

cat("\n Building symbol-level matrix for TEST ...\n")
expr_test_sym <- collapse_to_symbol(expr_test)
cat(" Test matrix:", dim(expr_test_sym), "\n")


# STAGE CLASSIFICATION, EVALUATE ONCE ON HELD-OUT TEST


cat("\n", rep("=", 70), "\n")
cat(" STAGE CLASSIFICATION MODEL\n")
cat(rep("=", 70), "\n\n")

step_results <- read.csv("../output/StepMiner_Q8_Final_Results.csv")

get_top_genes <- function(transition_name, step_df, n = 25) {
  genes <- step_df$SYMBOL[step_df$Transition == transition_name]
  genes <- unique(genes[!is.na(genes)])
  gene_fstat <- step_df[step_df$SYMBOL %in% genes, c("SYMBOL", "F_Stat")]
  gene_fstat <- gene_fstat[!duplicated(gene_fstat$SYMBOL), ]
  gene_fstat <- gene_fstat[order(-gene_fstat$F_Stat), ]
  return(head(gene_fstat$SYMBOL, n))
}

selected_genes <- unique(c(
  get_top_genes("Normal Colon → Adenoma (Polyp)", step_results, 25),
  get_top_genes("Adenoma (Polyp) → Primary Colorectal Cancer", step_results, 25),
  get_top_genes("Primary Colorectal Cancer → Metastasis (Liver/Lung)", step_results, 25)
))
selected_genes <- selected_genes[!is.na(selected_genes)]
cat(" Selected", length(selected_genes), "balanced genes (discovered on TRAIN only)\n")

common_genes <- intersect(selected_genes, intersect(rownames(expr_clean), rownames(expr_test_sym)))
cat(" Common genes available in both train & test:", length(common_genes), "\n")
if(length(common_genes) < 10) stop(" Too few common genes!")

#  TRAIN matrix
X_train <- as.data.frame(t(expr_clean[common_genes, rownames(pheno_clean), drop = FALSE]))
y_train <- droplevels(as.factor(pheno_clean$Stage))
colnames(X_train) <- make.names(colnames(X_train))

# TEST matrix 
test_common_samples <- intersect(colnames(expr_test_sym), rownames(pheno_test))
X_test <- as.data.frame(t(expr_test_sym[common_genes, test_common_samples, drop = FALSE]))
y_test  <- factor(pheno_test[test_common_samples, "Stage"], levels = levels(y_train))
colnames(X_test) <- make.names(colnames(X_test))

cat("\n Train:", nrow(X_train), "| Test:", nrow(X_test), "\n")
cat(" Sample overlap check (must be 0):",
    length(intersect(rownames(pheno_clean), rownames(pheno_test))), "\n")


# TRAIN RANDOM FOREST


set.seed(123)
rf_model <- randomForest(
  x = X_train, y = y_train,
  ntree = 500, importance = TRUE, na.action = na.omit
)




# EVALUATE ONCE ON HELD-OUT TEST


predictions <- predict(rf_model, newdata = X_test)
confusion <- table(Predicted = predictions, Actual = y_test)

cat("\n CONFUSION MATRIX:\n")
print(confusion)

accuracy <- sum(diag(confusion)) / sum(confusion)
cat("\n OVERALL ACCURACY (truly held-out):", round(accuracy * 100, 2), "%\n")

classes <- levels(y_test)
cat("\n PER-CLASS PERFORMANCE:\n")
for(c in classes) {
  tp <- if(c %in% rownames(confusion) & c %in% colnames(confusion)) confusion[c, c] else 0
  fp <- if(c %in% rownames(confusion)) sum(confusion[c, ]) - tp else 0
  fn <- if(c %in% colnames(confusion)) sum(confusion[, c]) - tp else 0
  sens <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  prec <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  f1 <- ifelse(sens + prec > 0, 2 * sens * prec / (sens + prec), 0)
  cat(sprintf("   %-25s | Sensitivity: %.3f | Precision: %.3f | F1: %.3f\n", c, sens, prec, f1))
}




set.seed(123)


predictions <- predict(rf_model, newdata = X_test)
confusion   <- table(Predicted = predictions, Actual = y_test)

classes <- levels(y_test)                      
cm <- as.matrix(confusion)                      
cm <- cm[classes, classes]                     
N  <- sum(cm)

cp_ci <- function(x, n) {                       # exact Clopper-Pearson 95% CI
  if (n == 0) return(c(NA_real_, NA_real_))
  as.numeric(binom.test(x, n)$conf.int)
}

# Class-specific metrics (one-vs-rest) with exact CIs 
per_class <- do.call(rbind, lapply(seq_along(classes), function(i) {
  tp <- cm[i, i]; fn <- sum(cm[, i]) - tp
  fp <- sum(cm[i, ]) - tp; tn <- N - tp - fn - fp
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1   <- if (is.na(prec) || (prec + sens) == 0) 0 else 2 * prec * sens / (prec + sens)
  se <- cp_ci(tp, tp + fn); sp <- cp_ci(tn, tn + fp)
  data.frame(
    Class = classes[i], Support = tp + fn, TP = tp, FN = fn, FP = fp, TN = tn,
    Sensitivity = round(sens, 3), `Sensitivity 95% CI` = sprintf("%.3f-%.3f", se[1], se[2]),
    Specificity = round(spec, 3), `Specificity 95% CI` = sprintf("%.3f-%.3f", sp[1], sp[2]),
    Precision = round(prec, 3), F1 = round(f1, 3), check.names = FALSE
  )
}))

overall_acc <- sum(diag(cm)) / N
acc_ci      <- cp_ci(sum(diag(cm)), N)
raw <- sapply(seq_along(classes), function(i) {
  tp <- cm[i, i]; fn <- sum(cm[, i]) - tp; fp <- sum(cm[i, ]) - tp
  sens <- tp / (tp + fn); prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  c(sens = sens, f1 = if (is.na(prec) || (prec + sens) == 0) 0 else 2 * prec * sens / (prec + sens))
})
balanced_acc <- mean(raw["sens", ])
macro_f1     <- mean(raw["f1", ])



#  Bootstrap 95% CI for balanced accuracy and macro-F1 
true_v <- as.character(y_test); pred_v <- as.character(predictions)
boot_metrics <- function(t, p) {
  m <- table(factor(p, levels = classes), factor(t, levels = classes))
  n <- sum(m)
  sens <- sapply(seq_along(classes), function(i) m[i, i] / sum(m[, i]))
  prec <- sapply(seq_along(classes), function(i) m[i, i] / sum(m[i, ]))
  f1   <- ifelse(is.na(prec) | is.na(sens) | (prec + sens) == 0, 0, 2 * prec * sens / (prec + sens))
  c(bal = mean(sens, na.rm = TRUE), f1 = mean(f1, na.rm = TRUE))
}
B <- 2000
bs <- t(replicate(B, { idx <- sample(seq_len(length(true_v)), replace = TRUE)
boot_metrics(true_v[idx], pred_v[idx]) }))
bal_ci <- quantile(bs[, "bal"], c(0.025, 0.975), na.rm = TRUE)
f1_ci  <- quantile(bs[, "f1"],  c(0.025, 0.975), na.rm = TRUE)

# Ordinal-aware summaries (no new model; uses the same predictions) 
ti <- as.integer(factor(true_v, levels = classes))
pr <- as.integer(factor(pred_v, levels = classes))
errors <- ti != pr
adjacent_share <- if (sum(errors) > 0) mean(abs(ti - pr)[errors] == 1) else NA_real_
mean_abs_stage_error <- mean(abs(ti - pr))
k <- length(classes)
W <- outer(seq_len(k), seq_len(k), function(a, b) (a - b)^2 / (k - 1)^2)
O <- cm / N
E <- outer(rowSums(cm), colSums(cm)) / N^2
qwk <- 1 - sum(W * O) / sum(W * E)

# save 
cat("\n FULL CONFUSION MATRIX (rows = Predicted, cols = Actual) \n")
print(addmargins(cm))
cat("\n CLASS-SPECIFIC METRICS \n")
print(per_class, row.names = FALSE)
cat(sprintf("\nOverall accuracy: %.3f (exact 95%% CI %.3f-%.3f)\n", overall_acc, acc_ci[1], acc_ci[2]))
cat(sprintf("Balanced accuracy: %.3f (bootstrap 95%% CI %.3f-%.3f)\n", balanced_acc, bal_ci[1], bal_ci[2]))
cat(sprintf("Macro-F1: %.3f (bootstrap 95%% CI %.3f-%.3f)\n", macro_f1, f1_ci[1], f1_ci[2]))
cat(sprintf("\nOrdinal-aware: %.1f%% of errors are to an adjacent stage; mean absolute stage error = %.3f; quadratic weighted kappa = %.3f\n",
            100 * adjacent_share, mean_abs_stage_error, qwk))

write.csv(as.data.frame.matrix(addmargins(cm)), "../output/Stage_Full_Confusion_Matrix.csv")
write.csv(per_class, "../output/Stage_Per_Class_Metrics_CI.csv", row.names = FALSE)
write.csv(data.frame(
  Metric = c("Accuracy", "Balanced accuracy", "Macro-F1", "Adjacent-stage share of errors",
             "Mean absolute stage error", "Quadratic weighted kappa"),
  Value = round(c(overall_acc, balanced_acc, macro_f1, adjacent_share, mean_abs_stage_error, qwk), 3),
  CI_95 = c(sprintf("%.3f-%.3f", acc_ci[1], acc_ci[2]),
            sprintf("%.3f-%.3f", bal_ci[1], bal_ci[2]),
            sprintf("%.3f-%.3f", f1_ci[1], f1_ci[2]), "", "", "")
), "../output/Stage_Overall_Metrics.csv", row.names = FALSE)







# FEATURE IMPORTANCE 


importance_df <- data.frame(
  Gene = rownames(rf_model$importance),
  Gini_Importance = rf_model$importance[, "MeanDecreaseGini"],
  Permutation_Importance = rf_model$importance[, "MeanDecreaseAccuracy"]
)
importance_df <- importance_df[order(-importance_df$Gini_Importance), ]

# Check rank concordance between the two importance measures
importance_df$Gini_Rank <- rank(-importance_df$Gini_Importance)
importance_df$Permutation_Rank <- rank(-importance_df$Permutation_Importance)
rank_correlation <- cor(importance_df$Gini_Rank, importance_df$Permutation_Rank, method = "spearman")

cat("\n TOP 10 IMPORTANT GENES (Gini):\n")
print(head(importance_df[order(-importance_df$Gini_Importance), c("Gene","Gini_Importance","Permutation_Importance")], 10))

cat("\n TOP 10 IMPORTANT GENES (Permutation):\n")
print(head(importance_df[order(-importance_df$Permutation_Importance), c("Gene","Gini_Importance","Permutation_Importance")], 10))

cat("\n Spearman rank correlation between Gini and Permutation importance:",
    round(rank_correlation, 3), "\n")

write.csv(importance_df, "../output/Stage_Classification_Importance.csv", row.names = FALSE)


saveRDS(rf_model, "../output/Stage_Classification_Final.rds")


saveRDS(selected_genes, "../output/selected_genes.rds")
saveRDS(X_train, "../output/X_train.rds")
saveRDS(y_train, "../output/y_train.rds")
saveRDS(X_test, "../output/X_test.rds")
saveRDS(y_test, "../output/y_test.rds")
saveRDS(common_genes, "../output/common_genes.rds") 

saveRDS(expr_train, "../output/expr_train.rds")
saveRDS(expr_test, "../output/expr_test.rds")
saveRDS(pheno_test, "../output/pheno_test.rds")
saveRDS(expr_clean, "../output/expr_clean.rds")
saveRDS(pheno_clean, "../output/pheno_clean.rds")

saveRDS(predictions, "../output/predictions.rds")

saveRDS(accuracy, "../output/accuracy.rds")






#Figure 


library(ggplot2)
library(reshape2)


# Load confusion matrix

conf_mat <- as.data.frame(confusion)
colnames(conf_mat) <- c("Predicted", "Actual", "Freq")


# Ploting heatmap

p1 <- ggplot(conf_mat, aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Stage Classification Confusion Matrix",
    x = "Actual Stage",
    y = "Predicted Stage",
    fill = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p1)

# Save
ggsave(
  "../output/Figure4A_Confusion_Matrix_Heatmap.png",
  plot = p1,
  width = 10,
  height = 8,
  dpi = 300
)




# Figure: Top 15 Feature Importance Plot



# Load importance

importance_df <- read.csv("../output/Stage_Classification_Importance.csv")

top15 <- importance_df[order(-importance_df$Gini_Importance), ][1:15, ]

# Remove duplicate genes if any
top15 <- top15[!duplicated(top15$Gene), ]

top15$Gene <- factor(top15$Gene, levels = rev(top15$Gene))


# Ploting feature importance

p2 <- ggplot(top15, aes(x = Gene, y = Gini_Importance)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 15 Stage Classification Biomarkers",
    x = "Gene",
    y = "Mean Decrease Gini"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p2)

# Save
ggsave(
  "../output/Figure4B_Top15_Feature_Importance.png",
  plot = p2,
  width = 10,
  height = 8,
  dpi = 300
)

