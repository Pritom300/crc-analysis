



# EXTERNAL VALIDATION ON GSE14333 


library(GEOquery)
library(randomForest)
library(hgu133plus2.db)
library(AnnotationDbi)
library(ggplot2)

rf_model <- readRDS("../output/Stage_Classification_Final.rds")
X_train <- readRDS("../output/X_train.rds")

cat("\n", rep("=", 70), "\n")
cat(" EXTERNAL VALIDATION ON GSE14333\n")
cat(rep("=", 70), "\n\n")






gse14333_file <- "../Data/input/GSE14333_series_matrix.txt.gz"

if(!file.exists(gse14333_file)) {
  cat("️ GSE14333 series matrix not found. Downloading from GEO...\n")
  gse14333 <- getGEO("GSE14333", destdir = "../output/", getGPL = FALSE)
  expr_gse14333 <- exprs(gse14333)
} else {
  gse14333 <- getGEO(filename = gse14333_file, getGPL = FALSE)
  expr_gse14333 <- exprs(gse14333)
}

cat(" Expression matrix dimensions:", dim(expr_gse14333), "\n")



# MAP PROBE ID TO GENE SYMBOL FOR GSE14333



probe_to_symbol_ext <- mapIds(
  hgu133plus2.db,
  keys = rownames(expr_gse14333),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

keep_ext <- !is.na(probe_to_symbol_ext)
expr_gse14333 <- expr_gse14333[keep_ext, ]
rownames(expr_gse14333) <- probe_to_symbol_ext[keep_ext]


expr_gse14333 <- aggregate(expr_gse14333, by = list(Gene = rownames(expr_gse14333)), FUN = max)
rownames(expr_gse14333) <- expr_gse14333$Gene
expr_gse14333$Gene <- NULL

cat(" After mapping:", dim(expr_gse14333), "\n")


# LOADInG CLINICAL DATA FOR GSE14333




clinical_path <- "../Data/input/GSE14333_clinical_data.csv"
clin_gse14333 <- read.csv(clinical_path, stringsAsFactors = FALSE)

# Extract GSM IDs
extract_gsm <- function(x) {
  match <- regexpr("GSM[0-9]+", x)
  if(match == -1) return(NA)
  return(substr(x, match, match + attr(match, "match.length") - 1))
}

clin_gse14333$Sample_ID <- sapply(clin_gse14333$Sample_geo_accession, extract_gsm)
clin_gse14333 <- clin_gse14333[!is.na(clin_gse14333$Sample_ID), ]




# GSE14333 has primary tumors only
clin_gse14333$Stage <- "Primary Colorectal Cancer"

cat(" Clinical samples:", nrow(clin_gse14333), "\n")


# MATCH SAMPLES AND GENES



# Common samples
common_samples_ext <- intersect(colnames(expr_gse14333), clin_gse14333$Sample_ID)
cat("   Common samples:", length(common_samples_ext), "\n")

# Common genes with our model
model_genes <- rownames(rf_model$importance)
common_genes_ext <- intersect(model_genes, rownames(expr_gse14333))
cat("   Common genes:", length(common_genes_ext), "\n")

if(length(common_genes_ext) < 10) {
  stop(" Too few common genes for validation!")
}

# Prepare validation data
X_valid <- t(expr_gse14333[common_genes_ext, common_samples_ext])
y_valid <- clin_gse14333[match(common_samples_ext, clin_gse14333$Sample_ID), "Stage"]
y_valid <- as.factor(y_valid)

cat("Validation set:", nrow(X_valid), "samples,", ncol(X_valid), "genes\n")


# PREDICT USING OUR TRAINED MODEL




# Ensure column names match
X_valid_df <- as.data.frame(X_valid)
colnames(X_valid_df) <- make.names(colnames(X_valid_df))

# Align to the exact feature set the model was trained on
train_genes <- colnames(X_train)
missing_genes <- setdiff(train_genes, colnames(X_valid_df))

if (length(missing_genes) > 0) {
  cat("️ Missing in GSE14333, imputing with train-set mean:", paste(missing_genes, collapse=", "), "\n")
  for (g in missing_genes) {
    X_valid_df[[g]] <- mean(X_train[[g]], na.rm = TRUE)
  }
}

X_valid_df <- X_valid_df[, train_genes, drop = FALSE]

# Predict
set.seed(123)
predictions <- predict(rf_model, newdata = X_valid_df)

# Since GSE14333 only has Primary Cancer, we can only validate that class
cat("\n Prediction distribution on GSE14333:\n")
print(table(predictions))

# Calculate accuracy (only for Primary Cancer class)
accuracy_primary <- mean(predictions == "Primary Colorectal Cancer")
cat("\n Accuracy for Primary Cancer class:", round(accuracy_primary * 100, 2), "%\n")


# SAVE 


validation_results <- data.frame(
  Dataset = "GSE14333",
  Samples = nrow(X_valid),
  Genes_Used = length(common_genes_ext),
  Primary_Cancer_Accuracy = round(accuracy_primary * 100, 2),
  Prediction_Distribution = paste(names(table(predictions)), table(predictions), collapse="; ")
)

write.csv(validation_results, "../output/External_Validation_GSE14333_Results.csv", row.names = FALSE)

cat("\n Saved: External_Validation_GSE14333_Results.csv\n")

cat("\n COMPLETE! \n")




#Figure

pred_dist <- as.data.frame(table(predictions))
colnames(pred_dist) <- c("Predicted_Class", "Count")
p6 <- ggplot(pred_dist, aes(x = reorder(Predicted_Class, -Count), y = Count, fill = Predicted_Class)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +
  theme_minimal(base_size = 12) +
  labs(title = "GSE14333 Prediction Distribution (n=290, all true Primary CRC)",
       x = "Predicted Class", y = "Count") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 20, hjust = 1))
print(p6)
ggsave("../output/Figure6_GSE14333_Prediction_Distribution.png", p6, width = 7, height = 6, dpi = 300)


#Figure 




saveRDS(accuracy_primary, "../output/accuracy_primary.rds")
saveRDS(predictions, "../output/predictions_GSE41333.rds")






