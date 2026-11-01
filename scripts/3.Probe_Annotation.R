

library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)

map_probes_q1 <- function(file_path) {
  
  cat("\n Processing:", file_path, "\n")
  
  if (!file.exists(file_path)) {
    cat(" File not found\n")
    return(NULL)
  }
  
  deg <- read.csv(file_path, stringsAsFactors = FALSE)
  

  if(!"ProbeID" %in% colnames(deg)){
    stop(" ProbeID column missing!")
  }
  
  deg$ProbeID <- as.character(deg$ProbeID)
  

  # Mapping

  mapping <- AnnotationDbi::select(
    hgu133plus2.db,
    keys = deg$ProbeID,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "PROBEID"
  )
  
  # Merge
  deg_mapped <- merge(deg, mapping,
                      by.x = "ProbeID",
                      by.y = "PROBEID",
                      all.x = TRUE)
  

  #  REMOVE NA GENES 

  deg_mapped <- deg_mapped[!is.na(deg_mapped$SYMBOL), ]
  
  
  # Keeping probe with highest 

  deg_mapped <- deg_mapped %>%
    group_by(SYMBOL) %>%
    slice_max(order_by = abs(logFC), n = 1) %>%
    ungroup()
  

  deg_mapped <- deg_mapped[order(deg_mapped$adj.P.Val), ]
  

  # Save
  
  new_name <- sub("\\.csv$", "_Q1Mapped.csv", file_path)
  write.csv(deg_mapped, new_name, row.names = FALSE)
  

  cat(" Final genes:", nrow(deg_mapped), "\n")
  cat(" Unique genes:", length(unique(deg_mapped$SYMBOL)), "\n")
  
  return(deg_mapped)
}


# main RUN



res1 <- map_probes_q1("../output/LIMMA_Results/DEG_Adenoma_vs_Normal.csv")
res2 <- map_probes_q1("../output/LIMMA_Results/DEG_Cancer_vs_Adenoma.csv")
res3 <- map_probes_q1("../output/LIMMA_Results/DEG_Metastasis_vs_Cancer.csv")

cat("\n gene lists created!\n")




