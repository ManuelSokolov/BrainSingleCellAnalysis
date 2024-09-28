library(scater)
library(SingleR)
library(readr)
library(dplyr)


#setwd("../..")

setwd("/Users/manuelsokolovr/Desktop/BrainSingleCellAnalysis")

reference <- readRDS("data/hiPSC/rds/hiPSC_combined_ref.rds")

# Verify the result
unique(reference$CellClass)

names(assays(reference))=c("counts")

reference <- logNormCounts(reference)

query <- readRDS("data/hiPSC/rds/jovanovic_FCDI_Astro_unnormalized.rds")
names(assays(query))=c("counts")


query1 <- readRDS("data/hiPSC/rds/jovanovic_Glu_Neurons_unnormalized.rds")
names(assays(query1))=c("counts")

query2 <- readRDS("data/hiPSC/rds/jovanovic_WA09_Astro_unnormalized.rds")
names(assays(query2))=c("counts")


query <- logNormCounts(query)
query1 <- logNormCounts(query1)
query2 <- logNormCounts(query2)


predictions <- SingleR(test=query, assay.type.test=1, 
                       ref=reference, labels=reference$CellClass, de.method="wilcox", de.n = 263)

predictions1 <- SingleR(test=query1, assay.type.test=1, 
                        ref=reference, labels=reference$CellClass, de.method="wilcox", de.n = 263)


predictions2 <- SingleR(test=query2, assay.type.test=1, 
                        ref=reference, labels=reference$CellClass, de.method="wilcox", de.n = 263)

names(assays(query))=c("logcounts")

res_1 <- runUMAP(query)
res_1$prediction <- predictions$labels
# Plot UMAP with custom colors
plotUMAP(res_1, colour_by="prediction")


plotDeltaDistribution(predictions)
plotScoreHeatmap(predictions)

names(assays(query1))=c("logcounts")
res_1 <- runUMAP(query1)
res_1$prediction <- predictions1$labels
#Plot UMAP with custom colors
plotUMAP(res_1, colour_by="prediction")

plotDeltaDistribution(predictions1)
plotScoreHeatmap(predictions1)

names(assays(query2))=c("logcounts")
res_1 <- runUMAP(query2)
res_1$prediction <- predictions2$labels
#Plot UMAP with custom colors
plotUMAP(res_1, colour_by="prediction")

plotDeltaDistribution(predictions2)
plotScoreHeatmap(predictions2)


# Initialize an empty named list (dictionary) to store gene counts
gene_counts <- list()

# List of all cell types to iterate through
cell_types <- names(predictions@metadata$de.genes)

# Loop through each cell type
for (cell_type in cell_types) {
  
  # Get the genes for the current cell type - Change $ for the specific cell type
  genes <- predictions@metadata$de.genes$Glioblast[[cell_type]]
  
  # Loop through each gene in the current cell type
  for (gene in genes) {
    if (gene %in% names(gene_counts)) {
      # If the gene is already in the dictionary, increment its count
      gene_counts[[gene]] <- gene_counts[[gene]] + 1
    } else {
      # If the gene is not in the dictionary, add it with a count of 1
      gene_counts[[gene]] <- 1
    }
  }
}

# Filter the dictionary to get only those genes that occur 8 times
genes_with_8_occurrences_rg <- names(gene_counts)[sapply(gene_counts, function(x) x == 8)]
genes_with_8_occurrences_gl <- names(gene_counts)[sapply(gene_counts, function(x) x == 8)]

genes_with_8_occurrences_gl


plotUMAP(res_2, colour_by="COL1A2", text_colour="red")









