#### Loading libraries ####

# tidyverse: collection of R packages for data manipulation, visualization and modeling 
library("tidyverse")

# DESeq2: package for differential gene expression analysis
library("DESeq2")

# pheatmap: package for creating heatmaps, which will be used to visualise the results
#install.packages("pheatmap") # To install the package missing in the current RStudio env
library("pheatmap")

# RColorBrewer: package for creating color palettes, which will be used to customise the heatmaps
library("RColorBrewer")

# ggrepel: package that provides geoms for ggplot2 to repel overlapping text labels in the plots
library("ggrepel")

# readr: package to read .tsv file
library(readr)

# Import the sample sheet with meta-information about the experiment
sampletable = read.csv("metadata_cross04.csv", header=TRUE)

# add the sample names as row names (it is needed for some of the DESeq functions)
rownames(sampletable) = sampletable$sample
sampletable$condition = as.factor(sampletable$condition)

# Import the count table from STAR
count_matrix <- read.delim("cross4_counts_matrix_header.tsv", header=T, sep="\t", row.names=1)

#### Check that sample names match in both files ####
all(colnames(count_matrix) %in% rownames(sampletable)) # Must be TRUE
all(colnames(count_matrix) == rownames(sampletable)) # Must be TRUE

#make star matrix
se_star_matrix <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sampletable, design = ~ condition)

# FROM HERE IT IS THE NF-CORE TUTORIAL (https://nf-co.re/rnaseq/3.18.0/docs/usage/differential_expression_analysis/de_rstudio)

# PRE-FILTERING base on NF-CORE TUTORIAL

# Number of genes before filtering:
nrow(se_star_matrix)

# Select a minimal number of samples = 3
smallestGroupSize <- 3

# Select genes with a sum counts of at least 10 in 3 samples
keep <- rowSums(counts(se_star_matrix) >= 10) >= smallestGroupSize

# Keep only the genes that pass the threshold
se_star_matrix <- se_star_matrix[keep,]

# Number of genes left after low-count filtering:
nrow(se_star_matrix)

# Fit statistical model
dds_final <- DESeq(se_star_matrix)


#### Transform normalised counts for data visualisation ####
# A user can choose among vst and rlog. In this tutorial we will work with rlog transformed data
rld <- rlog(dds_final, blind = TRUE)


#### PLOT PCA ####
pca_plot <- plotPCA(rld, intgroup = "condition")

# Save the plot
ggsave("de_results/pca_plot.png", plot = pca_plot, width = 6, height = 5, dpi = 300)


#### HIERARCHICAL CLUSTERING ####

# Extract the matrix of rlog-transformed counts from the rld object
sampleDists <- dist(t(assay(rld)))  # Calculate pairwise distances between samples using 
                                    #the dist() function with Euclidean distance as the default 
                                    #method. By transposing the matrix with t(), we ensure that 
                                    #samples become rows and genes become columns, so that the 
                                    #dist function computes pairwise distances between samples.

# Convert distances to a matrix
sampleDistMatrix <- as.matrix(sampleDists)

# Set the row and column names of the distance matrix
rownames(sampleDistMatrix) <- paste(rld$condition, rld$sample, sep = "_")
colnames(sampleDistMatrix) <- paste(rld$condition, rld$sample, sep = "_")

# Define a color palette for the heatmap
colors <- colorRampPalette(rev(brewer.pal(9, "Greens")))(255) # function from RColorBrewer package

# Create the heatmap
clustering_plot <- pheatmap(sampleDistMatrix,
                            clustering_distance_rows = sampleDists,
                            clustering_distance_cols = sampleDists,
                            col = colors,
                            fontsize_col = 8,
                            fontsize_row = 8)

# Save the plot
ggsave("de_results/clustering_plot.png", plot = clustering_plot, width = 6, height = 5, dpi = 300)





#### INSPECT THE NORMALISED COUNTS ####
# Display the first few rows of the normalised counts to inspect the data
head(counts(dds_final, normalized = TRUE))

# Display the first few rows of the raw counts (not normalised) to compare with the normalised counts
head(counts(dds_final))

# Convert the normalised counts from the DESeq2 object to a tibble
normalised_counts <- as_tibble(counts(dds_final, normalized = TRUE))

# Add a column for gene names to the normalised counts tibble
normalised_counts$gene <- rownames(counts(dds_final))

# Relocate the gene column to the first position
normalised_counts <- normalised_counts %>%
  relocate(gene, .before = CROSS3_CONTROL_REP1) #change CROSS3_CONTROL_REP1 according to the first sample in your count matrix!

# Save the normalised counts
write.csv(normalised_counts, file = "de_results/normalised_counts.csv")





#### EXTRACT RESULTS ####

# DESeq2 function to extract the name of the contrast
resultsNames(dds_final) #if it is not the correct contrast you need to reset it

# Command to set the contrast, if necessary
#res <- results(dds_final, contrast = c("design_formula", "condition_of_interest", "reference_condition"), pAdjustMethod = "BH")
res <- results(dds_final, contrast = c("condition", "DS", "WW"), alpha = 0.05, pAdjustMethod = "BH")

# Visualise the results
head(res)

# Summarise the results showing the number of tested genes (genes with non-zero total read count), the genes up- and down-regulated at the selected threshold (alpha) and the number of genes excluded by the multiple testing due to a low mean count
summary(res)

# Store the res object inside another variable because the original res file will be required for other functions
res_viz <- res

# Add gene names as a new column to the results table
res_viz$gene <- rownames(res)

# Convert the results to a tibble for easier manipulation and relocate the gene column to the first position
res_viz <- as_tibble(res_viz) %>%
  relocate(gene, .before = baseMean)

# Save the results table
write.csv(res_viz, file = "de_results/de_result_table.csv")





#### EXTRACT SIGNIFICANT DE GENES ####

# Filter the results to include only significantly DE genes with a padj less than 0.05 and a log2FoldChange of at least 1 or -1
resSig <- subset(res_viz, padj < 0.05 & abs(log2FoldChange) > 1)

# Convert the results to a tibble for easier manipulation and relocate the gene column to the first position
resSig <- as_tibble(resSig) %>%
  relocate(gene, .before = baseMean)

# Order the significant genes by their adjusted p-value (padj) in ascending order
resSig <- resSig[order(resSig$padj),]

# Display the final results table of significant genes
resSig

# Save the significant DE genes
write.csv(resSig, file = "de_results/sig_de_genes.csv")


