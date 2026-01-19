#!/usr/bin/env Rscript

# The script will first perform differential binding analysis based on a `DESeq2` method (`edgeR` is also available as an option). 
# A minimum of two replicates per condition or factor are necessary for differential analysis. 
# Once differential binding analysis is complete, a correlation matrix and scatter XY plot in `.pdf` format will be generated, along with a final report in `.csv` format. 

# As a final output, an `.RData` data object will also be generated that will contain all the necessary information and can be loaded into RStudio without the `sampleSheet` or final report.

# DiffBind package is required. 

# If install.packages("DiffBind") command fails in newer versions of R, uncomment and run the commands below:
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DiffBind")

# set working environment
working_dir <- "/home/tarta/Lavoro/NIND_2/diffBind_E_nindensis/FHvsRH1/"
sampleSheet <- "/home/tarta/Lavoro/NIND_2/diffBind_E_nindensis/FHvsRH1/samplesheet_FHvsRH1_diffBind.csv"
output_name <- "FHvsRH1"

setwd(working_dir)

library(DiffBind)

# read in sampleSheet
diff <- dba(sampleSheet=sampleSheet)
diff

# count reads in peaks
diff.counted <- dba.count(diff, summits=0, fragmentSize=0) 
diff.counted

# set up contrast
diff.counted <- dba.contrast(diff.counted, categories=DBA_CONDITION, minMembers=2, reorderMeta = list(Condition = "FH"))
diff.counted

# perform Deseq2-based differential analysis
diff.analysed <- dba.analyze(diff.counted)
diff.analysed

# save your work
save.image(paste(output_name, ".RData", sep=""))

# make cross correlation matrix
pdf(paste("heatmap_", output_name, ".pdf", sep=""))
plot(diff.analysed, contrast=1)
dev.off()

# make MAplot
pdf(paste("MAplot_", output_name, ".pdf", sep=""))
dba.plotMA(diff.analysed)
dev.off()

# make volcano plot
pdf(paste("Vplot_", output_name, ".pdf", sep=""))
dba.plotVolcano(diff.analysed)
dev.off()


# make PCA plot
pdf(paste("PCAplot_", output_name, ".pdf", sep=""))
dba.plotPCA(diff.analysed, contrast=1)
dev.off()


# write csv file of differentially bound regions
report <- dba.report(diff.analysed)
write.csv(report, paste(output_name, "_report.csv", sep=""))

# save your work again
save.image(paste(output_name, ".RData", sep=""))


