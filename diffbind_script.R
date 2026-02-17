#!/usr/bin/env Rscript

# load DiffBind library
library(DiffBind)

# set working environment and input dataset
working_dir <- "/home/jacopo/MOA_teff_paper/cross9/"
sampleSheet <- "/home/jacopo/MOA_teff_paper/cross9/samplesheet_cross9.csv"
output_name <- "cross9"
setwd(working_dir)


#### read in samplesheet ####
diff <- dba(sampleSheet=sampleSheet)
diff

#### consensus peaks sets ####

# genomic intervals that are identified as peaks in one or more of the samples
#olap.rate <- dba.overlap(diff, mode=DBA_OLAP_RATE)
#olap.rate

# make a plot of how many peaks overlap in how many samples 
#plot(olap.rate, xlab="Overlapping samples", ylab="Overlapping peaks", type="b")

# default consensus peak set includes merged peaks that overlap in at least 2 samples
# use "peaks=consensus.peaks" in dba.count() if you want to use this peakset
#consensus.peaks <- dba.peakset(diff, bRetrieve=TRUE)
#consensus.peaks[,0]

# count reads in peaks
# use summits=0 and fragmentSize=0 for MOA expriment
diff.counted <- dba.count(diff, summits=0, fragmentSize=0)
diff.counted


#### Examining and normalizing the binding matrix ####

# make a correlation heatmap
pdf(paste("Heatmap_notNorm_", output_name, ".pdf", sep=""))
plot(diff.counted)
dev.off()


# make a PCA plot
pdf(paste("PCA_notNorm_", output_name, ".pdf", sep=""))
dba.plotPCA(diff.counted, DBA_REPLICATE, label=DBA_CONDITION)
dev.off()


# Normalization (default)

# make MA plot to show how the fold changes between conditions are distributed as the mean concentration of reads counts increases
pdf(paste("MAplot_notNorm_", output_name, ".pdf", sep=""))
dba.plotMA(diff.counted, bNormalized=FALSE, sub="Non-Normalized", contrast=list(WW=diff.counted$masks$WW, DS=diff.counted$masks$DS))
dev.off()


# make normalization (default)
diff.norm <- dba.normalize(diff.counted)

pdf(paste("MAplot_norm_", output_name, ".pdf", sep=""))
dba.plotMA(diff.norm, sub="Normalized (Default)", contrast=list(WW=diff.norm$masks$WW, DS=diff.norm$masks$DS))
dev.off()


#### Modeling and testing ####

# set up the contrast
diff.contrast <- dba.contrast(diff.norm, categories=DBA_CONDITION, minMembers=2, reorderMeta=list(Condition="WW"))
diff.contrast

# Fitting and testing
diff.model <- dba.analyze(diff.contrast)
dba.show(diff.model,bContrasts=TRUE)


#### Examining analysis results ####
diff.db <- dba.report(diff.model)
diff.db

# site where log2foldC is greater than 0
sum(diff.db$Fold>0)

# site where log2foldC is smaller than 0
sum(diff.db$Fold<0)

# make MAplot
pdf(paste("MAplot_norm_", output_name, ".pdf", sep=""))
dba.plotMA(diff.model)
dev.off()

# make volcano plot
pdf(paste("Vplot_norm_", output_name, ".pdf", sep=""))
dba.plotVolcano(diff.model)
dev.off()

# make cross correlation matrix
pdf(paste("heatmap_norm_", output_name, ".pdf", sep=""))
plot(diff.model, contrast=1)
dev.off()

# make PCA plot
pdf(paste("PCAplot_norm_", output_name, ".pdf", sep=""))
dba.plotPCA(diff.model, contrast=1)
dev.off()


# write csv file of differentially bound regions
report <- dba.report(diff.model)
write.csv(report, paste(output_name, "_report.csv", sep=""))



#### Normalization comparison #### 

#diff.norm = dba.normalize(diff.counted, method=DBA_ALL_METHODS, normalize="lib", library="full", background=FALSE, offsets=FALSE)
#diff.norm = dba.normalize(diff.counted, method=DBA_ALL_METHODS, normalize="lib", library="RiP", background=FALSE, offsets=FALSE)
#diff.norm = dba.normalize(diff.counted, method=DBA_ALL_METHODS, normalize="RLE", library="RiP", background=FALSE, offsets=FALSE)
#diff.norm = dba.normalize(diff.counted, method=DBA_ALL_METHODS, normalize="TMM", library="RiP", background=FALSE, offsets=FALSE)
#diff.norm = dba.normalize(diff.counted, method=DBA_ALL_METHODS, normalize="RLE", library="background", background=TRUE, offsets=FALSE)
#diff.norm = dba.normalize(diff.counted, method=DBA_ALL_METHODS, normalize="TMM", library="background", background=TRUE, offsets=FALSE)

#diff.model <- dba.analyze(diff.norm)
#dba.plotMA(diff.model, bNormalized=TRUE)
#dba.plotMA(diff.norm, sub="Normalized (Default)", contrast=list(WW=diff.norm$masks$WW, DS=diff.norm$masks$DS))

