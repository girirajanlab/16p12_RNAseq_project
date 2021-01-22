#R pipeline for reading isoform quantification data (Rsem) into R, and performing pairwise comparison of isoforms between samples using DeSeq2.
#txImport pipeline modified from: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

library(tximport)
library(DESeq2)

dir <- "./output" #Set working directory
samples <- read.table("samples.txt") #Read in text file containing sample IDs
files <- file.path(dir, paste0(samples$V1, ".rsem.isoforms.results")) #Generate list of RSEM file locations

txi.rsem <- tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE) #Create txImport object for isoform-level quantification

#Create sample table with biological replicates in file order
sampleTable <- data.frame(condition = factor(c("SG001","SG003","SG003","SG006","SG007","SG007","SG002","SG011","SG011","SG025","SG022","SG022","SG024",
"SG021","SG021","SG023","SG026","SG026","SG155","SG030","SG030","SG027","SG031","SG069","SG001","SG001","SG003","SG006","SG006","SG007","SG002","SG011",
"SG025","SG025","SG022","SG024","SG024","SG021","SG023","SG026","SG155","SG155","SG030","SG027","SG027","SG031","SG069","SG069","SG040","SG040","SG039",
"SG038","SG038","SG037","SG045","SG045","SG046","SG041","SG041","SG042","SG043","SG043","SG044","SG151","SG151","SG149","SG148","SG152","SG150","SG152",
"SG148","SG002","SG040","SG039","SG039","SG038","SG037","SG037","SG045","SG046","SG046","SG041","SG042","SG042","SG043","SG044","SG044","SG151","SG149",
"SG149","SG152","SG150","SG150","SG148","SG031","SG023"))) 

rownames(sampleTable) <- colnames(txi.rsem$counts) #Align sample names in DeSeq2 table to txImport table
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition) #Import txImport object into DeSeq2

#Begin DeSeq2 pipeline
#Pre-filter to remove rows with <2 total counts; additional filtering will be applied later
dds<-dds[rowSums(counts(dds))>1,]

#Apply differential expression analysis binomial model to DeSeq2 object
dds<-DESeq(dds)

#Generate pairwise comparisions for each parent-child pair. P-value is set to 0.05, default in DESeq is 0.1

#PSU001
SG001_SG006<-results(dds, contrast=c("condition", "SG001", "SG006"), alpha=0.05)
SG001_SG007<-results(dds, contrast=c("condition", "SG001", "SG007"), alpha=0.05)
SG011_SG006<-results(dds, contrast=c("condition", "SG011", "SG006"), alpha=0.05)
SG011_SG007<-results(dds, contrast=c("condition", "SG011", "SG007"), alpha=0.05)
SG006_SG002<-results(dds, contrast=c("condition", "SG006", "SG002"), alpha=0.05)
SG006_SG003<-results(dds, contrast=c("condition", "SG006", "SG003"), alpha=0.05)

#PSU004
SG021_SG024<-results(dds, contrast=c("condition", "SG021", "SG024"), alpha=0.05)
SG021_SG025<-results(dds, contrast=c("condition", "SG021", "SG025"), alpha=0.05)
SG022_SG024<-results(dds, contrast=c("condition", "SG022", "SG024"), alpha=0.05)
SG022_SG025<-results(dds, contrast=c("condition", "SG022", "SG025"), alpha=0.05)
SG023_SG024<-results(dds, contrast=c("condition", "SG023", "SG024"), alpha=0.05)
SG023_SG025<-results(dds, contrast=c("condition", "SG023", "SG025"), alpha=0.05)

#PSU005
SG026_SG027<-results(dds, contrast=c("condition", "SG026", "SG027"), alpha=0.05)
SG026_SG069<-results(dds, contrast=c("condition", "SG026", "SG069"), alpha=0.05)
SG155_SG027<-results(dds, contrast=c("condition", "SG155", "SG027"), alpha=0.05)
SG155_SG069<-results(dds, contrast=c("condition", "SG155", "SG069"), alpha=0.05)
SG069_SG030<-results(dds, contrast=c("condition", "SG069", "SG030"), alpha=0.05)
SG069_SG031<-results(dds, contrast=c("condition", "SG069", "SG031"), alpha=0.05)

#PSU007
SG037_SG039<-results(dds, contrast=c("condition", "SG037", "SG039"), alpha=0.05)
SG037_SG040<-results(dds, contrast=c("condition", "SG037", "SG040"), alpha=0.05)
SG038_SG039<-results(dds, contrast=c("condition", "SG038", "SG039"), alpha=0.05)
SG038_SG040<-results(dds, contrast=c("condition", "SG038", "SG040"), alpha=0.05)
SG041_SG043<-results(dds, contrast=c("condition", "SG041", "SG043"), alpha=0.05)
SG041_SG044<-results(dds, contrast=c("condition", "SG041", "SG044"), alpha=0.05)
SG042_SG043<-results(dds, contrast=c("condition", "SG042", "SG043"), alpha=0.05)
SG042_SG044<-results(dds, contrast=c("condition", "SG042", "SG044"), alpha=0.05)
SG040_SG045<-results(dds, contrast=c("condition", "SG040", "SG045"), alpha=0.05)
SG040_SG046<-results(dds, contrast=c("condition", "SG040", "SG046"), alpha=0.05)
SG044_SG045<-results(dds, contrast=c("condition", "SG044", "SG045"), alpha=0.05)
SG044_SG046<-results(dds, contrast=c("condition", "SG044", "SG046"), alpha=0.05)

#PSU052
SG148_SG151<-results(dds, contrast=c("condition", "SG148", "SG151"), alpha=0.05)
SG148_SG152<-results(dds, contrast=c("condition", "SG148", "SG152"), alpha=0.05)
SG149_SG151<-results(dds, contrast=c("condition", "SG149", "SG151"), alpha=0.05)
SG149_SG152<-results(dds, contrast=c("condition", "SG149", "SG152"), alpha=0.05)
SG150_SG151<-results(dds, contrast=c("condition", "SG150", "SG151"), alpha=0.05)
SG150_SG152<-results(dds, contrast=c("condition", "SG150", "SG152"), alpha=0.05)

#Export results to text files
#PSU001
write.csv(as.data.frame(SG001_SG006), file="SG001_SG006_isoform_deseq2.csv")
write.csv(as.data.frame(SG001_SG007), file="SG001_SG007_isoform_deseq2.csv")
write.csv(as.data.frame(SG011_SG006), file="SG011_SG006_isoform_deseq2.csv")
write.csv(as.data.frame(SG011_SG007), file="SG011_SG007_isoform_deseq2.csv")
write.csv(as.data.frame(SG006_SG002), file="SG006_SG002_isoform_deseq2.csv")
write.csv(as.data.frame(SG006_SG003), file="SG006_SG003_isoform_deseq2.csv")

#PSU004
write.csv(as.data.frame(SG021_SG024), file="SG021_SG024_isoform_deseq2.csv")
write.csv(as.data.frame(SG021_SG025), file="SG021_SG025_isoform_deseq2.csv")
write.csv(as.data.frame(SG022_SG024), file="SG022_SG024_isoform_deseq2.csv")
write.csv(as.data.frame(SG022_SG025), file="SG022_SG025_isoform_deseq2.csv")
write.csv(as.data.frame(SG023_SG024), file="SG023_SG024_isoform_deseq2.csv")
write.csv(as.data.frame(SG023_SG025), file="SG023_SG025_isoform_deseq2.csv")

#PSU005
write.csv(as.data.frame(SG026_SG027), file="SG026_SG027_isoform_deseq2.csv")
write.csv(as.data.frame(SG026_SG069), file="SG026_SG069_isoform_deseq2.csv")
write.csv(as.data.frame(SG155_SG027), file="SG155_SG027_isoform_deseq2.csv")
write.csv(as.data.frame(SG155_SG069), file="SG155_SG069_isoform_deseq2.csv")
write.csv(as.data.frame(SG069_SG030), file="SG069_SG030_isoform_deseq2.csv")
write.csv(as.data.frame(SG069_SG031), file="SG069_SG031_isoform_deseq2.csv")

#PSU007
write.csv(as.data.frame(SG037_SG039), file="SG037_SG039_isoform_deseq2.csv")
write.csv(as.data.frame(SG037_SG040), file="SG037_SG040_isoform_deseq2.csv")
write.csv(as.data.frame(SG038_SG039), file="SG038_SG039_isoform_deseq2.csv")
write.csv(as.data.frame(SG038_SG040), file="SG038_SG040_isoform_deseq2.csv")
write.csv(as.data.frame(SG041_SG043), file="SG041_SG043_isoform_deseq2.csv")
write.csv(as.data.frame(SG041_SG044), file="SG041_SG044_isoform_deseq2.csv")
write.csv(as.data.frame(SG042_SG043), file="SG042_SG043_isoform_deseq2.csv")
write.csv(as.data.frame(SG042_SG044), file="SG042_SG044_isoform_deseq2.csv")
write.csv(as.data.frame(SG040_SG045), file="SG040_SG045_isoform_deseq2.csv")
write.csv(as.data.frame(SG040_SG046), file="SG040_SG046_isoform_deseq2.csv")
write.csv(as.data.frame(SG044_SG045), file="SG044_SG045_isoform_deseq2.csv")
write.csv(as.data.frame(SG044_SG046), file="SG044_SG046_isoform_deseq2.csv")

#PSU052
write.csv(as.data.frame(SG148_SG151), file="SG148_SG151_isoform_deseq2.csv")
write.csv(as.data.frame(SG148_SG152), file="SG148_SG152_isoform_deseq2.csv")
write.csv(as.data.frame(SG149_SG151), file="SG149_SG151_isoform_deseq2.csv")
write.csv(as.data.frame(SG149_SG152), file="SG149_SG152_isoform_deseq2.csv")
write.csv(as.data.frame(SG150_SG151), file="SG150_SG151_isoform_deseq2.csv")
write.csv(as.data.frame(SG150_SG152), file="SG150_SG152_isoform_deseq2.csv")
