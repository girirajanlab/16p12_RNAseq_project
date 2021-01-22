#R pipeline for reading gene quantification data (Rsem) into R, and performing weighted gene correlation network analysis using DeSeq2 and WGCNA
#txImport pipeline modified from: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

library("tximport")
library("DESeq2")
library("WGCNA")
library("sva") #Includes ComBat function

#Set environemnt options for WGCNA library
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

dir <- "./output" #Set working directory
samples <- read.table("samples.txt") #Read in text file containing sample IDs
files <- file.path(dir, paste0(samples$V1, ".rsem.genes.results")) #Generate list of RSEM file locations

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE) #Create txImport object for gene-level quantification

#Create sample table with biological replicates in file order
sampleTable <- data.frame(condition = factor(c("SG001","SG003","SG003","SG006","SG007","SG007","SG002","SG011","SG011","SG025","SG022","SG022","SG024",
"SG021","SG021","SG023","SG026","SG026","SG155","SG030","SG030","SG027","SG031","SG069","SG001","SG001","SG003","SG006","SG006","SG007","SG002","SG011",
"SG025","SG025","SG022","SG024","SG024","SG021","SG023","SG026","SG155","SG155","SG030","SG027","SG027","SG031","SG069","SG069","SG040","SG040","SG039",
"SG038","SG038","SG037","SG045","SG045","SG046","SG041","SG041","SG042","SG043","SG043","SG044","SG151","SG151","SG149","SG148","SG152","SG150","SG152",
"SG148","SG002","SG040","SG039","SG039","SG038","SG037","SG037","SG045","SG046","SG046","SG041","SG042","SG042","SG043","SG044","SG044","SG151","SG149",
"SG149","SG152","SG150","SG150","SG148","SG031","SG023"))) 

rownames(sampleTable) <- colnames(txi.rsem$counts) #Align sample names in DeSeq2 table to txImport table
txi.rsem$length[txi.rsem$length == 0] <- 1 #Fix bug where some gene lengths are set to 0, as per https://support.bioconductor.org/p/84304/
dds_wgcna <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition) #Import txImport object into DeSeq2

#Begin DeSeq2 pipeline: Transform/filter raw counts of biological replicate RNA-Seq data
#Filter to remove rows with <30 total counts (10 counts/replicate in >1 sample)
dds_wgcna<-dds_wgcna[rowSums(counts(dds_wgcna))>30,]

#Estimate size and dispersion factors for data
dds_wgcna<-estimateSizeFactors(dds_wgcna)
dds_wgcna<-estimateDispersions(dds_wgcna)

#Perform variance-stabilizing transformation on raw counts to produce recommended transformed/normalized counts
vsd <- getVarianceStabilizedData(dds_wgcna)

#Export transformed data
write.csv(as.data.frame(vsd), file="16p12_LCL_deseq2_biol_reps_vsd.csv")

#Export normalized counts to files for visualization
write.csv(counts(dds_wgcna, normalized=TRUE), file="norm_counts_deseq2.csv")

#Input data into WGCNA and perform pre-processing

#Load and transform data
wgcna_data<-read.csv("16p12_LCL_deseq2_biol_reps_vsd.csv",header=TRUE,row.names=1)
wgcna_data<-as.data.frame(t(wgcna_data))

#Filter genes for missing values; will likely be OK due to filtering in step 1
gsg<-goodSamplesGenes(wgcna_data,verbose=3)
gsg$allOK #Should be equal to TRUE
wgcna_data<-wgcna_data[gsg$goodSamples, gsg$goodGenes] #Only use this to remove genes with missing values


#Cluster samples to see if there are any outliers
sampleTree<-hclust(dist(wgcna_data), method="average")
pdf("hierarchal_clustering.pdf")
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5,cex.axis=1.5, cex.main=2)
dev.off()

#Use ComBat to account for family-specific batch effects
#Code and PlotPCA function from: https://peterlangfelder.com/2018/12/02/removal-of-unwanted-variation-based-on-a-subset-of-samples/

annotation2=data.frame(status=c("PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU004","PSU004","PSU004","PSU004",
"PSU004","PSU004","PSU004","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001",
"PSU004","PSU004","PSU004","PSU004","PSU004","PSU004","PSU004","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU007","PSU007","PSU007",
"PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU052","PSU052","PSU052","PSU052","PSU052","PSU052","PSU052",
"PSU052","PSU001","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU052","PSU052",
"PSU052","PSU052","PSU052","PSU052","PSU052","PSU005","PSU004"),batch=c("PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU004","PSU004","PSU004","PSU004",
"PSU004","PSU004","PSU004","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001",
"PSU004","PSU004","PSU004","PSU004","PSU004","PSU004","PSU004","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU007","PSU007","PSU007",
"PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU052","PSU052","PSU052","PSU052","PSU052","PSU052","PSU052",
"PSU052","PSU001","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU052","PSU052",
"PSU052","PSU052","PSU052","PSU052","PSU052","PSU005","PSU004"))

annotation=data.frame(status=c("1","1","1","1","0","0","0","1","1","1","1","1","0","1","1","0","1","1","1","1","1","0","0","1","1","1","1","1","1","0","0","1",
"1","1","1","0","0","1","0","1","1","1","1","0","0","0","1","1","1","1","0","0","0","1","1","1","0","1","1","1","0","0","1","1","1","0","1","0","0","0",
"1","0","1","0","0","0","1","1","1","0","0","1","1","1","0","1","1","1","0","0","0","0","0","1","0","0"),batch=c("PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU004","PSU004","PSU004","PSU004",
"PSU004","PSU004","PSU004","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001","PSU001",
"PSU004","PSU004","PSU004","PSU004","PSU004","PSU004","PSU004","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU005","PSU007","PSU007","PSU007",
"PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU052","PSU052","PSU052","PSU052","PSU052","PSU052","PSU052",
"PSU052","PSU001","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU007","PSU052","PSU052",
"PSU052","PSU052","PSU052","PSU052","PSU052","PSU005","PSU004"))

#Generate PCA plots for clustering of samples by carrier status and family
pdf("pca_by_carrier_status.pdf")
plotPCA(data, annotation, main = "Merged data before batch correction")
dev.off()
pdf("pca_by_family.pdf")
plotPCA(data, annotation2, main = "Merged data before batch correction")
dev.off()

#Apply ComBat using family as batch effect, including disease status (carrier vs non-carrier) as covariant
data.ComBat = t(ComBat(t(data), batch = annotation$batch, mod = model.matrix(~status, data = annotation)))

#Generate PCA plost post-correction
pdf("Combat_corr_pca_by_carrier_status.pdf")
plotPCA(data.ComBat, annotation, main = "Corrected application of ComBat")
dev.off()
pdf("Combat_corr_pca_by_family.pdf")
plotPCA(data.ComBat, annotation2, main = "Corrected application of ComBat")
dev.off()

#Repeat hierarchal clustering post-ComBat
sampleTree<-hclust(dist(data.ComBat), method="average")
pdf("hierarchal_clustering_combat.pdf")
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5,cex.axis=1.5, cex.main=2)
dev.off()


#Automatic network construction and module detection

#Generate a set of soft-threshold powers (parameters) to try with the data
powers<-c(c(1:10), seq(from = 12, to=20, by=2))
#Use network topology analysis function to select best data parameter
sft<-pickSoftThreshold(data.ComBat, powerVector = powers, verbose = 5)

#Plot network topology results (scale independence and mean connectivity). Save result as PDF.
pdf("scale_independence_combat.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],type="n",xlab="Soft Threshold (power)",
	ylab="Scale Free Topology Model Fit,signed R^2",main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.9,col="red")
abline(h=0.90,col="red")
dev.off()

pdf("mean_connectivity_combat.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5],type="n",xlab="Soft Threshold (power)",ylab="Mean Connectivity",main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()
#Choose the power that, on the scale independence curve (1st plot), is the lowest at the "plateau" of the curve
#For this dataset, 8 was selected as the power threshold value

#Run WGCNA using the selected power above
net<-blockwiseModules(data.ComBat,power=8,networkType="signed hybrid",TOMType="unsigned",minModuleSize = 30,maxBlockSize=30000,reassignThreshold=0,
	mergeCutHeight=0.25,numericLabels=TRUE,pamRespectsDendro=FALSE,saveTOMs=FALSE,verbose=3,corType="bicor",maxPOutliers=0.05,method="hybrid",deepSplit=2)

#View number of modules and number of genes/module
table(net$colors)
#Plot dendrogram of gene clustering. Save as PDF
mergedColors<-labels2colors(net$colors)
pdf("wgcna_clusters_dendrogram_combat.pdf")
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()

#Step 4: Visualize network data

#Make heatmap and clustering of all modules/"eigengenes"
#See tutorial for how to separate these plots
MEs<-moduleEigengenes(wgcna_data, mergedColors)$eigengenes
pdf("eigengene_networks_combat.pdf")
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()

#Make heatmap of all genes using above power variable; save as PNG
options(device=png)
dissTOM<-1-TOMsimilarityFromExpr(wgcna_data,power=8)
plotTOM<-dissTOM^6
diag(plotTOM)=NA
moduleColors<-labels2colors(net$colors)
geneTree<-net$dendrograms[[1]]
sizeGrWindow(9,9)
png("wgcna_heatmap_combat.png")
TOMplot(plotTOM,geneTree,moduleColors,main ="Network heatmap plot, all genes")
dev.off()

#Step 5: Export modules as lists of genes

#Get list of all modules
moduleColorNames<-names(table(mergedColors))

for(module in moduleColorNames)
{
	modGenes<-(mergedColors==module)
	modGeneNames<-names(wgcna_data[modGenes])
	fileName<-paste("16p12_LCL_wgcna_",module,"_module_combat.txt",sep="")
	write.table(as.data.frame(modGeneNames),file=fileName,row.names=FALSE,col.names=FALSE)
}