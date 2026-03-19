########################################################################
### DESeq2 RNA Seq Analysis Script for multiple factors / treatments ###
########################################################################

# Author: Hannes Reinwald
# Maintainer: Fabian Essfeld
# Contact: fabian.essfeld@ime.fraunhofer.de

# README ### -----------------------------------------------------------
# This script is desgined to run DESeq2 normalization and statisticalsourc testing on RNAseq experiments
# in Ecotox testing for multiple concentrations. This script will automatically adopt to the k numbers
# of samples and N numbers of conditions (High,Mid,Low or Treat1, Treat2, Treat3, ... , Control) and 
# will generate Log2FC values with respect to control groups.

# Requiered Input: 
# - CountMatrix (txt or csv is fine, if csv read.csv2 is used)
# - coldata.csv

# Main Analysis Steps are: 
# 1) RLE normalization (DESeq2) across all samples element of CountMatrix
# 2) Calc LFC (with respect to control) and LFcutOff (upper 90% quantile of abs(LFC values))
# 3) apeglm shrinking on LFC

# 4) Multiple t-testing with Benjamin-Hochberg correction (padj < 0.05) 
#    and independet hypothesis weighing (IHW) to identify DEGs
#    for LFC and apeglm(LFC) values for H0: LFC = 0; apeglm(LFC) = 0
# 4.1) Annotation of Genelists via org.Dr.eg.db [https://www.bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html]
# 4.2) DEGs with LFC / apeglm(LFC) > < abs(LFcutOff) are kept as potential molecular marker features

# 5) Plotting
# 5.1 MA- & Vulcano plots for N-1 conditions for LFC [LFcut & padj cut]
# 5.2 MA- & Vulcano plots for N-1 conditions for apeglm(LFC) [LFcut & padj cut]
# 5.3 DataQC
#     - Correlation bio. replicates (rlog(meanCounts))
#     - MeanSD plots (rlog(meanCounts))
#     - RLE_normalisation
#     - normCounts_rlogTransformation

# 6) Output: 
# 6.1 DESeq2 Results (xcel; still working on html output though ... )
#     - N-1 DESeq2 result files for LFC
#     - N-1 DESeq2 result files for apeglm(LFC)
# 6.2 DEGs (padj < 0.05 and LFcut for: 
#     - N-1 files for LFC
#     - N-1 files for apeglm(LFC)


# LOAD PACKAGES --------------------------------------------------------------
require(DESeq2) 
require(IHW)    
require(apeglm) 
require(qvalue)
#require(locfdr)
#require(RColorBrewer)
require(ggplot2)
require(ggpubr)
require(hexbin)


# LOAD FUNCTIONS --------------------------------------------------------------

# Funtion files are loaded from S: or mint. If working locally, make sure to
# copy this directory in your working directory.

# Define the potential directories
directories <- c(
  "fun"
)

# Initialize an empty list for files
files <- character()

# Loop through directories and gather files
for (dir in directories) {
  if (dir.exists(dir)) {
    files <- list.files(dir, full.names = TRUE)
    if (length(files) > 0) {
      break  # Exit loop if files are found
    }
  }
}

# Check if files were found
if (length(files) == 0) {
  stop("No files found in any of the specified directories.")
}

# Load the files using source()
for (file in files) {
  source(file)
}

rm(files, file, directories)

# DATA IMPORT & Color settings ---------------------------------------------------
message("\nImporting CountMatrix and coldata files ...")

## coldata import ##
tmp = list.files(full.names = T, pattern = "[Cc]oldata[.]")
if(grepl(".csv",tmp)==T){
  coldata = read.csv2(file = tmp, row.names = 1, header = T)
  # In case csv file was not differently exported run read.csv
  if (ncol(coldata) <= 1) {
    coldata = read.csv(file = tmp, row.names = 1, header = T)
  }
} else {
  # if not ending with csv import with read.delim2
  coldata = read.delim2(file = tmp, row.names = 1, header = T)
  if (ncol(coldata) <= 1) {
    # if file was differently formated use read.delim
    coldata = read.delim(file = tmp, row.names = 1, header = T)
  }
}
coldata <- droplevels(coldata[which(coldata$Condition != ""), # filtering for non empty rows
                              which(!(grepl("X.",colnames(coldata))) == T)]) # filtering cols
coldata$Condition <- gsub(" ","",coldata$Condition) #removing any spaces in Condition levels
coldata$Tank <- gsub("-","_",coldata$Tank)

# set factor levels
coldata$Tank      <- factor(coldata$Tank)
coldata$Condition <- factor(coldata$Condition)
coldata$Substance <- factor(coldata$Substance)
substance         <- levels(coldata$Substance) #Extract the name(s) of the tested substance(s)

# relevel factors in Condition in right order:
condition = levels(coldata$Condition)
tmp = as.character(coldata$Condition[grep("[Cc]ontrol",coldata$Condition)][1]) #extract control factor level
if(all(grepl("[Ee]xposure",condition[2:4]))==T & length(condition)==4) {# Control, LE, ME, HE
  coldata$Condition <- factor(coldata$Condition, levels = condition[c(1,3,4,2)])
} else if (all(grepl("[Ee]xposure",condition[2:3]))==T & length(condition)==3) { # Control, LE, HE
  coldata$Condition <- factor(coldata$Condition, levels = condition[c(1,3,2)])
} else { # just specifying the reference level --> Control
  coldata$Condition <- relevel(coldata$Condition, ref = tmp)
}
message(paste0("Reference level:\t",levels(coldata$Condition)[1]),"\nSorting ",length(condition)-1," Treatments:\t",paste(levels(coldata$Condition)[-1], collapse = " "))
condition <- levels(coldata$Condition)

## CountMatrix import ##
# Import all provided count mtx files and merge them together;
# then filter according to the provided sample ID in coldata

cmtx.ls = list()
tmp = list.files(full.names = F, pattern = "[Cc]ount[Mm]atrix")
for(i in tmp){ cmtx.ls[[sub("[Cc]ount[Mm]atrix","",i)]] <- importCountMtx(i) }
if(length(cmtx.ls) == 1){
  count.matrix = cmtx.ls[[1]]
} else {
  count.matrix <- countMtxMerge(cmtx.ls)
}

# subset cols in count.matrix after rows in coldata. So there is only the need to trim the coldata file manually
count.matrix <- count.matrix[,rownames(coldata)] 
if(all(sort(row.names(coldata)) == sort(colnames(count.matrix))) == F) {
  stop("\nCount matrix colnames and coldata rownames do not match up!\nPlease check your input files!")
}
stopifnot(all(rownames(coldata) == colnames(count.matrix))) # Must be TRUE!
message("Done!\n")
rm(i,tmp, cmtx.ls)

### Set annotation colors for downstream plotting ###
col.Cond = colorRampPalette(RColorBrewer::brewer.pal(n=8, name="YlOrRd"))(length(levels(coldata$Condition)))
col.Cond[1] <- "#56B4E9" #Defines color for control
col.Tank = colorRampPalette(c("gray95","gray50"))(length(levels(coldata$Tank)))
ann_colors = list(Tank = setNames(col.Tank, levels(coldata$Tank)),
                  Condition = setNames(col.Cond, levels(coldata$Condition)))
### Setup output environment ###
dir.create("DESeq2_Pairwise", showWarnings = F)
setwd("DESeq2_Pairwise")
home = getwd()


#  DESeq2   -------------------------------------------------------------------

## Count Matrix filtering -----------------------------------------------------
message(paste0("\nRemoving low abundant counts from count matrix.",
               "\nMinimum number of gene counts per row: \t\t",ncol(count.matrix)))

## Create condition groups ID list
id.ls = list()
for(i in condition) {id.ls[[i]] = rownames(coldata[coldata$Condition %in% i,])}
stopifnot(length(id.ls) > 1) # Stop script here if there is only one condition provided
# We will need this object later ... 

## NEW - Threshold criteria ----------------------------------------
# We will use the CPM filtering threshold ("relevance threshold") suggested by Verheijen et al 2022.
# excludes all genes that do not have at least 75% (65% in case of only 3 relicates)of their replicates expressed at 1 CPM in
# atleast one of the experimental conditions. This will remove low abundant gene counts prior DESeq2
# DGEA which can speed up things

#CPM.matrix = my_CPM(count.matrix)
#boxplot(log10(CPM.matrix+1))
#barplot(apply(CPM.matrix,2,sum))


## Filter input matrix 
minCPM <- 1 # Minimum counts per million in some samples

count.matrix.cl = relevance.filter(count.matrix, id.ls, minCPM)


## Create DESeq2 object -------------------------------------------------------
p  = 0.05 # padj cutoff
Qt = 0.75  # %-Quantile for effect size cut off (i.e. 0.9 => 90% = top 10% abs(LFC))
message(paste0("\n Starting DESeq2 Analysis - Pairwise Wald's t-test and IHW \n padj < ",p,"\n Tested Substance: \t",substance,"\n"))

# Check if we have a balanced test design.
t1 = all(table(coldata$Condition) == table(coldata$Condition)[1])
t2 = all(table(coldata$Tank) == table(coldata$Tank)[1])
t3 = length(unique(coldata$Tank)) #to check if there are more than one tank condition.
if(all(c(t1,t2) == T) & t3 > 1){
  # if TRUE we have a balanced test design between Conditions and Tanks.
  # Hence we can apply a multifactor model as: ~ Tank + Condition
  message("Test design of experiment is BALANCED.\nApplying multi-factor model:\t~ Tank + Condition")
  dds <- DESeqDataSetFromMatrix(countData = count.matrix.cl, colData = coldata,
                                design = ~ Tank + Condition) 
  # Multifactor model controlling for variability in Tanks
} else {
  # if FALSE test design is not balanced and we can not directly apply the batch effect correction here.
  # Hence the model design in this case is only ~ Tank
  message("Test design of experiment is NOT BALANCED.\nApplying uni-factor model:\t~ Condition")
  dds <- DESeqDataSetFromMatrix(countData = count.matrix.cl, colData = coldata,
                                design = ~ Condition)
}
rm(t1,t2,t3)


# Get a list of different disp estimates for later QC plotting
fitType = c("parametric", "local") #additonally "mean" possible
estDisp.ls <- list()
for(i in fitType){
  message(paste("\nEstimate Dispersions for fitType: ",i))
  x <- estimateSizeFactors(dds)
  x <- estimateDispersions(x, fitType = i)
  estDisp.ls[[i]] <- x
}
rm(x,i,fitType)

## RUN DESeq() ----------------------------------------------------------------
message("\nRunning DESeq2 pairwise comparison ...")
dds <- DESeq(dds, test = "Wald")                   # Pairwise comparison
#dds <- DESeq(dds, test = "LRT", reduced = ~ Tank) # ANOVA-like approach

## Extract norm. & transformed count.matrix ------------------------------------
normMtx =list()
normMtx[["norm"]]   = round(counts(dds, normalized = T),3) # [[1]]  normalized read counts
normMtx[["ntd"]]    = round(assay(normTransform(dds)),3)   # [[2]] (n+1)log2 transformed norm counts
normMtx[["nrl"]]    = round(assay(rlog(dds, blind = F)),3) # [[3]] rlog transformed mean read counts
normMtx[["nrl.bl"]] = round(assay(rlog(dds, blind = T)),3) # [[4]] rlog transformed mean read counts; BLIND!

message("\nExporting DESeq2 normalized count matrix. This might take a while...\nSaving to txt table ...")
write.table(normMtx[[1]], paste0(substance,"_DESeqNormCounts.txt")) # Extract DESEq2 normalized count matrix
#write.table(df.rld, paste0(substance,"_rlogDESeqNormCounts.txt")) # Extract rlog(norm counts) matrix
message("Done!\n")

gc()

## Extract DESeq2 result tables to list objects --------------------------------
resNames = resultsNames(dds)[grepl("^[Cc]ondition_",resultsNames(dds))]
## Non-shrunk results
res.ls <- list() #create an empty list to store objects in
for(i in resNames) {
  x = gsub("[Cc]ondition_","", i)
  message(paste0("Extracting DESeq2 results for: ",x," [IHW, non-shrunk Log2FC] ..."))
  res.ls[[gsub("_.+","",x)]] <- results(dds, name = i, filterFun = ihw, alpha = p)
  message("Done!\n")
  gc()
}
## apeglm shunk results
resLfs.ls <- list() # save log2FC shrunk results in different list
for(i in resNames) {
  x = gsub("[cC]ondition_","", i)
  message(paste0("Extracting DESeq2 results for: ",x," [IHW, apeglm shrunk Log2FC] ..."))
  resLfs.ls[[gsub("_.+","",x)]] <- lfcShrink(dds, coef = i, type = "apeglm", 
                                             res = res.ls[[gsub("_.+","",x)]])
  message("Done!\n")
  gc()
}
rm(x,i)

## Append qvalue to res.ls / resLfs.ls -----------------------------------------
message("\nCalculating and appending qvalues to DESeq2 results ...")
# res.ls
res.ls = lapply(res.ls, function(x){
  x$qval = qvalue(x$pvalue)$qvalue
  x })
# resLfs.ls
resLfs.ls = lapply(resLfs.ls, function(x){
  x$qval = qvalue(x$pvalue)$qvalue
  x })
message("Done!\n")


##   biomaRt   ------------------------------------------------------------------

### Create annotation object for DESeq2 result tables --------------------------
message("\nAnnotating ENSEMBL Gene IDs in result tables using biomaRt ...")
# if rerio mart object found locally load it, create new mart object from host 
if(file.exists("S:/bioinformatics/biomaRt/drerio_mart.Robj")) {
  message("\nDanio rerio mart object found locally. \nLoading 'drerio_mart.Robj' into R session ...")
  load("S:/bioinformatics/biomaRt/drerio_mart.Robj")
} else {
    message("\nCould not find 'drerio_mart.Robj'. Creting new mart object from 'www.ensembl.org'",
            "\nMake sure you have a working interent connection. Otherwise this will fail!\nConnecting to server ...")
    rerio = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset="drerio_gene_ensembl",host="https://www.ensembl.org")
}

### Build gene annotation object & merge ----------------------------------------
id     = rownames(res.ls[[1]])
idType = "ensembl_gene_id"
attr   = c("ensembl_gene_id","external_gene_name","description","gene_biotype",
           "entrezgene_id")
message("\nAnnotating Ensembl gene IDs ...")
GeneAnno = biomaRt::getBM(attributes = attr, mart = rerio, uniqueRows = F,
                          filters = idType, values = id, useCache = F)
stopifnot(length(intersect(GeneAnno$ensembl_gene_id,id)) == length(id)) #Check Sum

# !!! In a few cases there can be more than one entrez id for a single ensembl gene id !!!
# Here we remove duplicated entries (~ 360 out of ~ 25 000, so it's really a minor fraction)
# As it turns out higher integer values in the entrez ID correspond to a deeper level of organization
# i.e. ID 999 -> Protein A, ID 1000348 -> Protein A subunit a1. For our purpose it is good enough 
# to retrieve the overall gene / protein information.
# => We sort GeneAnno by entrezgene_id and then remove the duplicates with !duplicate().
# That way, in case multiple entrez ids are assigned to a single ensembl id the lower entrez id
# will be kept for downstream analysis. 
GeneAnno = dplyr::arrange(GeneAnno, entrezgene_id)
GeneAnno = GeneAnno[!duplicated(GeneAnno$ensembl_gene_id),] #keep only non duplicated
stopifnot(length(id) == nrow(GeneAnno)) #Check Sum
row.names(GeneAnno) = GeneAnno$ensembl_gene_id
# To ensure compatability with downstream GSEA scripts rename external_gene_name & entrezgene_id
colnames(GeneAnno) = c("ensembl_gene_id","SYMBOL","description","biotype","ENTREZID")

# Now merge GeneAnno with DESeq2 result tables
AnnoFun = function(x){
  tmp = merge(as.data.frame(x), GeneAnno, by=0)[,-1]
  row.names(tmp) = tmp$ensembl_gene_id
  tmp[,-which(colnames(tmp) %in% "ensembl_gene_id")]
}
res.ls    = lapply(res.ls, FUN = AnnoFun)
resLfs.ls = lapply(resLfs.ls, FUN = AnnoFun)
# Finally sort res.ls / resLfs.ls after the padj value
res.ls    = lapply(res.ls, function(x){dplyr::arrange(x, padj)})
resLfs.ls = lapply(resLfs.ls, function(x){dplyr::arrange(x, padj)})
rm(id, rerio)
message("Done!\n")



##  EXPORT DESeq2 RESULTS  -------------------------------------------------------
message("Exporting DESeq2 result tables. This might take a while :)")
dir.create("Results", showWarnings = F)
for(k in names(res.ls)){
  n = gsub("[Ee]xposure","",k)
  message(paste("Saving DESeq2 result table",k,"[IHW, non-shrunk Log2FC] to csv ..."))
  write.csv2(res.ls[[k]], file = paste0("Results/",substance,"_res_",n,".csv"))
  message(paste("Saving DESeq2 result table",k,"[IHW, apeglm shrunk Log2FC] to csv ..."))
  write.csv2(resLfs.ls[[k]], file = paste0("Results/",substance,"_reslfs_",n,".csv"))
}
message("Done!\n")

## Select only DEGs for export and downstream plotting ##
deg.ls    <- lapply(res.ls, function(x){subset(x, padj <= p)}) # pcut
degCut.ls <- lapply(res.ls, function(X){                       # pcut & log2FC cutoff 
  lfcut = quantile(abs(X$log2FoldChange), Qt) #(90% quantile of abs(log2FC))
  if(lfcut > 1){lfcut = 1} # in case lfcut greater 1, set to 1
  x = subset(X, padj <= p)
  subset(x, abs(log2FoldChange) >= lfcut)
})
# Selecting only DEGs with padj < p & abs(apeglm shrunk Log2FC) > lfcut (90% quantile from abs(non shrunk log2FC))
degLfs.ls = list()
for(k in names(resLfs.ls)){
  # determine lfc cut off based on non-shrunk values
  lfcut = quantile(abs(res.ls[[k]]$log2FoldChange), Qt)
  if(lfcut > 1){lfcut = 1} # in case lfcut greater 1, set to 1
  message(paste("Log2FC cut off for:\t",k,"\t----->\t",round(lfcut,2),"(Top",(1-Qt)*100,"% Quantile)"))
  # Select padj < p from resLfs.ls and subset by log2FC
  degLfs.ls[[k]] <- subset(subset(resLfs.ls[[k]], padj <= p), abs(log2FoldChange) >= lfcut)
}

## Exporting DEGs in separate result csv tables ##
dir.create("DEGs", showWarnings = F)
out = c("unshrink","unshrinkFCcut","lfsFCcut") # output dir. for different levels of filtering
for(k in out){dir.create(paste0("DEGs/",k), showWarnings = F)}
# Export to respective dir
for(k in names(deg.ls)) {
  message(paste("Saving selection of DEGs for",k,"to csv ..."))
  n = gsub("[Ee]xposure","",k)
  write.csv2(deg.ls[[k]],paste0("DEGs/",out[1],"/",substance,"_pcut_",n,".csv"))   # unshrink
  write.csv2(degCut.ls[[k]],paste0("DEGs/",out[2],"/",substance,"_pcut_LFcut_",n,".csv"))# unshrink + FCcut
  write.csv2(degLfs.ls[[k]],paste0("DEGs/",out[3],"/",substance,"_lfs_pcut_LFcut_",n,".csv"))# shrink + FCcut
}
message("Done!\n")
rm(n,k,lfcut,out)



#PLOTTING   ----------------------------------------------------------------
dir.create(paste0(home,"/Plots"), showWarnings = F)
dir.create(paste0(home,"/Plots/DataQC"), showWarnings = F)
message("\nStarting Data Composition and QC plotting ...")

## QC - DESeq2 norm. & count distr. & countMtx filtering  ### -----------
message("DESeq2 normalization & gene count distribution ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_CountNorm_mtxFiltering.pdf"),
    width = 9, height = 11.7)
par(mfrow=c(4,2))

## Normalization bar plots ##
barplot(colSums(counts(dds)),
        ylab= "Total gene counts",
        main= "Raw counts",
        las= 3, #rotating sample labels 90
        ylim= c(0,1.2*max(colSums(counts(dds)))))
barplot(colSums(counts(dds, normalized=T)),  #plot mean normalized counts
        ylab= "Total gene counts",
        main= "DESeq2 norm. counts",
        las= 3, #rotating sample labels 90
        ylim= c(0,1.2*max(colSums(counts(dds)))))
boxplot(log10(counts(dds)+1),
        ylab= "Log10(Gene counts + 1)",
        las = 3, #rotating sample labels 90
        main= "Raw counts")
boxplot(log10(counts(dds, normalized=T)+1),
        ylab= "Log10(Gene counts + 1)",
        las = 3, #rotating sample labels 90
        main= "DESeq2 norm. counts")

## Gene count distr. & countMtx filtering ##

geneCountDistribution(rawCounts = count.matrix,
                      filteredCounts = count.matrix.cl,
                      breaks = 50)

dev.off()


## QC - DESeq2's dispersion estimate models ### ---------------------------------------------
message("DESeq2's dispersion estimate models (estimateDispersions) ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_DispEstimates.pdf"),
    width = 7, height = 5.9, onefile = T, bg = "transparent")
for(i in names(estDisp.ls)){
  plotDispEsts(estDisp.ls[[i]], main = paste(substance,'- DESeq2 fit type:',i))
}
dev.off()
rm(estDisp.ls, i)
gc() #free some memory 


## QC - pvalue and LFC distribution ### -----------------------------------------------------
message("pvalue and LFC distribution ...")
while (!is.null(dev.list()))  dev.off()
n = length(condition)-1
pdf(file = paste0("Plots/DataQC/",substance,"_pvalue_LFC_distr.pdf"),
    width = 3*n, height = 13.33, onefile = T, bg = "transparent")
#png(filename = paste0("Plots/DataQC/",substance,"_pvalue_LFC_distr.png"),
#    width = 6*n, height = 18, units = "cm", bg = "white", pointsize = 7, res = 450)
par(mfrow=c(4,n))
## pval distr
for(i in names(res.ls)){
  hist(res.ls[[i]]$pvalue, main= paste(substance,i),col = "gray50", 
       border = "gray50", xlab = "p-value", breaks = 40)
}
## pval vs padj & qval

for(i in names(res.ls)) { pVSpadjPlot(res.ls[[i]], title = i) }
## Log2FC distr.

for(i in names(res.ls)) { lfcDistPlot(res.ls[[i]], title = i) }
for(i in names(res.ls)) {
  lfcut = quantile(abs(res.ls[[i]]$log2FoldChange),Qt)
  lfcDistPlot(resLfs.ls[[i]], title = i, LFcut = lfcut, Ylim = c(0,200), xlab = "apgelm(lfc)")
  }
dev.off()
rm(n,i,lfcut)


## QC - Norm. gene count transformation -------------------------------------------------
message("Norm. count matrix transformation ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_NormCount_Transformation.pdf"),
    width = nrow(coldata)-2, height = 4, onefile = T, bg = "transparent")
par(mfrow=c(1,3))
boxplot(normMtx$norm , notch = TRUE,las = 3, #rotating sample labels 90
        main = "Normalized read counts", ylab = "norm. read counts")
boxplot(normMtx$ntd, notch = TRUE, las = 3, #rotating sample labels 90
        main = "log2 Transformation", ylab = "log2 (norm. read counts + 1)")
boxplot(normMtx$nrl.bl, notch = TRUE, las = 3, #rotating sample labels 90
        main = "rlog Transformation", ylab = "rlog (norm. read counts)")
dev.off()


## QC - PearsonCor of biol. replicates ----------------------------------------------------
message("Biological replicates correlation ...")
# get all possible combinations of samples for each condition
comb = lapply(id.ls, function(x){gtools::combinations(length(x),2,x)})
# plot dim.
w = nrow(comb[[1]])   #width
h = length(condition) #height


## log2 ##
gg = multiCorPlot(normMtx$ntd, label = "log2(norm.counts)")
while (!is.null(dev.list()))  dev.off()
png(filename = paste0("Plots/DataQC/",substance,"_Correlation_log2.png"),
    width = 6.33*w, height = 6.6*h, units = "cm", bg = "white", 
    pointsize = .5, res = 350)
print(ggpubr::ggarrange(plotlist = gg, ncol = w, nrow = h, 
                        labels = LETTERS[1:length(gg)]))
dev.off()
## rlog ##
gg = multiCorPlot(normMtx$nrl.bl, label = "rlog(norm.counts)")
while (!is.null(dev.list()))  dev.off()
png(filename = paste0("Plots/DataQC/",substance,"_Correlation_rlog.png"),
    width = 6.33*w, height = 6.6*h, units = "cm", bg = "white", 
    pointsize = .5, res = 350)
print(ggpubr::ggarrange(plotlist = gg, ncol = w, nrow = h, 
                        labels = LETTERS[1:length(gg)]))
dev.off()
# clear env
rm(w,h,gg,comb)
gc()


## QC - SD varriance plot  -------------------------------------------------


message("Plotting standard deviation composition ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_meanSD.pdf"), #print directly to pdf
    width = 17, height = 9.6)
my_meanSdPlot(normMtx)
dev.off()
message("Done!\n")


## Sample Distance Matrix  --------------------------------------------------
message("Sample distance matrix plotting ...")

mtx    = normMtx$nrl.bl #Input for disMtx
topVar = c(2500,5000,nrow(mtx)) # Number top N varying genes/proteins

while(!is.null(dev.list()))  dev.off()
pdf(file = paste0(home,"/Plots/",substance,"_SampleDistMtx.pdf"), width = 8, height = 6.5)
for(dis in c("euclidean","canberra","manhattan")){
  message(paste0("\nComputing Sample Dist with dist measure:\t",dis))
  for(i in topVar){ myDismtx(mtx, method=dis, top=i, title=substance) }
}
dev.off()
rm(i,dis,topVar,mtx)
message("Done!\n")



## PCA & t-SNE  -----------------------------------------------------------------
message("PCA and t-SNE clustering ...")
# Add other k-mean cluster related type => See ENSEMBL Course


topVar = c(500,2000,10000,nrow(normMtx$nrl)) # Number top N varying genes/proteins
gg = list() #to store plots in
for(i in topVar){ # ref.norm - ref channel & global mean normalized
  gg[[paste0("pca",i)]]  <- myPCA(normMtx$nrl.bl, top=i, title=paste(substance,"rlog[norm].bl"))
  gg[[paste0("tsne",i)]] <- mytSNE(normMtx$nrl.bl, top=i, title=paste(substance,"rlog[norm].bl")) }
gg1 = list() #to store plots in
for(i in topVar){ # norm - only global mean normalized
  gg1[[paste0("pca",i)]]  <- myPCA(normMtx$nrl, top=i, title=paste0(substance,"rlog[norm]"))
  gg1[[paste0("tsne",i)]] <- mytSNE(normMtx$nrl, top=i, title=paste0(substance,"rlog[norm]")) }
gg2 = list()
for(i in topVar){
  gg2[[paste0("pca",i)]]  <- myPCA(normMtx$nrl.bl, top=i, title=paste(substance,"rlog[norm].bl"), IDlabs=T)
  gg2[[paste0("tsne",i)]] <- mytSNE(normMtx$nrl.bl, top=i, title=paste(substance,"rlog[norm].bl"), IDlabs=T) }
gg3 = list()
for(i in topVar){
  gg3[[paste0("pca",i)]]  <- myPCA(normMtx$nrl, top=i, title=paste(substance,"rlog[norm]"), IDlabs=T)
  gg3[[paste0("tsne",i)]] <- mytSNE(normMtx$nrl, top=i, title=paste(substance,"rlog[norm]"), IDlabs=T) }
# Export to pdf
while(!is.null(dev.list())) dev.off()
pdf(file = paste0(home,"/Plots/",substance,"_PCA.tSNE.rlogNormCounts.pdf"),
    width = 12, height = length(topVar)*4.2 )
print(ggpubr::ggarrange(plotlist=gg, ncol=2, nrow=length(topVar)))
print(ggpubr::ggarrange(plotlist=gg2, ncol=2, nrow=length(topVar)))
print(ggpubr::ggarrange(plotlist=gg1, ncol=2, nrow=length(topVar)))
print(ggpubr::ggarrange(plotlist=gg3, ncol=2, nrow=length(topVar)))
dev.off()
message("Done!\n")
rm(gg,gg1,gg2,gg3,topVar,i)



## MA plots & Vulcano plots -----------------------------------------------------
message("Plotting MA and Vulcano plots ...")


# plot
ggMA <- list()
ggVU <- list()
while (!is.null(dev.list()))  dev.off()
#pdf(file = paste0("Plots/",substance,"_MAplots_Vulcano.pdf"), width = 5*length(res.ls), height = 10)

# res - non-shrunk
for(k in names(res.ls)){
  ggMA[[k]] <- MAfun(res.ls[[k]], title = paste(substance,k))
  ggVU[[k]] <- VulcFun(res.ls[[k]], title = paste(substance,k))
}
message("Printing non-shrunk results to png ...")
png(file = paste0("Plots/",substance,"_MAplots_Vulcano.png"),
    width = 8.5*length(res.ls), height = 20, units = "cm", res = 350)
print(ggpubr::ggarrange(plotlist = c(ggMA,ggVU),ncol = length(res.ls), nrow =2))
dev.off()

# resLfs - lfc shrunk
stopifnot(all(names(res.ls) == names(resLfs.ls)))
for(k in names(resLfs.ls)){
  ggMA[[k]] <- MAfun(resLfs.ls[[k]], title = paste(substance,k),
                     LFcut = quantile(abs(res.ls[[k]]$log2FoldChange),Qt))
  ggVU[[k]] <- VulcFun(resLfs.ls[[k]], title = paste(substance,k),
                       LFcut = quantile(abs(res.ls[[k]]$log2FoldChange),Qt))
}
message("Printing apeglm shrunk results to png ...")
png(file = paste0("Plots/",substance,"_MAplots_Vulcano.lfs.png"),
    width = 8.5*length(resLfs.ls), height = 20, units = "cm", res = 350)
print(ggpubr::ggarrange(plotlist = c(ggMA,ggVU),ncol = length(res.ls), nrow =2))
dev.off()

rm(ggMA, ggVU, k)
message("Done!\n")


## DEG Correlation pairwise comparison -----------------------------------------

# Add a gene column & replace 'log2FoldChange' with 'log2FC'; to fitt upper function
mod_res <- function(res.ls){
  RES <- lapply(res.ls, function(res){
    x = res %>% as.data.frame()
    x$Gene <- row.names(x)
    colnames(x)[grep("log2FoldChange", colnames(x))] <- "log2FC"
    return(x)
  })
  return(RES)
}

## Run the next lines only when more than one condition is present!!!! ##
if(length(condition) > 2){
  message("\nStart shiny DEG Correlation & Venn plotting ...")
  
  gg <- list() #to store plots in
  RES = mod_res(res.ls)
  RESlfs = mod_res(resLfs.ls)
  
  # No LFcut
  for(k in names(RES)[-length(RES)]){
    #n.x = gsub("[.].+$","", tail(names(RES),1))
    #n.y = gsub("[.].+$","", k)
    n.x = tail(names(RES),1)
    n.y = k
    message(paste("Plotting\t", n.y, "vs", n.x))
    tmp = DE.Cor.Venn(RES[[length(RES)]], RES[[k]], LFcut.x=0, LFcut.y=0, title = paste0("padj<",p),
                      typ.x = n.x, typ.y = n.y)
    gg[[paste0(n.y,".VE")]]  <- tmp$VE
    gg[[paste0(n.y,".COR")]] <- tmp$COR
  }
  
  # LFcut
  for(k in names(RES)[-length(RES)]){
    #n.x = gsub("[.].+$","", tail(names(RES),1))
    #n.y = gsub("[.].+$","", k)
    n.x = tail(names(RES),1)
    n.y = k
    message(paste("Plotting\t", n.y, "vs", n.x))
    tmp = DE.Cor.Venn(RES[[length(RES)]], RES[[k]], title = paste0("padj<",p," & LFcut (",(1-Qt)*100,"%Q)"),
                      typ.x = n.x, typ.y = n.y)
    gg[[paste0(n.y,".VEcut")]]  <- tmp$VE
    gg[[paste0(n.y,".CORcut")]] <- tmp$COR
  }
  
  # LFcut on apeglm shrunk lfc values
  for(k in names(RES)[-length(RES)]){
    #n.x = gsub("[.].+$","", tail(names(RES),1))
    #n.y = gsub("[.].+$","", k)
    n.x = tail(names(RES),1)
    n.y = k
    message(paste("Plotting\t", n.y, "vs", n.x))
    LFcut.x = quantile(abs(RES[[length(RES)]]$log2FC),Qt)
    LFcut.y = quantile(abs(RES[[k]]$log2FC),Qt)
    X = RESlfs[[length(RESlfs)]]
    Y = RESlfs[[k]]
    s.x = droplevels(subset(X[X$padj <= p, ], abs(log2FC) >= LFcut.x))[,"Gene"]
    s.y = droplevels(subset(Y[Y$padj <= p, ], abs(log2FC) >= LFcut.y))[,"Gene"]
    tmp = DE.Cor.Venn(RES[[length(RES)]], RES[[k]], title = paste0("padj<",p," & LFcut (",(1-Qt)*100,"%Q) on lfs"),
                      LFcut.x = LFcut.x, LFcut.y = LFcut.y,
                      typ.x = n.x, typ.y = n.y, lfsSelect = T,
                      lfsSet.x = s.x, lfsSet.y = s.y)
    gg[[paste0(n.y,".VEcut.LFS")]]  <- tmp$VE
    gg[[paste0(n.y,".CORcut.LFS")]] <- tmp$COR
  }
  
  # LFcut on apeglm shrunk lfc values with shrunk lfc
  for(k in names(RES)[-length(RES)]){
    #n.x = gsub("[.].+$","", tail(names(RES),1))
    #n.y = gsub("[.].+$","", k)
    n.x = tail(names(RES),1)
    n.y = k
    message(paste("Plotting\t", n.y, "vs", n.x))
    tmp = DE.Cor.Venn(RESlfs[[length(RESlfs)]], RESlfs[[k]], title = paste("padj <",p,"with lfs values"),
                      LFcut.x = 0, LFcut.y = 0,
                      typ.x = n.x, typ.y = n.y)
    gg[[paste0(n.y,".VEcut.LFS2")]]  <- tmp$VE
    gg[[paste0(n.y,".CORcut.LFS2")]] <- tmp$COR
  }
  
  if(length(gg) > 0){
    while (!is.null(dev.list()))  dev.off()
    pdf(file = paste0("Plots/",substance,"_DEG_CorrVenn.pdf"),
        width = 12, height = 6.1*(length(res.ls)-1) )
    print( ggpubr::ggarrange(plotlist = gg, ncol = 2, nrow = length(res.ls)-1) )
    dev.off()
    message("Done!\n")
    rm(RES,RESlfs,tmp,gg,k,n.x,n.y,LFcut.x,LFcut.y,X,Y,s.x,s.y)
    gc()
  } else {
    message("No DEG correlation plot was plotted as no DEGs were found!\n")
  }
}


## Multi Venn Plot  --------------------------------------------------------


if(length(condition) > 2){
  if(length(condition)-1 > 6){
    warning("Dataset contains ",length(condition)-1," treatment conditions.\n",
            "Multi venn diagrams can only be drawn for a maximum of 6 treatments.\nSkipping multi venn plotting.")
  } else {
    message("Multi Venn plotting ...")
    # use deg.ls ; degCut.ls & degLfs.ls as input
    while (!is.null(dev.list()))  dev.off()
    pdf(file = paste0("Plots/",substance,"_VennPlots.pdf"),
        width = 13, height = 11)
    par(mfrow = c(2,2))
    gg <- list()
    for(i in c("deg.ls","degCut.ls","degLfs.ls")){
      deg <- lapply(get(i), row.names)
      if(i == "deg.ls"){n = paste0("padj<",p)}
      if(i == "degCut.ls"){n = paste0("padj<",p," & LFcut(",(1-Qt)*100,"% Qt)")}
      if(i == "degLfs.ls"){n = paste0("padj<",p," & LFcut(",(1-Qt)*100,"% Qt) on lfs")}
      # Venn1
      venn::venn(deg, ilab=TRUE, zcolor = ann_colors$Condition[-1], lty = 0, 
                 ilcs = 1, sncs = 1, box = F)
      text(0,1000, labels = paste(substance,"DEGs [",n,"]"), pos = 4)
      # Venn 2
      gg[[i]] = myVenn(deg, paste(substance,"DEGs\n",n))
    }
    print(ggpubr::ggarrange(plotlist = gg, ncol = 2, nrow = 2))
    dev.off()
    rm(i,n, gg, deg)
  }
  message("Done!\n")
}


# Session Information & save Rdata  --------------------------------------
sink(paste0(home,"/DESeq2_SessionInfo.txt"))
print(date())
print(devtools::session_info())
sink()
message("Saving RData object. This might take a while ...")
save.image(paste0(home,"/DESeq2_pairwise.RData"))
message(paste0("\n\nFinished Wald's pairwise testing DESeq2's DEG analysis & data plotting.",
               "\nAll outputs were saved in:\n",home,
               "\n\nWhat a ride! Finally END of SCRIPT! :)\nJ.A.R.V.I.S. over and out\n"))
setwd("../")
gc()
###    END OF SCRIPT   ####
