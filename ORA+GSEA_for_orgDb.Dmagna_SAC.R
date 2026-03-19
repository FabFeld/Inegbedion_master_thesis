#######################################################################
###  ORA & GSEA for GO terms via ClusterProfiler on DESeq2 Results  ###
#######################################################################

# Author:   Hannes Reinwald
# Contact:  hannes.reinwald@bayer.com

### README: ### -------------------------------------------------------------------------------
# This script performs Over Represenation Analysis (ORA) with a set of predefined diff. expressed
# genes / proteins against the background of all commonly detected genes/proteins &
# Gene Set Enrichment Analysis (GSEA) based on the sorted log2-FCs. 
# Here we use ClusterProfiler & ReactomePA for this purpose. For details please see: 
# https://yulab-smu.github.io/clusterProfiler-book/index.html
# To increase the power of GSEA we recommend to use apeglm shrunk Effect sizes (log2-FC) for the analysis

# More details about the principals of GSEA:
# https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/

# Input ORA:
# Data table list of DEGs/DEPs with ENSEMBL/PROT IDs, apeglm(LFC), ENTREZ, padj (=output tables of DESeq2 => lfsFCcut)
# Universe - set of background genes (all detected genes in the experiment) => This background can be used for GSEA then (=any result table of DESeq2 output)

# Input GSEA: 
# Data table with ENSEMBL / PROT IDs, apeglm(LFC), ENTREZ, padj 
# Transcriptomics - Ideally the result table from the DESeq2 with apeglm shrunk LFC (named with _reslfs) and annotated ENTREZ IDs
# Proteomics - I recommend to shrinking here as well on the LFC tables.

###############

install.packages("devtools", type = "binary")
###  PACKAGES  ### --------------------------------------------------------------------------
require(dplyr)
require(clusterProfiler)
require(AnnotationDbi)
require(ggplot2)
require(ggridges)
require(ggnewscale)
##################


###  PARAMETERS  ### -------------------------------------------------------------
## Working directory = project directory containing the DESeq2 output folder
#setwd("C:/Users/GNVKY/OneDrive - Bayer/Projects/F-IME_daphniaOMICs/data/Dmagna_RNAseq/fipronil/") # <-- rmv this line later
HOME = getwd()

## OrgDb 
ORG = "org.Dmagna.eg.db" # <-- specify the species orgDb you are working with.
#ORG = "org.Cdipterum.eg.db"

require(ORG, character.only = T) # <-- load respective orgDb package
KEY = "GID" #<-- for custom orgDb packages it is called "GID"
ONT = c("BP","MF","CC") # <-- Gene Ontologies to analyse. 
padj.cutoff = 0.1       # <-- padj cutoff for filtering the enriched results
####################


###  DATA IMPORT  ###
# Import DEG tables (apgelm shrunk lfcut & pcut) for ORA ---------------------------------------
## import DEGs for the three different criteria levels
files = list.files(path = "./DESeq2_Pairwise/DEGs", pattern = "[.]csv$", full.names = T, recursive = T)
# extract substance name/s
subs = gsub("^.+[/]","",files) %>% sub("_.+$","",.) %>% unique %>% paste(., collapse = "_")

inputORA = list()
for(i in c("[/]lfsFCcut[/]","[/]unshrinkFCcut[/]","[/]unshrink[/]")){
  f = grep(i,files, value = T)
  if(length(f) == 2){f = f[c(2,1)]}
  if(length(f) == 3){f = f[c(2,3,1)]}
  ls = lapply(f, read.csv2, header=T, row.names=1)
  names(ls) = gsub("^.+[/]","",f) %>% sub("_.+_","_",.) %>% sub("[.]csv$","",.)
  inputORA[[gsub("\\[/\\]", "", i)]] = ls }
rm(i,f,ls)
inputORA = unlist(inputORA, recursive = F)

# Import complete res tables (apgelm shrunk lfcut & pcut) for GSEA & universe definition -------
files = list.files(path = "./DEseq2_Pairwise/Results/", pattern = "_reslfs_.+[.]csv$", full.names = T)
if(length(files) == 2){files = files[c(2,1)]}
if(length(files) == 3){files = files[c(2,3,1)]}
inputGSEA = list()
for(i in files){
  x = read.csv2(i,header = T, row.names = 1)
  x = x[which(!(is.na(x$pvalue))),]# rmv genes with no computed pvalue
  name = gsub("^.+[/]","",i) %>% gsub(".csv","",.) %>% gsub("_reslfs_","_",.)
  inputGSEA[[name]] = x }
rm(name,x,i)

# Unify tables from protein and transcript output ---------------------------------------
renameFun = function(ls,string="",replace=""){
  tmp = lapply(ls, function(x){
    colnames(x)[which(colnames(ls[[1]]) %in% string)] <- as.character(replace)
    x })
  tmp }

# Replace "log2FoldChange" with "log2FC"
if(any(colnames(inputGSEA[[1]]) %in% "log2FoldChange")){
  inputGSEA = renameFun(inputGSEA,"log2FoldChange","log2FC")
  inputORA  = renameFun(inputORA,"log2FoldChange","log2FC") }
# Replace "adj.pvalue" with "padj"
if(any(colnames(inputGSEA[[1]]) %in% "adj.pvalue")){
  inputGSEA = renameFun(inputGSEA,"adj.pvalue","padj")
  inputORA  = renameFun(inputORA,"adj.pvalue","padj") }
# To ensure a better compatibility between AnnotaionDbi and biomaRt annotations we will
# change the attributes label accordingly. 
# Replace "external_gene_name" with "SYMBOL"
if(any(colnames(inputGSEA[[1]]) %in% "external_gene_name")){
  inputGSEA = renameFun(inputGSEA,"external_gene_name","SYMBOL")
  inputORA  = renameFun(inputORA,"external_gene_name","SYMBOL") }
# Replace "external_gene_name" with "SYMBOL"
if(any(colnames(inputGSEA[[1]]) %in% "entrezgene_id")){
  inputGSEA = renameFun(inputGSEA,"entrezgene_id","ENTREZID")
  inputORA  = renameFun(inputORA,"entrezgene_id","ENTREZID") }


# Define common gene sets = universe (background for ORA) -------------------------------
# Common gene set definition only really needed when comparing results from different experiments
# as identified proteins / genes could differ between experiments
# But we will still do it anyways to make this script applicable for a wider usage.
univ = lapply(inputGSEA, function(x){row.names(x)}) %>% Reduce(intersect, .)


# Create FINAL INPUT tables GSEA / ORA --------------------------------------------------
# (with common gene sets (univ)) for GSEA & ORA (name & log2FC)
# 1) select common set (element of universe)
# 2) sort by pvalue then remove duplicates (needed for entrezid)
#    (That way the most signif. result is keept for downstream analysis)
# 3) sort bei log2FC and export lfc values with names into ls object
InputGenerator = function(df, entrezID = F){
  if(entrezID == T){
    x = df[df$ENTREZID %in% univ.entrez,]
    x = dplyr::arrange(x, padj)
    x = x[!duplicated(x$ENTREZID),] #rmv duplicated entrez ids
  } else { # for ensembl gene IDs
    x = df[row.names(df) %in% univ,]
  }
  x = dplyr::arrange(x, desc(log2FC))
  tmp = x$log2FC
  if(entrezID == T){
    names(tmp) = as.character(x$ENTREZID)
  } else {
    names(tmp) = as.character(row.names(x))
  }
  tmp = na.omit(tmp)
  tmp[names(tmp) == ""] <- NA
  tmp = tmp[!is.na(names(tmp))] #remove any entry without an ID
  setNames(tmp, names(tmp))
}

# ENSEMBL / ENSEMBLPROT for GO
GSEA = lapply(inputGSEA, InputGenerator)
ORA  = lapply(inputORA,  InputGenerator) #%>% lapply(., names)

# split ORA ls back into different sets
ls = list()
for(i in unique(sub("[.].+$","",names(ORA))) ){ 
  ls[[i]] = ORA[grep(paste0(i,"[.]"),names(ORA))] }
ORA = ls
rm(ls,i)
######################


###############
###   ORA   ###
###############
## Run ORA clusterComparison ## -----------------------------------------------------------------
# function with function parameters
my_ORA = function(ora.ls, univ, ORG, KEY="GID", ont="BP", sim.method="Wang", simpl=T, cutoff=.8, 
                  pcut=.1, qcut=.1, padjm="BH", MIN=10, MAX=500){
  message("\nStarting ORA via 'enrichGO()' for GO.",ont," terms ...")
  # check if input is only character list or a named list of log2-fc values
  if(!is.character(ora.ls[[1]][1])) {ora.ls = lapply(ora.ls, names)}
  
  GO = clusterProfiler::compareCluster(
    ora.ls, OrgDb = ORG, keyType = KEY, ont = ont, fun = "enrichGO", universe = univ, 
    pvalueCutoff = pcut, qvalueCutoff = qcut, pAdjustMethod = padjm,
    minGSSize = MIN, maxGSSize = MAX)
  
  if(!is.null(GO)){
  #if(!is.null(GO) & nrow(GO@compareClusterResult) > 0){ # <-- this bugs if GO == NULL
    if(nrow(GO@compareClusterResult) > 0){
      # Compute semantic similarities among GO terms:
      message("Computing semantic similarities among GO terms. This may take some time ... ")
      d = GOSemSim::godata(ORG, ont=ont, computeIC=FALSE, keytype = KEY)
      GO = enrichplot::pairwise_termsim(GO, method = sim.method, semData = d)
      if(simpl){
        message("Simplifying ORA result's GO terms. This may take some time ... ")
        # Rmv GO terms with redundant biological information <-- this takes forever!!! 
        GO = simplify(GO,cutoff) }
      # resort results after pvalues
      GO@compareClusterResult = GO@compareClusterResult[order(GO@compareClusterResult$pvalue),]
      message("Done!\n")
      return(GO)
    }else{ message("No enriched terms were identified.\n") }
  }
  
}


ORA.res = lapply(ORA, function(x){ 
  res = list()
  for(ont in ONT){
    res[[ont]] = my_ORA(x, univ, ORG, KEY, ont = ont, qcut = padj.cutoff) }
  return(res)
})
n = lapply(ORA.res, length) %>% unlist > 0 # rmv all objects from list which are empty
ORA.res = ORA.res[n]
gc() # clear memory

## Export ORA results ## ------------------------------------------------------
# Create Output folders in "clusterProfiler" 
dir.create("clusterProfiler", showWarnings = F)
dir.create("clusterProfiler/ORA", showWarnings = F)
# retrieve data frames to a new list object; we might need for plotting later
ORA.resDf = list()
for(i in names(ORA.res)){
  ORA.resDf[[i]] = lapply(ORA.res[[i]], function(res){ res@compareClusterResult })
}
# write result tables to csv
message("Exporting ORA result tables ...")
for(i in names(ORA.resDf)){
  out = paste0("clusterProfiler/ORA/",i)
  dir.create(out, showWarnings = F)
  for(k in names(ORA.resDf[[i]])){
    write.csv2(ORA.resDf[[i]][[k]], paste0(out,"/","GO.",k,"_",subs,"_",i,"_compareCluster.csv"), row.names = F)
  }
}
message("Done!\n")

## Plot ORA results ## --------------------------------------------------------
message("Start plotting ORA results ...")

## EMAP-PLOT - Enrichment map ##
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html?q=emapplot#enrichment-map
for(i in names(ORA.res)){
  ls = ORA.res[[i]]
  for(k in names(ls)){
    message("Plotting ORA Enrichment Map via 'emapplot()' for: ",i," - GO.",k)
    res = ls[[k]]
    gg = list()
    for(n in c(10,20,30)){# plot the top 10, 20, 30 categories
      x = length(unique(res@compareClusterResult$Description))
      if(x > n){ x = n }
      set.seed(42)
      gg[[paste0("top.",x)]] = emapplot(
        res, x, shadowtext=F, repel=T, pie="Count", legend_n=3, layout = "kk") +
        labs(subtitle = paste0(subs," - ORA Enrichment Map \nGO.",k," - Top: ",x))
      
      gg[[paste0("top.",x,".gr")]] = emapplot(
        res, x, shadowtext=F, repel=T, pie="Count", legend_n=3, layout = "kk",
        group_category=T, group_legend=T) +
        labs(subtitle = paste0(subs," - ORA Enrichment Map \nGO.",k," - Top: ",x))
    }
    ## rmv buggy plots
    message("Checking plots ...")
    tmp = list()
    for(p in names(gg)){
      z = try(print(gg[p]), silent = T)
      if(inherits(z, "try-error")){tmp[[p]] = ggplot()} else {tmp[[p]] = gg[[p]]}
    }
    
    if(length(tmp) > 0){
      while (!is.null(dev.list()))  dev.off() # plotting ...
      paste0(HOME,"/clusterProfiler/ORA/",i,"/GO.",k,"_",subs,"_",i,"_emapplot.pdf") %>%
        pdf(., width = 20, height = 9)
      print( ggpubr::ggarrange(plotlist = tmp, nrow = 1, ncol = 2) )
      dev.off()
    }
  }
}

## CNET-PLOT - Gene-Concept Network ##
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html?q=cnetplot#compare-cnetplot 
for(i in names(ORA.res)){
  ls = ORA.res[[i]]
  fc = ORA[[i]]
  for(k in names(ls)){
    message("Plotting ORA Gene-Concept Network via 'cnetplot()' for: ",i," - GO.",k)
    res = ls[[k]]
    gg = list()
    for(n in c(5,10,20)){ # plot the top 5, 10, 15 categories
      x = length(unique(res@compareClusterResult$Description))
      if(x > n){ x = n }
      set.seed(42)
      
      gg[[paste0("top.",x)]] = cnetplot(
        res, x, shadowtext="none", repel=T, pie="Count",
        legend_n=3, layout = "kk") + 
        labs(subtitle = paste0(subs," - ORA Gene-Concept Network \nGO.",k," - Top: ",x))
      
      gg[[paste0("top.",x,".2")]] = cnetplot(
        res, x, shadowtext="none", repel=T, pie="Count",
        legend_n=3, layout = "kk", node_label="category") +
        labs(subtitle = paste0(subs," - ORA Gene-Concept Network \nGO.",k," - Top: ",x))
    }
    
    if(length(gg) > 0){
      while (!is.null(dev.list()))  dev.off() # plotting ... 
      paste0(HOME,"/clusterProfiler/ORA/",i,"/GO.",k,"_",subs,"_",i,"_cnetplot.pdf") %>%
        pdf(., width = 25, height = 11)
      print( ggpubr::ggarrange(plotlist = gg, nrow = 1, ncol = 2) )
      dev.off()
    }
  }
}

## DOT-PLOT ## 
# Custom bade function based on the idea from dotplot()
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html?q=dotplot#compare-dotplot

# The below function will create a dotplot from your ORA analysis results created with
# compareCluster(). The function will:
# 1) Filter the results table for signif results based on the provided pcut threshold for
#    for either pvalue or p.adjust.(DEFAULT: 0.05 for p.adjust)
# 2) Select the top N significantly enriched results for each Cluster. Top N selection either 
#    based on pvalue or p.adjust (DEFAULT: top 10 for pvalue)
# 3) Will sort the order of terms plotted after significance for each Cluster. This gives a
#    much better visualization of the results as it can be instantly seen if different clusters
#    have distinct enriched terms or more overlap among them.
# NOTE: The 'res' object must be a data frame: res <- as.data.frame(compClustRes)
ggDot4ORA = function(res, top=10, topBy="pvalue", pcut=.1, minCount=3, minGeneR = .02, 
                     filter="p.adjust", title="", clust.order = NULL, show.enrichment = T){
  # clust.order: A vector describing the order in which to display the gene group clusters from left to right. 
  require(ggplot2)
  if(filter != "p.adjust" & filter != "pvalue"){ stop("'filter' must be specified as 'p.adjust' or 'pvalue'!") }
  if(topBy != "p.adjust" & topBy != "pvalue"){ stop("'filter' must be specified as 'p.adjust' or 'pvalue'!") }
  stopifnot(is.data.frame(res))
  
  # subset df based on filter features 
  # For details of GeneRatio and BgRatio see: https://www.biostars.org/p/220465/ 
  df = subset(res, res[,filter] <= pcut) %>% subset(., Count >= minCount) #%>% subset(., GeneRatio >= minGeneR)
  if(nrow(df) > 0){
    # compute "Bg ratio"
    M = as.integer(sub("/.+$","",df$BgRatio)) 
    N = as.integer(sub("^.+/","",df$BgRatio))
    stopifnot(length(M) == length(N))
    df$BgR = M/N # BgR = M/N
    # compute "gene ratio"
    k = df$Count
    n = as.integer(sub("^.+/","",df$GeneRatio))
    stopifnot(length(k) == length(n))
    df$GeneRatio.ori = df$GeneRatio # stil keep the initial gene ratio format
    df$GeneRatio = k/n # GeneR = k/n
    # get Nbr of all DEGs identified for a GO term set in a cluster
    df$Cluster = paste0(df$Cluster," (",n,")")
    
    ## Add Enrichment column (percent enrichment = k / M) ##
    df$Enrichment = (k/M) * 100 
    
    # Last subset
    df = subset(df, df$GeneRatio >= minGeneR)
  }
  
  if(nrow(df) > 0){
    # subset df back into cluster groups and pick the top N GO terms from each Cluster
    ls = split(df, df$Cluster)
    # # This will give an order of sample group clusters with decreasing Nbr of enriched terms
    if(is.null(clust.order)){ clust.order = lapply(ls, nrow) %>% unlist %>% sort(., decreasing = T) %>% names }
    ls = ls[clust.order]
    ls = ls[unlist(lapply(ls, nrow)) > 0] %>% lapply(., function(x){ if(nrow(x)>0){x[order(x[,topBy]),]} }) %>%
      # Select top N terms from each cluster
      lapply(., function(x){ if(top > nrow(x)){x = x[1:nrow(x),]}else{x = x[1:top,]} }) %>%
      # sort df in ls after pvalue or GeneRatio to have a nicer order for plotting of description terms
      #lapply(., function(x) x[order(x$pvalue),])
      lapply(., function(x) x[order(x$GeneRatio, decreasing = T),])
    
    # Combine to single vector for ordering the enriched terms in plot
    DescrOrd = lapply(ls, function(x) x$Description) %>% unlist() %>% unique() %>% rev()
    
    ## Build final df for plotting!
    df1 = df[df$Description %in% DescrOrd, ]
    df1$Description = factor(stringr::str_wrap(df1$Description, 70), 
                             levels =  stringr::str_wrap(DescrOrd, 70))
  
    # compute GeneRatio over bg-ratio
    #df1$genR.bgR <- df1$GeneRatio / df1$bgR
    g = ggplot(df1, aes(x=GeneRatio, y=Description, color=df1[,filter], size=Count)) + geom_point() +
      #ggplot(df1, aes(x=geneR, y=Description, color=df1[,filter], size=Count)) + geom_point() +
      scale_colour_gradient2(low = "firebrick2", mid = "yellow", high = "#0D0887FF", midpoint = .05) +
      labs(col= paste0(filter)) +
      facet_grid(~Cluster) + theme_light() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      labs(title = paste("ORA -",title,":",length(unique(df1$ID)),"terms"),
           subtitle = paste(filter,"<",pcut,"/ Top:",top,"by",topBy,
                            " \nmin.GeneRatio:",minGeneR,"/ min.Count:",minCount)) + ylab(NULL)
    if(show.enrichment){ g = g + aes(x=Enrichment) + xlab("GeneSet Enrichment [%]")}
    return(g)
  } else { message("Nothing to plot based on the provdied filter criteria :(") }
}

message("Creating dotplot ...")
gg = list()
for(i in names(ORA.resDf)){
  tmp = ORA.resDf[[i]]
  for(k in names(tmp)){
    res = tmp[[k]]
    if(nrow(res) < 25){ top = nrow(res) }else{ top = 25 }
    gg[[paste0(i,".",k)]] = ggDot4ORA(res,top, pcut = .05, title = paste0("GO.",k,"_",subs," (",i,")"))
    gg[[paste0(i,".",k,"2")]] = ggDot4ORA(res,top, pcut = padj.cutoff, title = paste0("GO.",k,"_",subs," (",i,")"))
  }
}

if(length(gg) > 0){
  gg = gg[sort(names(gg))]
  for(i in unique(sub("[.].+$","",names(gg)))){
    g = gg[grep(paste0(i,"[.]"),names(gg), value = T)]
    while (!is.null(dev.list()))  dev.off() # plotting ... 
    paste0(HOME,"/clusterProfiler/ORA/",i,"/GO.ALL_",subs,"_",i,"_dotplot.pdf") %>%
      #paste0(HOME,"/clusterProfiler/ORA/GO.ALL_",subs,"_dotplot.pdf") %>%
      pdf(., width = 9, height = 7)
    print( g )
    dev.off()
  }
}

rm(ls,res,gg,x,i,k,n,tmp,p) # clean up 
gc()
################


################
###   GSEA   ###
################
## Run GSEA ## -----------------------------------------------------------------
# function with function parameters
my_gseGO = function(gse.ls, ORG, KEY="GID", ont="BP", sim.method="Wang", simpl=T, cutoff=.8, 
                    qcut=.2, padjm="BH", MIN=10, MAX=500, perm = 1000){
  message("\nStarting GSEA via 'gseGO()' for GO.",ont," terms ...")
  GO = gseGO(gse.ls, OrgDb = ORG, keyType = KEY, ont = ont, minGSSize = MIN, maxGSSize = MAX,
             pvalueCutoff = qcut, pAdjustMethod = padjm, seed = T, by = "fgse", nPermSimple = perm) # <- large number here returns more results but longer computation time!
  
  if(!is.null(GO)){
    #if(!is.null(GO) & nrow(GO@result) > 0){ # <-- this bugs if GO == NULL
    if(nrow(GO@result) > 0){
      # Compute semantic similarities among GO terms:
      message("Computing semantic similarities among GO terms. This may take some time ... ")
      d = GOSemSim::godata(ORG, ont=ont, computeIC=F, keytype = KEY)
      GO = enrichplot::pairwise_termsim(GO, method = sim.method, semData = d)
      if(simpl){
        message("Simplifying ORA result's GO terms. This may take some time ... ")
        # Rmv GO terms with redundant biological information <-- this takes forever!!! 
        GO = simplify(GO,cutoff) }
      # resort results after pvalues
      GO@result = GO@result[order(GO@result$pvalue),]
      message("Done!\n")
      return(GO) }
  }else{ message("No enriched terms were identified.\n") }
} 

# Performing GSEA with the above custom function
GSEA.res = lapply(GSEA, function(x){ 
  res = list()
  for(ont in ONT){
    # perm: <- large number here returns more results but needs longer computation time! Ideal range: 1000 - 100000
    res[[ont]] = my_gseGO(x, ORG, KEY, ont=ont, qcut=padj.cutoff, perm=100000) }
  return(res)
})
n = lapply(GSEA.res, length) %>% unlist > 0
GSEA.res = GSEA.res[n]
gc() # clear memory


## Export ORA results ## -------------------------------------------------------
## Aggregate results to Df 
# retrieve data frames to a new list object; we might need for plotting later
ls = unlist(GSEA.res, recursive = F) %>% lapply(., function(res){ res@result })
# add cluster column 
for(i in names(ls)){ls[[i]]$Cluster = i}#sub("[.][BCM][PCF]$","",i) }
# merge respective GO ONT term results
GSEA.resDf = list() # to store output Df in
for(i in ONT){
  GSEA.resDf[[paste0("GO.",i,"_",subs)]] = ls[ grep(paste0("[.]",i,"$"), names(ls)) ] %>% 
    rlist::list.rbind(.) }
# add Gene Counts Column to the GSEA Output Table
GSEA.resDf = lapply(GSEA.resDf, function(x){
  x$Count = strsplit(x$core_enrichment, "/") %>% lapply(., length) %>% unlist
  return(x) })

# Create Output folders in "clusterProfiler" and write tables in there
message("Exporting GSEA result tables ...")
dir.create("clusterProfiler/GSEA", showWarnings = F)
for(i in names(GSEA.resDf)){
  paste0(HOME,"/clusterProfiler/GSEA/",i,"_gsea.csv") %>% write.csv2(GSEA.resDf[[i]], ., row.names = F)
}
message("Done!\n")
rm(ls,i)

## Plot GSEA results ## --------------------------------------------
message("Start plotting GSEA results ...")

# simpel export plot function 
export.PDF = function(gg, file.out, ONT=c("BP","MF","CC"), w = 8, h = 7){
  l = lapply(ONT, function(x){ gg[grep(x,names(gg))] %>% length(.) }) %>% unlist %>% max
  pdf(file.out, width = w*l, height = h)
  for(i in ONT){
    g = gg[grep(i,names(gg))]
    if(length(g) > 0){ print(ggpubr::ggarrange(plotlist = g, nrow = 1, ncol = l)) }
  }
  dev.off()
  message("Done!\n") }

#  dotplot --------------------------------------------------------
ggDot4GSE = function(res, top=10, topBy="pvalue", pcut=.1, minCount=3, clust.order = NULL, minEnrichment = 2, # min 2% enrichment
                     filter="p.adjust", title=""){
  if(filter != "p.adjust" & filter != "pvalue"){ stop("'filter' must be specified as 'p.adjust' or 'pvalue'!") }
  if(topBy != "p.adjust" & topBy != "pvalue"){ stop("'filter' must be specified as 'p.adjust' or 'pvalue'!") }
  stopifnot(is.data.frame(res))
  
  # If no 'Count' column present in df simply add it
  if(!any(colnames(res) == "Count")){ res$Count = strsplit(res$core_enrichment, "/") %>% lapply(., length) %>% unlist }
  # Compute percent Enrichment for a given annotated gene set
  if(!any(colnames(res) == "Enrichment")){ res$Enrichment = (res$Count/res$setSize)*100 } # percent of enrichment for a given annotated gene set
  # subset df based on filter features 
  df = subset(res, res[,filter] <= pcut) %>% subset(., Count >= minCount) %>% subset(., Enrichment >= minEnrichment)
  
  if(nrow(df) > 0){
    
    # subset df back into cluster groups and pick the top N-signif. GO terms from each Cluster
    ls = list()
    for(i in unique(df$Cluster)){
      x = df[df$Cluster %in% i,]
      # select topN from x
      if(nrow(x) > 0){
        x = x[order(x[,topBy]),]
        if(top > nrow(x)) { x = x[1:nrow(x),] } else { x = x[1:top,] }
      }
      ls[[i]] = x }
    # This will give an order of sample group clusters with decreasing Nbr of enriched terms
    if(is.null(clust.order)){ clust.order = lapply(ls, nrow) %>% unlist %>% sort(., decreasing = T) %>% names }

    # rearrange order of description terms for plotting
    ls = lapply(ls[clust.order], function(x) x[order(-x$Enrichment),])
    t1 = lapply(ls, function(x){ x[x$NES >= 0,] }) %>% lapply(., function(x) x[order(-x$Enrichment),"Description"]) %>%
      unlist() %>% unique() %>% rev()
    t2 = lapply(ls, function(x){ x[x$NES < 0, ] }) %>% lapply(., function(x) x[order(-x$Enrichment),"Description"]) %>%
      unlist() %>% unique() %>% rev()
    DescrOrd = unique(union(t2,t1)) # Combine to single vector for ordering the enriched terms in plot
    
    ## Build final df for plotting!
    df1 = df[df$Description %in% DescrOrd, ]
    df1$Description = factor(stringr::str_wrap(df1$Description, 70), 
                             levels =  stringr::str_wrap(DescrOrd, 70))
    df1$Cluster = factor(df1$Cluster, levels = names(ls))
    
    g = ggplot(df1, aes(x=Enrichment, y=Description, color=NES, size=Count )) + geom_point() +
      scale_colour_gradient2(low = "mediumblue", mid = "white", high = "red2", midpoint = 0) +
      facet_grid(~Cluster) + theme_light() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      #aes(Count, reorder(stringr::str_wrap(Description, 80), NES))
      labs(title = paste("GSEA -",title,":",length(unique(df1$ID)),"terms"),
           subtitle = paste(filter,"<",pcut,"/ Top:",top,"by",topBy,
                            " \nmin.Enrichment:",minEnrichment,"/ min.Count:",minCount)) + ylab(NULL) + xlab("GeneSet Enrichment [%]")
    return(g)
  } else {
    message("Nothing to plot based on the provdied filter criteria :(")
  }
}

message("Creating dotplot ...")
gg = list()
for(i in names(GSEA.resDf)){
  res = GSEA.resDf[[i]]
  if(nrow(res) < 25){ top = nrow(res) }else{ top = 25 }
  gg[[paste0(i)]] = ggDot4GSE(res,top, pcut = .05, title = i)
  gg[[paste0(i,"2")]] = ggDot4GSE(res,top, pcut = padj.cutoff, title = i)
}
if(length(gg) > 0){
  gg = gg[sort(names(gg))]
  while (!is.null(dev.list()))  dev.off() # plotting ... 
  paste0(HOME,"/clusterProfiler/GSEA/GO.ALL_",subs,"_dotplot.pdf") %>%
    pdf(., width = 9, height = 7)
  print( gg )
  dev.off()
}


#  gseaplot ---------------------------------------------------
message("Creating gseaplot for QC ...")
gg = list()
for(i in names(GSEA.res)){
  tmp = GSEA.res[[i]]
  for(k in names(tmp)){
    res = tmp[[k]]
    if(nrow(res@result) < 6){x = nrow(res@result)} else {x=6}
    n = paste0("GSEA GO.",k," - ",i,"\nTop ",x," enriched Terms")
    gg[[paste0(k,".",i)]] = enrichplot::gseaplot2(
      res, 1:x, pvalue_table = T, ES_geom = "dot", title = n)
  }
}
if(length(gg) > 0){
  gg = gg[sort(names(gg))]
  while (!is.null(dev.list()))  dev.off() # plotting ... 
  paste0(HOME,"/clusterProfiler/GSEA/GO.ALL_",subs,"_gseaplot.pdf") %>%
    export.PDF(gg,., w=11, h=8)
}


#  ridgeplot ---------------------------------------------------
message("Creating GSEA ridgeplot ...")
gg = list()
for(i in names(GSEA.res)){
  tmp = GSEA.res[[i]]
  for(k in names(tmp)){
    res = tmp[[k]]
    if(nrow(res@result) < 25){x = nrow(res@result)} else {x=25}
    n = paste0("GSEA GO.",k," - ",i)
    gg[[paste0(k,".",i)]] = ridgeplot(res, showCategory = x) + xlab("log2(fold-change)") +
      labs(title=n, subtitle = paste0("Top ",x," enriched Terms"))
  }
}
if(length(gg) > 0){
  gg = gg[sort(names(gg))]
  while (!is.null(dev.list()))  dev.off() # plotting ... 
  paste0(HOME,"/clusterProfiler/GSEA/GO.ALL_",subs,"_ridgeplot.pdf") %>%
    export.PDF(gg,., w=10, h=15)
}


#  upsetplot -------------------------------
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html?q=upsetplot#upset-plot 
require(ggupset)
message("Creating GSEA upsetplot ...")
gg = list()
for(i in names(GSEA.res)){
  tmp = GSEA.res[[i]]
  for(k in names(tmp)){
    res = tmp[[k]]
    if(nrow(res@result) < 25){x = nrow(res@result)} else {x=25}
    n = paste0("GSEA GO.",k," - ",i)
    if(nrow(res@result) > 1) {
      gg[[paste0(k,".",i)]] = enrichplot::upsetplot(res, x) + ylab("log2(fold change)") + 
        geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
        labs(title=n, subtitle = paste0("Top ",x," enriched Terms"))
    }
  }
}
if(length(gg) > 0){
  gg = gg[sort(names(gg))]
  while (!is.null(dev.list()))  dev.off() # plotting ... 
  paste0(HOME,"/clusterProfiler/GSEA/GO.ALL_",subs,"_upsetplot.pdf") %>%
    export.PDF(gg,., w=12, h=9)
}

#  heatplot ------------------------------------------------------------------
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html?q=heatplot#heatmap-like-functional-classification 
message("Creating GSEA heatplot ...")
gg = list()
for(i in names(GSEA.res)){
  tmp = GSEA.res[[i]]
  for(k in names(tmp)){
    res = tmp[[k]]
    if(nrow(res@result) < 25){x = nrow(res@result)} else {x=25}
    n = paste0("GSEA GO.",k," - ",i)
    if(nrow(res@result) > 1) {
      gg[[paste0(k,".",i)]] = heatplot(res, x, foldChange = GSEA[[i]]) +
        scale_fill_gradient2(low = "mediumblue", mid = "white", high = "red2", midpoint = 0, guide = "colorbar") +
        labs(title=n, subtitle = paste0("Top ",x," enriched Terms")) }
  }
}
if(length(gg) > 0){
gg = gg[sort(names(gg))]
while (!is.null(dev.list()))  dev.off() # plotting ... 
paste0(HOME,"/clusterProfiler/GSEA/GO.ALL_",subs,"_heatplot.pdf") %>%
  export.PDF(gg,., w=11, h=7)
}

#  treeplot ---------------------------------------------------------------
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html?q=heatplot#tree-plot
message("Creating GSEA treeplot ...")
gg = list()
for(i in names(GSEA.res)){
  tmp = GSEA.res[[i]]
  for(k in names(tmp)){
    res = tmp[[k]]
    if(nrow(res@result) < 25){x = nrow(res@result)} else {x=25}
    n = paste0("GSEA GO.",k," - ",i)
    nk = ceiling(x/5)# nbr of clusters
    if(nrow(res@result) > 1) {
      gg[[paste0(k,".",i)]] = enrichplot::treeplot(
        res, hclust_method = "average", nCluster = nk) +
        labs(title=n, subtitle = paste0("Top ",x," enriched Terms")) }
  }
}
if(length(gg) > 0){
gg = gg[sort(names(gg))]
while (!is.null(dev.list()))  dev.off() # plotting ... 
paste0(HOME,"/clusterProfiler/GSEA/GO.ALL_",subs,"_treeplot.pdf") %>%
  export.PDF(gg,., w=10, h=7)
}

#  emapplot ------------------------------------
message("Creating GSEA emapplot ...")
gg = list()
val = c(10,20,30) # print the top 5, 10, 15, 20 enriched terms in a Network 
for(i in names(GSEA.res)){
  tmp = GSEA.res[[i]]
  for(k in names(tmp)){
    res = tmp[[k]]
    if(nrow(res@result) > 1){
      for(v in val){
        if(nrow(res@result) < v){v = nrow(res@result)}
        
        n = paste0("GSEA GO.",k," - ",i)
        set.seed(42)
        gg[[paste0(k,".",i,".",v)]] = emapplot(
          res, v, shadowtext=F, repel=T, layout = "kk") +
          labs(title=n, subtitle = paste0("Top ",v," enriched Terms")) 
        
        gg[[paste0(k,".",i,".",v,".2")]] = emapplot(
          res, v, shadowtext=F, repel=T, layout = "kk", group_category=T, group_legend=T) +
          labs(title=n, subtitle = paste0("Top ",v," enriched Terms")) }
    }
  }
}
if(length(gg) > 0){
  while (!is.null(dev.list()))  dev.off() # plotting ... 
  paste0(HOME,"/clusterProfiler/GSEA/GO.ALL_",subs,"_emapplot.pdf") %>%
    pdf(., width = 20, height = 9)
  for(k in ONT){
    for(i in names(GSEA)){
      g = paste0(k,".",i) %>% grep(., names(gg)) %>% gg[.]
      if(length(g) > 1){
        try( print( ggpubr::ggarrange(plotlist = g, nrow = 1 , ncol = 2) ),
             silent = T)
      }
    }
  }
}
dev.off()
message("Done!\n")


#  cnetplot -----------------------------------
message("Creating GSEA cnetplot ...")
gg = list()
val = c(5,10,20) # print the top 5, 10, 15 enriched terms in a Network 
for(i in names(GSEA.res)){
  fc = GSEA[[i]]
  tmp = GSEA.res[[i]]
  for(k in names(tmp)){
    res = tmp[[k]]
    if(nrow(res@result) > 1){
      for(v in val){
        if(nrow(res@result) < v){v = nrow(res@result)}
        n = paste0("GSEA GO.",k," - ",i)
        set.seed(42)
        
        gg[[paste0(k,".",i,".",v)]] = cnetplot(
          res, x, foldChange = fc, shadowtext="none", repel=T, layout = "kk") +
          labs(title=n, subtitle = paste0("Top ",v," enriched Terms")) 
        
        gg[[paste0(k,".",i,".",v,".2")]] = cnetplot(
          res, x, foldChange = fc, shadowtext="none", repel=T, layout = "kk", node_label="gene", colorEdge = F) +
          labs(title=n, subtitle = paste0("Top ",v," enriched Terms")) }
    }
  }
}
if(length(gg) > 0){
  while (!is.null(dev.list()))  dev.off() # plotting ... 
  paste0(HOME,"/clusterProfiler/GSEA/GO.ALL_",subs,"_cnetplot.pdf") %>%
    pdf(., width = 25, height = 11)
  for(k in ONT){
    for(i in names(GSEA)){
      g = paste0(k,".",i) %>% grep(., names(gg)) %>% gg[.]
      if(length(g) > 1){
        try( print( ggpubr::ggarrange(plotlist = g, nrow = 1 , ncol = 2) ),
             silent = T)
      }
    }
  } 
}
dev.off()
message("Done!\n")

# clean up 
rm(i,k,n,nk,v,val,x,files, fc, res,g, gg, tmp)
gc()
######################


### Save R Image & Session Info ### --------------------------
message(paste0("Saving R image under:\n",HOME,"/clusterProfiler/Pathway.RData\nThis might take a while ..."))
save.image(paste0(HOME,"/clusterProfiler/clusterProfiler.RData"))

sink(paste0(HOME,"/clusterProfiler/SessionInfo_clusterProfiler.txt"))
print(date())
print(devtools::session_info())
sink()

message("\nAll Done! END OF ORA & GSEA SCRIPT\nJ.A.R.V.I.S over and out! :)\n")
####### END OF SCRIPT #######