# Nadine Schrode
# 08/20/19
# input: Differential expression results list from DE-limma script and contrasts.fit limma object


#### LIBRARIES ####
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(wesanderson)
library(reshape2)
library(plyr)
library(plotrix)
library(fmsb)
library(calibrate)
library(variancePartition)
library(edgeR)
library(doParallel)
library(lme4)
library(ggpubr)
library(grid)
library(VennDiagram)
library(rlist)
library(gProfileR)
library(biomaRt)
library(GSEABase)
library(WebGestaltR)
library(gridExtra)
library(synapseClient)


#-----------------------------------------------------------------------#


#### FUNCTIONS ####

multiMerge.ensembl = function(df1, df2){
  df = merge(df1, df2, by = "ensembl")
  return(df)
}


#-----------------------------------------------------------------------#


#### LOADING DATA ####

# load result list and contrasts.fit object from DE analysis #
# res.list #
# fit.cont #


#-----------------------------------------------------------------------#


#### SYNERGY CATEGORIES 1 ####

# calculate log2FC standard errors and mean
SE = sqrt(fit.cont$s2.post) * fit.cont$stdev.unscaled
meanSE=mean(SE)

# create table of contrasts' log2FCs
names(res.list)
log2FC.matrix = Reduce(multiMerge.ensembl, list(res.list$AdditiveModel[,c(1,2,4)],
                                                res.list$ALLvsCTRLs[,c(1,4)], 
                                                res.list$SYNERGY[,c(1,4)]))
colnames(log2FC.matrix) = c("ensembl", "Gene_name", "Additive.logFC", "Combinatorial.logFC", "Synergistic.logFC")
rownames(log2FC.matrix) = log2FC.matrix$ensembl

# create columns for synergistic effect categories
log2FC.matrix$syn.sign = NA         #synergistic effect sign
log2FC.matrix$magnitude.diff = NA   #synergistic effect magnitude

# for all genes ...
for (i in 1:length(log2FC.matrix$Gene_name)){
  
  # assign negative and positive synergy, based on Synergy logFC sign
  if (log2FC.matrix$Synergistic.logFC[i] > meanSE){
    log2FC.matrix$syn.sign[i] = "pos"
  } 
  else if (log2FC.matrix$Synergistic.logFC[i] < -meanSE){
    log2FC.matrix$syn.sign[i] = "neg"
  } 
  else log2FC.matrix$syn.sign[i] = "same"
  
  # assign magnitude expression, based on additive model logFC sign and dependent on Synergy logFC sign
  if (log2FC.matrix$syn.sign[i] == "same"){
    log2FC.matrix$magnitude.diff[i] = "same"
  } 
  else if (log2FC.matrix$syn.sign[i] == "neg"){
    if (log2FC.matrix$Additive.logFC[i] > meanSE){
      log2FC.matrix$magnitude.diff[i] = "less"
    } else log2FC.matrix$magnitude.diff[i] = "more"
  }
  else if (log2FC.matrix$syn.sign[i] == "pos"){
    if (log2FC.matrix$Additive.logFC[i] < -meanSE){
      log2FC.matrix$magnitude.diff[i] = "less"
    } else log2FC.matrix$magnitude.diff[i] = "more"
  }
  
  # assign "opposite" magnitude, if additive and combinatorial models have opposing logFC signs
  if (log2FC.matrix$syn.sign[i] == "neg" && log2FC.matrix$Additive.logFC[i] > meanSE && log2FC.matrix$Combinatorial.logFC[i] < -meanSE){
    log2FC.matrix$magnitude.diff[i] = "opposite"
  }
  if (log2FC.matrix$syn.sign[i] == "pos" && log2FC.matrix$Additive.logFC[i] < -meanSE && log2FC.matrix$Combinatorial.logFC[i] > meanSE){
    log2FC.matrix$magnitude.diff[i] = "opposite"
  }
  
}

#### SYNERGY CATEGORIES 2 ####

# copy magnitude column
log2FC.matrix$magnitude.diff2 = log2FC.matrix$magnitude.diff

for (i in 1:length(log2FC.matrix$Gene_name)){
  
  # assign more specific magnitude expressions
  if (log2FC.matrix$syn.sign[i] == "same"){
    if (log2FC.matrix$Additive.logFC[i] > meanSE){
      log2FC.matrix$magnitude.diff2[i] = "same.up"
    } 
    else if (log2FC.matrix$Additive.logFC[i] < -meanSE){
      log2FC.matrix$magnitude.diff2[i] = "same.down"
    }
    else log2FC.matrix$magnitude.diff2[i] = "same"
  } 
  else if (log2FC.matrix$syn.sign[i] == "neg"){
    if (log2FC.matrix$Additive.logFC[i] > meanSE){
      log2FC.matrix$magnitude.diff2[i] = "less.up"
    } else log2FC.matrix$magnitude.diff2[i] = "more.down"
  }
  else if (log2FC.matrix$syn.sign[i] == "pos"){
    if (log2FC.matrix$Additive.logFC[i] < -meanSE){
      log2FC.matrix$magnitude.diff2[i] = "less.down"
    } else log2FC.matrix$magnitude.diff2[i] = "more.up"
  }
  if (log2FC.matrix$syn.sign[i] == "neg" && log2FC.matrix$Additive.logFC[i] > meanSE && log2FC.matrix$Combinatorial.logFC[i] < -meanSE){
    log2FC.matrix$magnitude.diff2[i] = "opposite.down"
  }
  if (log2FC.matrix$syn.sign[i] == "pos" && log2FC.matrix$Additive.logFC[i] < -meanSE && log2FC.matrix$Combinatorial.logFC[i] > meanSE){
    log2FC.matrix$magnitude.diff2[i] = "opposite.up"
  }
  
}

#### SYNERGY CATEGORIES 3 ####

# simplify opposite expression in yet another column
log2FC.matrix$magnitude.diff3 = log2FC.matrix$magnitude.diff2
for (i in 1:length(log2FC.matrix$Gene_name)){
  if (log2FC.matrix$magnitude.diff2[i] == "opposite.down"){
    log2FC.matrix$magnitude.diff3[i] = "more.down"
  }
  if (log2FC.matrix$magnitude.diff2[i] == "opposite.up"){
    log2FC.matrix$magnitude.diff3[i] = "more.up"
  }
}