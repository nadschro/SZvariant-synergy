# Nadine Schrode
# 08/20/19
# code adapted from Gabriel Hoffman: https://github.com/GabrielHoffman/COS_public_release
# input: Differential expression results list from DE-limma script and DE tables of different disorders


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


#### LOADING DATA ####

# load result list from DE analysis #
# res.list #

# load disorder DE tables #
# resDE_ETOH #
# resDE_MDD #
# resDE_BD #
# resDE_SCZ #
# resDE_ASD #


#-----------------------------------------------------------------------#


#### CORRELATION ####

# create vector to store final sample names 
# (retaining order the samples have in res.list)
final.names = c("TSNARE CRISPRa", 
                "SNAP91 CRISPRa", 
                "CLCN3 CRISPRa", 
                "FURIN RNAi", 
                "Expected Additive Model", 
                "Measured Combinatorial Perturbation", 
                "Synergistic effect")

# compute correlation stats between all current and disorder data 
dis.table = data.frame(matrix(ncol = 5, nrow = 0))
for (i in 1:length(res.list)){
  concordanceDF = matrix(NA, ncol = 5, nrow = 0)
  resDE_CRISPR = res.list[[i]]
  name = final.names[i]
  # ETOH
  # logFC
  df = merge(resDE_ETOH, resDE_CRISPR,by.x = "genes",by.y = "ensembl")  
  corValues = with(df, cor.test(logFC.x, logFC.y, method = "spearman", exact = F))
  concordanceDF = rbind(concordanceDF, c("ETOH", name, "logFC", corValues$estimate,  corValues$p.value))
  # t-statistics
  df = merge(resDE_ETOH, resDE_CRISPR,by.x = "genes",by.y = "ensembl")
  corValues = with(df, cor.test(t.x, t.y, method = "spearman", exact = F))
  concordanceDF = rbind(concordanceDF, c("ETOH", name, "t", corValues$estimate,  corValues$p.value))
  # MDD
  # logFC
  df = merge(resDE_MDD, resDE_CRISPR,by.x = "genes",by.y = "ensembl")  
  corValues = with(df, cor.test(logFC.x, logFC.y, method = "spearman", exact = F))
  concordanceDF = rbind(concordanceDF, c("MDD", name, "logFC", corValues$estimate,  corValues$p.value))
  # t-statistics
  df = merge(resDE_MDD, resDE_CRISPR,by.x = "genes",by.y = "ensembl")
  corValues = with(df, cor.test(t.x, t.y, method = "spearman", exact = F))
  concordanceDF = rbind(concordanceDF, c("MDD", name, "t", corValues$estimate,  corValues$p.value))
  # BD
  # logFC
  df = merge(resDE_BD, resDE_CRISPR,by.x = "genes",by.y = "ensembl")  
  corValues = with(df, cor.test(logFC.x, logFC.y, method = "spearman", exact = F))
  concordanceDF = rbind(concordanceDF, c("BD", name, "logFC", corValues$estimate,  corValues$p.value))
  # t-statistics
  df = merge(resDE_BD, resDE_CRISPR,by.x = "genes",by.y = "ensembl")
  corValues = with(df, cor.test(t.x, t.y, method = "spearman", exact = F))
  concordanceDF = rbind(concordanceDF, c("BD", name, "t", corValues$estimate,  corValues$p.value))
  # SCZ
  # logFC
  df = merge(resDE_SCZ, resDE_CRISPR,by.x = "genes",by.y = "ensembl")  
  corValues = with(df, cor.test(logFC.x, logFC.y, method = "spearman", exact = F))
  concordanceDF = rbind(concordanceDF, c("SCZ", name, "logFC", corValues$estimate,  corValues$p.value))
  # t-statistics
  df = merge(resDE_SCZ, resDE_CRISPR,by.x = "genes",by.y = "ensembl")
  corValues = with(df, cor.test(t.x, t.y, method = "spearman", exact = F))
  concordanceDF = rbind(concordanceDF, c("SCZ", name, "t", corValues$estimate,  corValues$p.value))
  # ASD (scatterplots)
  # logFC
  df = merge(resDE_ASD, resDE_CRISPR,by.x = "genes",by.y = "ensembl")  
  corValues = with(df, cor.test(logFC.x, logFC.y, method = "spearman", exact = F))
  concordanceDF = rbind(concordanceDF, c("ASD", name, "logFC", corValues$estimate,  corValues$p.value))
  # t-statistics
  df = merge(resDE_ASD, resDE_CRISPR,by.x = "genes",by.y = "ensembl")
  corValues = with(df, cor.test(t.x, t.y, method = "spearman", exact = F))
  concordanceDF = rbind(concordanceDF, c("ASD", name, "t", corValues$estimate,  corValues$p.value))
  
  # prep data frame and save
  concordanceDF = data.frame(concordanceDF, stringsAsFactors = FALSE)
  colnames(concordanceDF) = c("dataset1", "dataset2", "stat", "rho", "pvalue")
  concordanceDF$rho = as.numeric(concordanceDF$rho)
  concordanceDF$pvalue = as.numeric(concordanceDF$pvalue)
  #write.csv(concordanceDF,file = paste0(name,"Disorders_table.csv"))
  colnames(dis.table) = colnames(concordanceDF)
  dis.table = rbind(dis.table, concordanceDF)
  write.csv(dis.table, "Disorders_summary_table.csv")
}


#-----------------------------------------------------------------------#


#### PLOT ####

# merged disorders barplot 
dis.table$dataset2 = factor(dis.table$dataset2,
                      levels = rev(c("SNAP91 CRISPRa", "TSNARE CRISPRa",
                                     "FURIN RNAi", "CLCN3 CRISPRa",
                                     "Measured Combinatorial Perturbation",
                                     "Synergistic effect")))
dis.table$dataset1 = factor(dis.table$dataset1, 
                      levels = c("SCZ","BD","MDD","ASD","ETOH"))

g = ggplot(subset(dis.table, stat == statistic),
           aes(dataset2, rho, fill = dataset2)) +
    geom_bar(stat = 'identity', position = position_dodge()) +  
    xlab('') + 
    ylab(bquote(Correlation~(Spearman~rho))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1,legend.position = "none") +
    coord_flip() + 
    facet_grid(~dataset1)

pdf("CRISPR.data.sets.in.disorders.pdf")
print(g)
dev.off()


#-----------------------------------------------------------------------#

