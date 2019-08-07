# Nadine Schrode
# 07/29/19
# input: meta data matrix, counts matrix, gene annotation table


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
mds = function(normDGE, pos, metacol){
  col = rainbow(length(levels(metacol)), 1, 0.8, alpha = 0.5)[metacol]
  plotMDS(normDGE, col = col, pch = 16, cex = 2)
  if (pos == "top"){
    legend("top", fill = rainbow(length(levels(metacol)), 1, 0.8), legend = levels(metacol), xpd = T, horiz = T, inset = c(-0.1, -0.1), bty = "n")
  }else{
    legend(pos, fill = rainbow(length(levels(metacol)), 1, 0.8), legend = levels(metacol), bty = "n")
  }
}


#-----------------------------------------------------------------------#


#### LOADING DATA, CLEAN-UP and DGE transformation ####

# load meta data, count matrices and gene annotation table #
# counts #
# meta #
# anno #

# sort meta data to match counts sample order
meta = meta[match(colnames(counts), meta$samples), ]

# remove lowly expressed genes
keep = rowSums(cpm(counts) > 0.75) >= 2
gExpr = counts[keep, ]
dim(gExpr)

# DGE transformation
y = DGEList(gExpr)
y = calcNormFactors(y)
y$samples$group = meta$mod.gene
y$samples$line = meta$line
y$samples$donor = meta$donor
# add gene annotation
anno = anno[match(rownames(y), anno$ensembl), ]
y$genes = subset(anno, select = c("ensembl", "Gene_name"))


#-----------------------------------------------------------------------#


#### MDS plots ####
pdf("mds.pdf")
plotMDS(y)
mds(y, "top", factor(meta$donor))
mds(y, "top", factor(meta$VPR))
mds(y, "top", factor(meta$shCtrl))
mds(y, "top", factor(meta$sgCtrl))
mds(y, "top", factor(meta$sgTSNARE1))
mds(y, "top", factor(meta$sgSNAP91))
mds(y, "top", factor(meta$sgCLCN3))
mds(y, "top", factor(meta$line))
mds(y, "top", factor(meta$mod))
mds(y, "top", factor(meta$mod.gene))
dev.off()


#-----------------------------------------------------------------------#


#### DIFFERENTIAL EXPRESSION ####
# voom transform and fit linear model
design = model.matrix(~ 0 + mod.gene + line + donor, meta)
colnames(design) = gsub("mod.gene", "", colnames(design))
v = voom(y, design, plot = TRUE)
fit = lmFit(v, design)

# contrasts
cont.matrix = makeContrasts(TSNAREvsCTRL = tsnare1 - ctrl,
                            SNAP91vsCTRL = sanp91 - ctrl,
                            CLCN3vsCTRL = clcn3 - ctrl, 
                            FURINvsCTRL = furin - ctrl,
                            AdditiveModel = tsnare1 + sanp91 + clcn3 + furin - 4*ctrl,
                            ALLvsCTRLs = all - all.ctrl,
                            SYNERGY = all - tsnare1 - sanp91 - clcn3 - furin - all.ctrl + 4*ctrl,
                            levels = design)
# calculate coefficients
fit.cont = contrasts.fit(fit, cont.matrix)
# compute moderated t-statistics, and log-odds by empirical Bayes moderation
fit.cont = eBayes(fit.cont)
summa.fit = decideTests(fit.cont, adjust.method = "fdr")
summary(summa.fit)

# create and save contrast matrix heatmap
cont.p = t(cont.matrix[c(1:7), ])
h = pheatmap(cont.p,
             display_numbers = T, number_format = "%.0f",
             breaks = seq(-3, 1, by = 0.5),
             color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(12),
             cluster_cols = F, cluster_rows = F)
pdf("contrast-heatmap.pdf")
print(h)
dev.off()


#-----------------------------------------------------------------------#


# create variables for each contrast and save in list object for further analysis
res.list = list()
for (i in 1:length(colnames(fit.cont$contrasts))){
  x = topTable(fit.cont, coef = i, sort.by = "p", n = Inf)
  res.list[[i]] = x
  names(res.list)[i] = colnames(fit.cont$contrasts)[i]
  write.csv(x, paste0(colnames(fit.cont$contrasts)[i], ".csv"))
}
list2env(res.list, globalenv())


#-----------------------------------------------------------------------#


#### VOLCANO and MD plots ####
# create and save plots for all contrasts
pdf("md-volcano.pdf")
par(mfrow = c(1, 2))
for (i in 1:length(colnames(fit.cont$contrasts))){
  plotMD(fit.cont, coef = i, status = summa.fit[, i], values = c(-1, 1))
  volcanoplot(fit.cont, coef = i, highlight = 10, names = fit.cont$genes$Gene_name)
} 
dev.off()


#-----------------------------------------------------------------------#
# res.list object will be used in further analysis
#-----------------------------------------------------------------------#
