# Info --------------------------------------------------------------------
# Liam Hendrikse
# Apr 08, 2021
# Taylor lab

# Notes -------------------------------------------------------------------
# Code to create similar plot to Figure 5h and Extended Data Figure 11j of 
# Hendrikse et al., Nature, 2022 https://doi.org/10.1038/s41586-022-05215-w

# Libraries ---------------------------------------------------------------
library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(RColorBrewer)
library(mgcv)
library(gam)
library(pheatmap)
library(tradeSeq)
library(reticulate)
library(destiny)
library(slingshot)
library(ggplot2)
library(Matrix)
library(future)
library(glue)
library(tools)
options(stringsAsFactors = F)

# Functions ---------------------------------------------------------------

'%!in%' = function(x,y)!('%in%'(x,y))

Intersect = function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

scale_rows = function (x) {
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}


# Inputs ------------------------------------------------------------------
options(future.globals.maxSize = 10000 * 1024^2, stringsAsFactors = F)
plan("multiprocess", workers = 8)
args = commandArgs(trailingOnly=TRUE)
so_path = args[1] #path to seurat object
id = args[2] #name of seurat object which has been normalized, clustered, annotated, etc. 

colors_dutch = c('#B53471','#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67','#F79F1F','#A3CB38','#D980FA','#EE5A24','#9980FA','#EA2027','#5758BB')
colors_spanish = c('#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2','#2c2c54','#474787','#aaa69d','#227093','#218c74','#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79','#b33939','#cd6133','#84817a','#cc8e35','#ccae62')
colours = c(colors_dutch, colors_spanish)

# Directories -------------------------------------------------------------
out_dir = "/outs/"
plot_dir = "/outs/plots/"
dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)
setwd(glue("{so_path}/"))

# Imports -----------------------------------------------------------------
so = readRDS(id)
sce = as.SingleCellExperiment(so, assay = "SCT")
sce$seurat_clusters = as.numeric(sce$seurat_clusters)

# Run slingshot -----------------------------------------------------------
sce = slingshot(sce, clusterLabels = 'new_cell_type', reducedDim = 'UMAP', start.clus = "RL-VZ")

# Visualize
colors = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) 
plotcol = colors[cut(sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce)$UMAP, col=plotcol, pch=20, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

# Add results back to seurat object
so$Sling1 = sce$slingPseudotime_1
so$Sling2 = sce$slingPseudotime_2

# Gene of interest expression along pseudotime
gg_data = data.frame(Pseudotime = so$Sling1, 
                     Expression = so@assays$SCT@data[which(rownames(so@assays$SCT@data) == "OTX2"),],
                     Cell_type = so$new_cell_type)
gg_data = subset(gg_data, is.na(gg_data$Pseudotime) == F)
gg_data = gg_data[order(gg_data$Pseudotime),]

ggplot(gg_data , aes(x = Pseudotime, y = Expression, colour = Pseudotime)) + 
  geom_point() + 
  ylab("OTX2 Expression") +
  scale_colour_viridis(option="magma") +
  geom_smooth(method = "lm", formula = (y ~ x + exp(x)), colour = "black") + 
  theme_classic()

# Run tradeSeq ------------------------------------------------------------

sce = fitGAM(sce)
ATres = associationTest(sce, lineages = T)

# Determine signficant lineage genes
ATres_subset_1 = subset(ATres, ATres$pvalue_1 < 0.05)
ATres_subset_2 = subset(ATres, ATres$pvalue_2 < 0.05)
sig_lineage_genes = unique(c(ATres_subset_1, ATres_subset_2))

# Grab pseudotime values and cell type annotations for each cell
plot_data = data.frame(Pseudo = so$Sling1,
                       Cell_type = so$Cell_type)
plot_data = subset(plot_data, is.na(plot_data$Pseudo) == F )

# Density of annotated cell types along pseudotime 
dens = ggplot(data=plot_data, aes(x=Pseudo, group=Cell_type, fill=Cell_type)) +
  geom_density(adjust=2, alpha=0.6) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Pseudotime") +
  ylab("Density") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + NoLegend()

# Differential expression testing to find cell type specific markers
markers = FindAllMarkers(so, test.use = "MAST", only.pos = T)

# Designate thresholds to select top cell type specific markers for each cell type
topgenes = c(subset(markers, markers$cluster == "Early VZ" & markers$pct.2 < 0.3)$gene[1:20],
             subset(markers, markers$cluster == "Prolif. VZ" & markers$pct.2 < 0.3)$gene[1:20]) #... etc 
             
topgenes = unique(topgenes[!is.na(topgenes)]) #In case there are not 20 sig. markers for a cluster
topgenes = Intersect(list(topgenes, sig_lineage_genes)) # Ensure cell type specific marker genes are also pseudotime correlated

# Create smoothened expression heatmap of significant lineage markers
yhatSmooth = predictSmooth(sce, gene = topgenes, nPoints = 50, tidy = F) 
yhatSmooth = yhatSmooth[!duplicated(yhatSmooth),]
heatSmooth = pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE, 
                       cluster_rows = F, 
                       treeheight_row = 0,
                       show_rownames = T,
                       border_color = NA,
                       fontsize_row = 4,
                       #gaps_row = c(10, 24, 35),
                       show_colnames = FALSE)

# Stack density and smoothened expression plots to show significant lineage associated genes correlate with changes in cell type density and differentiation
heatSmooth = as.grob(heatSmooth)
grid.arrange(dens, heatSmooth, nrow = 2, ncol=3,
             widths = c(0.28,3.5, 0.74),
             heights = c(0.75, 3),
             layout_matrix = rbind(c(1, 1, 3),
                                   c(3, 2, 2)))


sessionInfo()
rm(list=ls())




