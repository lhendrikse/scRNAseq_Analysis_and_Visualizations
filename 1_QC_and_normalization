# Info --------------------------------------------------------------------
# Liam Hendrikse
# May 14, 2020
# Taylor lab

# Notes -------------------------------------------------------------------
#

# Libraries ---------------------------------------------------------------
library(dplyr)
library(Seurat)
library(reticulate)
library(harmony)
library(ggplot2)
library(Matrix)
library(future)
library(glue)
library(tools)

# Inputs ------------------------------------------------------------------
options(future.globals.maxSize = 10000 * 1024^2, stringsAsFactors = F)
plan("multiprocess", workers = 8)
args = commandArgs(trailingOnly=TRUE)
cr_path = args[1] #path to output from cellranger
id = args[2] #sample IDs

# Directories -------------------------------------------------------------
out_dir = "/outs/"
plot_dir = "/outs/plots/"
dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)
setwd(glue("{cr_path}/outs/filtered_feature_bc_matrix/"))

# Imports -----------------------------------------------------------------
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes

so = Read10X(data.dir = ".", gene.column = 2, unique.features = TRUE)

so.list = list(so)

# Create Seurat Objects ---------------------------------------------------
for (i in 1:length(so.list)) {
        # Initialize Seurat object 
        so.list[[i]] = CreateSeuratObject(counts = so.list[[i]], project = id[i], min.cells = 10, min.features = 500)
        so.list[[i]][["percent.mt"]]  <- PercentageFeatureSet(so.list[[i]], pattern = "^MT-") # pattern = "^mt-" for mice
        
        # Pre-QC diagnostic plots
        pdf(paste0(plot_dir, id, "_Pre_QC_Violins.pdf")
        print(VlnPlot(so.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
        dev.off()
        
        # QC cutoffs per sample
        topfeature = (quantile(so.list[[i]]$nFeature_RNA,.75)+(1.5*IQR(so.list[[i]]$nFeature_RNA)))
        lowfeature = (quantile(so.list[[i]]$nFeature_RNA,.25)-(1.5*IQR(so.list[[i]]$nFeature_RNA)))
        topmito =(quantile(so.list[[i]]$percent.mt,.75)+(1.5*IQR(so.list[[i]]$percent.mt)))
        lowmito = (quantile(so.list[[i]]$percent.mt,.25)-(1.5*IQR(so.list[[i]]$percent.mt)))
        
        # Remove low quality cells
        so.list[[i]] = subset(so.list[[i]], subset = nFeature_RNA < topfeature & nFeature_RNA > lowfeature & percent.mt < topmito & percent.mt > lowmito)
        
        # Perform scoring of cell cycle genes in order to do linear regression during normalization steps
        so.list[[i]] = CellCycleScoring(so.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
        saveRDS(so.list[[i]], file = paste0(out_dir, id[i], "_scRNAseq_raw.rds"))

        # Post-QC diagnostic plots
        pdf(paste0(plot_dir, id, "_Post_QC_Violins.pdf")
        print(VlnPlot(so.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
        dev.off()
}

# Read raw seurat objects for normalization after checking QC
for (i in 1:length(id)) {
	assign(id[i],readRDS(paste0(out_dir, id[i], "_scRNAseq_raw.rds")))
}

# Run SCTransform normalization -------------------------------------------

SCT.list = list(mget(ls(pattern = id)))

# Run SCT normalization
for (i in 1:length(SCT.list)) {
	SCT.list[[i]] = SCTransform(SCT.list[[i]], vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),variable.features.n=3000,return.only.var.genes=FALSE)
	
	SCT.list[[i]] = RunPCA(SCT.list[[i]],assay="SCT")
	
  # Determine significant PCs
  pdf(paste0(plot_dir, id, "_Elbow.pdf")
	print(ElbowPlot(SCT.list[[i]]))
	dev.off()

	# Jackstraw plots should be assessed in addition to eblow plots for sig PCs
	SCT.list[[i]] = JackStraw(SCT.list[[i]], num.replicate = 100)
	SCT.list[[i]] = ScoreJackStraw(SCT.list[[i]], dims = 1:20)
	
  pdf(paste0(plot_dir, id, "_Jackstraw.pdf")
	print(JackStrawPlot(SCT.list[[i]], dims = 1:20))
	dev.off()

  # Save normalized object
	saveRDS(SCT.list[[i]], file = paste0(out_dir, id[i], "_scRNAseq_sct.rds"))

}

# Run Clustering and dimensionality reduction -----------------------------

# Inspect elbow and jackstraw plots and input
pc_values = list(c(1:20))
res = 0.4

for (i in 1:length(SCT.list)) {
	SCT.list[[i]] = FindNeighbors(SCT.list[[i]],reduction = "pca", dims = pc_values[[i]])
	SCT.list[[i]] = FindClusters(SCT.list[[i]],resolution = res, algorithm=2) 
	SCT.list[[i]] = RunUMAP(SCT.list[[i]],reduction = "pca", dims = pc_values[[i]]) 
  
	saveRDS(SCT.list[[i]], file = paste0(out_dir, id[i], "_scRNAseq_sct.rds"))
  
}

sessionInfo()
rm(list = ls())











