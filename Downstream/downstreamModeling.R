#!/home/bdecato/miniconda3/envs/r4/lib/R/bin/Rscript
#
# downstreamModeling.R
#
# Reproduce some plots from Vince and Jinfeng's analysis and dig further into
# downstream modeling for development of resistance in CAMA-1 cells.
#
# @author Ben Decato

################################################################################
# Libraries
################################################################################
rm(list=ls())

library(Seurat)
library(tidyverse)
options(ggrepel.max.overlaps = Inf)

print(paste("Working directory:", getwd()))

# Set up my ggplot theme for use with all visualizations
myTheme <- theme_bw() +
  theme(text = element_text(size=16), axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.placement = "outside")

################################################################################
# Read and integrate data
################################################################################
#seu <- readRDS("Feline2_integrated_downsampled.rds")
seu <- readRDS("~/Downloads/VG_CAMA1_D11_ALL_10x_Seurat_2kgenes_vst_cc.rds")

DimPlot(seu, reduction = "umap", 
        group.by = "seurat_clusters", 
        label=FALSE, repel=TRUE, pt.size = 0.25, raster = TRUE) + 
  myTheme






















