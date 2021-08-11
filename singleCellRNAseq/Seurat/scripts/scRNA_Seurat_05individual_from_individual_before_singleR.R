library(Seurat)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)


dimplot_png <- function(cells, method, filename, figdir, group, size=8, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")    
    png(filename, width = width, height = height)
    font_size_text = 16
    font_size_axis = 42
    DimPlot(cells, reduction = method, group.by = group, pt.size=size, label.size=font_size_text, label=TRUE) + theme(plot.title=element_text(size=font_size_axis, face="bold"), axis.text=element_text(size=font_size_axis, face="bold"), axis.title=element_text(size=font_size_axis, face="bold"), legend.text=element_text(size=font_size_axis, face="bold"), legend.title=element_text(size=font_size_axis, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = font_size_text)))
}

featureplot_png <- function(cells, features, filename, figdir, size=12, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size_text = 26
    font_size_axis = 42
    FeaturePlot(cells, features = features) + theme(plot.title=element_text(size=font_size_axis, face="bold"), axis.text=element_text(size=font_size_axis, face="bold"), axis.title=element_text(size=font_size_axis, face="bold"), legend.text=element_text(size=font_size_axis, face="bold"), legend.title=element_text(size=font_size_axis, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape
 = guide_legend(override.aes = list(size = 25)))
    #dev.off()
}

volinplot_png <- function(cells, features, filename, figdir, size=12, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 42
    VlnPlot(cells, features = features) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 12)))
    #dev.off()
}

doheatmap_png <- function(cells, features, filename, figdir, width=1200, height=1800){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 26
    DoHeatmap(pbmc, features = features) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    #dev.off()
}


args <- commandArgs(trailingOnly = TRUE)
#sample=args[1]
#platform=args[2]
#patient=args[3]
#sample="FEL011"
#platform="10x"
#prefix=paste(sample, platform, sep="_")
prefix=args[1]
figdir=paste(prefix, "seurat_figures_png", sep="_")

pdf(paste(prefix, "Seurat_2kgenes_vst_cc.pdf", sep="_"), height=7, width=10)

#load cell types for UMAP clusters from singleR
pbmc <- readRDS(paste(prefix, "Seurat_2kgenes_vst_cc.raw.rds", sep="_"))
#celltype <- fread(paste(prefix, "_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt", sep=""))
#rownames(celltype) <- celltype$Cell.ID
#pbmc <- AddMetaData(pbmc, celltype)
#
pbmc <- subset(pbmc, subset = nFeature_RNA < 7000 & nFeature_RNA > 3000 & nCount_RNA > 2000 & nCount_RNA < 60000 & percent.mt < 30)
#pbmc <- subset(pbmc, subset = orig.ident == patient)
#
##############################################################################################
###NormalizeData and ScaleData
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
#scaling
#pbmc <- ScaleData(pbmc)
#Regress during scaling
#skip regression first to make cell markers are present then do regression to generate results
if ( TRUE ) {
   pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("percent.mt", "nCount_RNA","S.Score", "G2M.Score"))
   #pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("percent.mt"))
}
##############################################################################################

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca", group.by = "Phase")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = c(2.0, 1.2, 0.5, 0.8))

pbmc <- RunTSNE(pbmc, dims = 1:20)
pbmc <- RunUMAP(pbmc, dims = 1:20)

DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.0.5", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.0.8", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.1.2", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.2", label=TRUE, pt.size = 1)
pbmc <- BuildClusterTree(pbmc)
PlotClusterTree(pbmc)

#plot features on top of cluster and cells
FeaturePlot(object = pbmc, features = c("nCount_RNA"), pt.size = 4)
FeaturePlot(object = pbmc, features = c("nFeature_RNA"), pt.size = 4)

#pdf plots
#plot TSNE
DimPlot(pbmc, reduction = "tsne", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "tsne", group.by = "Sample", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Encode_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Encode_main_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "hpca_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "hpca_main_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Encode_main_cluster", label=TRUE, pt.size = 1)
#+ NoLegend() to remove legend
#plot UMAP
DimPlot(pbmc, reduction = "umap", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "umap", group.by = "Sample", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "Encode_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "Encode_main_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "hpca_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "hpca_main_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "Encode_main_cluster", label=TRUE, pt.size = 1)

#pdf plots
#Immune cell
#FeaturePlot(pbmc, features = c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"))
#VlnPlot(pbmc, features = c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"))
#Fibroblast cell
#FeaturePlot(pbmc, features = c("PDPN", "THY1", "CD34", "CDH11"))
#VlnPlot(pbmc, features = c("PDPN", "THY1", "CD34", "CDH11"))
#epithelial cell
#FeaturePlot(pbmc, features = c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"))
#VlnPlot(pbmc, features = c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"))

dev.off()

#png plots
fig_tsne = paste(prefix, 'umap', "ident", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "ident")
fig_tsne = paste(prefix, 'umap', "Sample", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Sample")
#fig_tsne = paste(sample, platform, patient, 'umap', "Encode_main_type", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Encode_main_type")
#fig_tsne = paste(sample, platform, patient, 'umap', "hpca_main_type", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "hpca_main_type")
#fig_tsne = paste(sample, platform, patient, 'umap', "Encode_type_cluster", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Encode_main_cluster")

fig_tsne = paste(prefix, 'tsne', "ident", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "ident")
fig_tsne = paste(prefix, 'tsne', "Sample", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Sample")
#fig_tsne = paste(sample, platform, patient, 'tsne', "Encode_main_type", sep="_")
#dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Encode_main_type")
#fig_tsne = paste(sample, platform, patient, 'tsne', "hpca_main_type", sep="_")
#dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "hpca_main_type")
#fig_tsne = paste(sample, platform, patient, 'tsne', "Encode_type_cluster", sep="_")
#dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Encode_main_cluster")

#png plots
#Immune cell
#fig_immune_feature = paste(prefix, "Immunecell_feature", sep="_")
#featureplot_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"), fig_immune_feature, figdir)
#fig_immune_violin = paste(prefix, "Immunecell_violin", sep="_")
#volinplot_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"), fig_immune_violin, figdir)
#Fibroblast cell
#fig_fibroblast_feature = paste(prefix, "Fibroblast_feature", sep="_")
#featureplot_png(pbmc, c("PDPN", "THY1", "CD34", "CDH11"), fig_fibroblast_feature, figdir)
#fig_fibroblast_violin = paste(prefix, "Fibroblast_violin", sep="_")
#volinplot_png(pbmc, c("PDPN", "THY1", "CD34", "CDH11"), fig_fibroblast_violin, figdir)
#epithelial cell
#fig_epithelial_feature = paste(prefix, "Epithelial_feature", sep="_")
#featureplot_png(pbmc, c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"), fig_epithelial_feature, figdir)
#fig_epithelial_violin = paste(prefix, "Epithelial_violin", sep="_")
#volinplot_png(pbmc, c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"), fig_epithelial_violin, figdir)
#TEM
#fig_epithelial_feature = paste(prefix, "MacrophagesTEM_feature", sep="_")
#featureplot_png(pbmc, c("CD200R1", "TNFSF10", "HGF", "ANGPT1"), fig_epithelial_feature, figdir)
#fig_epithelial_violin = paste(prefix, "MacrophagesTEM_violin", sep="_")
#volinplot_png(pbmc, c("CD200R1", "TNFSF10", "HGF", "ANGPT1"), fig_epithelial_violin, figdir)


#find marker genes
#cells.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#pdf
#DoHeatmap(object = pbmc, features = top10$gene)
#png
#fig_markers = paste(prefix, "top10_genes", sep="_")
#doheatmap_png(pbmc, top10$gene, fig_markers, figdir)


###Write metadata
out <- data.table(Cell.ID = colnames(x=pbmc))
out_meta <- cbind(out, pbmc@meta.data)
#add one columen with cluster# name (numbers are not working with "merge function")
#paste("Cluster", Annot.df$seurat_clusters, sep="")
#cluster <- as.data.frame(paste("Cluster", out_meta$seurat_clusters, sep=""))
#colnames(cluster) <- "seurat_clusters2"
#cluster$Cell.ID <- out_meta$Cell.ID
#out_meta <- merge(out_meta, cluster, by="Cell.ID")
fwrite(out_meta, paste(prefix, "cell_metadata.UMAPcluster.txt", sep="_"), sep="\t", quote=F, col.names=T)

###Write count
counts <- GetAssayData(object = pbmc, slot = "counts")
counts <- as.matrix(counts)
write.table(counts, paste(prefix, "gene_symbols.filtered.counts.txt", sep="_"), sep="\t", quote = F)


saveRDS(pbmc, file = paste(prefix, "Seurat_2kgenes_vst_cc.rds", sep="_"))
#dev.off()

