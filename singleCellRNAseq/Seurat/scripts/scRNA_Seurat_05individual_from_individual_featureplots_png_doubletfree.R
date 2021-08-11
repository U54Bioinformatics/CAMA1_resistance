library(Seurat)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)


cell_number_png <- function(cells, filename, figdir, size=4, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    y <- summary(cells$Sample)
    df <- data.frame(Sample=names(y), Cell_numbers=y)
    df <- df[df$Cell_numbers!=0,]
    df 
    p <- ggplot(data=df, aes(x = Sample, y = Cell_numbers, width=.5)) + labs(title="Cell numbers", x="", y="Cell numbers") + geom_bar(stat="identity") + geom_text(aes(x = Sample, y = Cell_numbers*1.1, label = Cell_numbers), size= font_size/2) + theme(plot.title=element_text(size=font_size,face="bold", hjust = 0.5), axis.text=element_text(size=font_size, face="bold"), axis.text.x = element_text(angle = 45, hjust = 1), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    p

}

dimplot_png <- function(cells, method, filename, figdir, group, size=4, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")    
    png(filename, width = width, height = height)
    font_size = 46
    DimPlot(cells, reduction = method, group.by = group, pt.size=size, label.size=12, label=TRUE) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    #dev.off()
}

featureplot_png <- function(cells, features, filename, figdir, size=4, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    FeaturePlot(cells, features = features) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    #dev.off()
}

volinplot_png <- function(cells, features, group, filename, figdir, size=2, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    VlnPlot(cells, features = features, group.by = group) + labs(title=features, x="", y=features) +  theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    #dev.off()
}

doheatmap_png <- function(cells, features, filename, figdir, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 26
    DoHeatmap(pbmc, features = features) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    #dev.off()
}


args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
platform=args[2]
#sample="FEL011"
#platform="10x"
prefix=paste(sample, platform, sep="_")
figdir=paste(sample, platform, "seurat_doubletfree_png", sep="_")

#pdf(paste(prefix, "_Seurat_2kgenes_vst_cc.featureplots.pdf", sep=""), height=7, width=10)

#load cell types for UMAP clusters from singleR
pbmc <- readRDS(paste(prefix, "_Seurat_2kgenes_vst_cc.doublet_free.rds", sep=""))
#celltype <- fread(paste(prefix, "_cell_metadata.UMAPcluster.SingleR_anno_cluster.txt", sep=""))
#rownames(celltype) <- celltype$Cell.ID
#pbmc <- AddMetaData(pbmc, celltype)
#pbmc <- AddMetaData(pbmc, celltype$Encode_type, col.name="Encode_type")
#pbmc <- AddMetaData(pbmc, celltype$Encode_main_type, col.name="Encode_main_type")
#pbmc <- AddMetaData(pbmc, celltype$hpca_type, col.name="hpca_type")
#pbmc <- AddMetaData(pbmc, celltype$hpca_main_type, col.name="hpca_main_type")
#pbmc <- AddMetaData(pbmc, celltype$Encode_type_cluster, col.name="Encode_type_cluster")
#pbmc <- AddMetaData(pbmc, celltype$seurat_clusters2, col.name="seurat_clusters2")
#Encode_type     Encode_main_type        hpca_type       hpca_main_type

#Load exp and wellList
#counts.df <-  read.table(paste(sample, platform, "gene_symbols.filtered.counts.txt", sep="_"), sep="\t")
#counts.df <- as.data.frame(counts.df)
##counts.df <- subset(counts.df, counts.df$'Gene Symbol' != "")
##counts.df <- unique(counts.df, by='Gene Symbol')
##counts.df <- as.data.frame(counts.df)
##rownames(counts.df) <- counts.df$`Gene Symbol`
##counts <- counts.df[,4:ncol(counts.df)]
##write.table(counts, paste(prefix, "gene_symbols.counts.txt", sep="_"), sep="\t", quote = F)
##Annotation
#Annot.df <- read.table(paste(sample, platform, "cell_metadata.UMAPcluster.singleR_cell_cluster_anno.txt", sep="_"), header=TRUE, sep = "\t")
#rownames(Annot.df) <- Annot.df$Cell.ID
#Annot.df <- Annot.df[,2:ncol(Annot.df)]
#Annot.df  <- Annot.df[match(colnames(counts.df), rownames(Annot.df)),]
##Create Seurat Object
#pbmc <- CreateSeuratObject(counts.df,  project = sample,  meta.data = Annot.df,  min.cells = 3, min.features = 300)

#
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
#DimPlot(pbmc, reduction = "pca", group.by = "Phase")
#DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(pbmc)


#pbmc <- FindNeighbors(pbmc, dims = 1:20)
#pbmc <- FindClusters(pbmc, resolution = 0.5)

#pbmc <- RunTSNE(pbmc, dims = 1:20)
#pbmc <- RunUMAP(pbmc, dims = 1:20)


#plot features on top of cluster and cells
#FeaturePlot(object = pbmc, features = c("nCount_RNA"), pt.size = 4)
#FeaturePlot(object = pbmc, features = c("nFeature_RNA"), pt.size = 4)

#pdf plots
#plot TSNE
#DimPlot(pbmc, reduction = "tsne", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Sample", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Encode_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Encode_main_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "hpca_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "hpca_main_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Encode_type_cluster", label=TRUE, pt.size = 1)
#+ NoLegend() to remove legend
#plot UMAP
#DimPlot(pbmc, reduction = "umap", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "Sample", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "Encode_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "Encode_main_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "hpca_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "hpca_main_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "umap", group.by = "Encode_type_cluster", label=TRUE, pt.size = 1)

#FeaturePlot(object = pbmc, features = c("nCount_RNA"), pt.size = 4)
#FeaturePlot(object = pbmc, features = c("nFeature_RNA"), pt.size = 4)


fig_cell_num = paste(sample, platform, "Cell_number", sep="_")
cell_number_png(pbmc, fig_cell_num, figdir)

fig_rna = paste(sample, platform, "nCount_RNA", sep="_")
volinplot_png(pbmc, c("nCount_RNA"), "Sample", fig_rna, figdir)
fig_feature = paste(sample, platform, "nFeature_RNA", sep="_")
volinplot_png(pbmc, c("nFeature_RNA"), "Sample", fig_feature, figdir)
fig_mt = paste(sample, platform, "percent.mt", sep="_")
volinplot_png(pbmc, c("percent.mt"), "Sample", fig_mt, figdir)


#fig_tsne = paste(sample, platform, 'umap', "ident", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "ident")
#fig_tsne = paste(sample, platform, 'umap', "Sample", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Sample")
#fig_tsne = paste(sample, platform, 'umap', "Encode_main_type", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Encode_main_type")
#fig_tsne = paste(sample, platform, 'umap', "hpca_main_type", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "hpca_main_type")
#fig_tsne = paste(sample, platform, 'umap', "Encode_type_cluster", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Encode_main_cluster")

#fig_tsne = paste(sample, platform, 'tsne', "ident", sep="_")
#dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "ident")
#fig_tsne = paste(sample, platform, 'tsne', "Sample", sep="_")
#dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Sample")
#fig_tsne = paste(sample, platform, 'tsne', "Encode_main_type", sep="_")
#dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Encode_main_type")
#fig_tsne = paste(sample, platform, 'tsne', "hpca_main_type", sep="_")
#dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "hpca_main_type")
#fig_tsne = paste(sample, platform, 'tsne', "Encode_type_cluster", sep="_")
#dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Encode_main_cluster")


#pdf plots
#Immune cell
#FeaturePlot(pbmc, features = c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A"))
#VlnPlot(pbmc, features = c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A"))
#Fibroblast cell
#FeaturePlot(pbmc, features = c("PDPN", "THY1", "CD34", "CDH11"))
#VlnPlot(pbmc, features = c("PDPN", "THY1", "CD34", "CDH11"))
#epithelial cell
#FeaturePlot(pbmc, features = c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"))
#VlnPlot(pbmc, features = c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"))

#png plots
#Immune cell
#fig_immune_feature = paste(sample, platform, "Immunecell_feature", sep="_")
#featureplot_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A"), fig_immune_feature, figdir)
#fig_immune_violin = paste(sample, platform, "Immunecell_violin", sep="_")
#volinplot_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A"), "seurat_clusters", fig_immune_violin, figdir)
#Fibroblast cell
#fig_fibroblast_feature = paste(sample, platform, "Fibroblast_feature", sep="_")
#featureplot_png(pbmc, c("PDPN", "THY1", "CD34", "CDH11"), fig_fibroblast_feature, figdir)
#fig_fibroblast_violin = paste(sample, platform, "Fibroblast_violin", sep="_")
#volinplot_png(pbmc, c("PDPN", "THY1", "CD34", "CDH11"), "seurat_clusters", fig_fibroblast_violin, figdir)
#epithelial cell
#fig_epithelial_feature = paste(sample, platform, "Epithelial_feature", sep="_")
#featureplot_png(pbmc, c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"), fig_epithelial_feature, figdir)
#fig_epithelial_violin = paste(sample, platform, "Epithelial_violin", sep="_")
#volinplot_png(pbmc, c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"), "seurat_clusters", fig_epithelial_violin, figdir)

#find marker genes
#cells.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#pdf
#DoHeatmap(object = pbmc, features = top10$gene)
#png
#fig_markers = paste(sample, platform, "top10_genes", sep="_")
#doheatmap_png(pbmc, top10$gene, fig_markers, figdir)

###Write metadata
#out <- data.table(Cell.ID = colnames(x=pbmc))
#out_meta <- cbind(out, pbmc@meta.data)
#add one columen with cluster# name (numbers are not working with "merge function")
#paste("Cluster", Annot.df$seurat_clusters, sep="")
#cluster <- as.data.frame(paste("Cluster", out_meta$seurat_clusters, sep=""))
#colnames(cluster) <- "seurat_clusters2"
#cluster$Cell.ID <- out_meta$Cell.ID
#out_meta <- merge(out_meta, cluster, by="Cell.ID")
#fwrite(out_meta, paste(prefix, "_cell_metadata.UMAPcluster.SingleR_anno_cluster.featureplots_confirm.txt", sep=""), sep="\t", quote=F, col.names=T)

#saveRDS(pbmc, file = paste(prefix, "_Seurat_2kgenes_vst_cc.featureplots.rds", sep=""))

#dev.off()

