library(Seurat)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)


dimplot_pdf <- function(pbmc, meta, method, filename, figdir, group, size=4, width=1200, height=1200){
    library(Seurat)
    library(data.table)
    library(dplyr)
    library(RColorBrewer)
    library(ggplot2)

    filename <- paste(figdir, '/', filename, '.ggplot.pdf', sep="")
    pdf(filename, width=8, height=7)
  
    coords <- as.data.frame(Embeddings(object =pbmc, reduction = method))
    celltype <- meta
    rownames(celltype) <- celltype$Cell.ID
    pbmc <- AddMetaData(pbmc, celltype)
    Idents(pbmc) <- pbmc$Encode_main_cluster
    
    head(coords) 
    head(celltype)

    pbmc$UMAP_1 <- coords$UMAP_1
    pbmc$UMAP_2 <- coords$UMAP_2
    #x <- cbind(pbmc$UMAP_1, pbmc$UMAP_2, pbmc$Encode_main_cluster)
    #x <- cbind(pbmc$UMAP_1, pbmc$UMAP_2, pbmc$seurat_clusters)
    x <- cbind(pbmc$UMAP_1,pbmc$UMAP_2,Idents(pbmc))
    colnames(x) <- c("UMAP_1", "UMAP_2", "Cell.Type")
    x <- as.data.frame(x)

    #cell type -> color 
    c <- read.table("cell_type_color.txt", sep="\t", check.names=F,  comment.char = "!", header=T)
    rownames(c) <- c$Cell.Type

    #cell type -> ident number
    cell_type_rank <- cbind(pbmc$Encode_main_cluster, Idents(pbmc))
    rownames(cell_type_rank) <- cell_type_rank[,1]
    cell_type_rank <- unique(cell_type_rank)
    colnames(cell_type_rank) <- c("Cell.Type", "Ident")

     #ident -> color 
    color_data <- merge(c, cell_type_rank)
    color_data <- color_data[order(color_data$Ident),]
    color_flag <- as.vector(color_data$Color)
    color_label<- as.vector(color_data$Cell.Type)
    #color_flag <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")
    #color_label <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")

    #pdf("test.pdf", height=7, width=8)
    #ggplot(x, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(fill=factor(Cell.Type)), size=2, shape=21, stroke=0)
    ggplot(x, aes(x= UMAP_1, y= UMAP_2)) + geom_point(aes(fill=factor(Cell.Type)), size=2, shape=21, stroke=0) + scale_fill_manual(values=color_flag, labels=color_label) + labs(fill="", colour="", stroke=0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.key = element_rect(fill = NA, colour = NA, size = 1)) + theme(text = element_text(size=14)) + guides(fill = guide_legend(override.aes = list(size=4)))
    #dev.off()
}

cell_number_png <- function(cells, filename, figdir, size=4, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    y <- summary(as.factor(cells$Sample))
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

featureplot_umap_png <- function(cells, features, filename, figdir, size=4, width=1200, height=1200){
    filename_umap <- paste(figdir, '/', filename, '.umap.png', sep="")
    png(filename_umap, width = width, height = height)
    font_size = 46
    FeaturePlot(cells, features = features, reduction='umap', ncol = 3) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    #dev.off()
}

featureplot_tsne_png <- function(cells, features, filename, figdir, size=4, width=1200, height=1200){
    filename_tsne <- paste(figdir, '/', filename, '.tsne.png', sep="")
    png(filename_tsne, width = width, height = height)
    font_size = 46
    FeaturePlot(cells, features = features, reduction='tsne', ncol = 3) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    #dev.off()
}


volinplot_png <- function(cells, features, group, title, filename, figdir, size=2, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    VlnPlot(cells, features = features, group.by = group, pt.size = 0.01) + labs(title=title, x="", y=features) +  theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    #dev.off()
}

doheatmap_png <- function(cells, features, filename, figdir, width=2500, height=3800){
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
figdir=paste(sample, platform, "seurat_after_singler_png", sep="_")

#pdf(paste(prefix, "_Seurat_2kgenes_vst_cc.featureplots.pdf", sep=""), height=7, width=10)

#load cell types for UMAP clusters from singleR
pbmc <- readRDS(paste(prefix, "_Seurat_2kgenes_vst_cc.rds", sep=""))
fig_cell_num = paste(sample, platform, "Cell_number", sep="_")
cell_number_png(pbmc, fig_cell_num, figdir)

celltype <- fread(paste(prefix, "_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt", sep=""))
rownames(celltype) <- celltype$Cell.ID
pbmc <- AddMetaData(pbmc, celltype)
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
#sample
#platform
#fig_cell_num = paste(sample, platform, "Cell_number", sep="_")
#cell_number_png(pbmc, fig_cell_num, figdir)

#fig_tsne = paste(sample, platform, 'tsne', "Encode_type_cluster", sep="_")
#dimplot_pdf(pbmc, celltype, 'tsne', fig_tsne, figdir, "Encode_main_cluster")

#fig_tsne = paste(sample, platform, 'umap', "Encode_type_cluster", sep="_")
#dimplot_pdf(pbmc, celltype, 'umap', fig_tsne, figdir, "Encode_main_cluster")


fig_rna = paste(sample, platform, "nCount_RNA", sep="_")
volinplot_png(pbmc, c("nCount_RNA"), "Sample", "", fig_rna, figdir)
fig_feature = paste(sample, platform, "nFeature_RNA", sep="_")
volinplot_png(pbmc, c("nFeature_RNA"), "Sample", "", fig_feature, figdir)
fig_mt = paste(sample, platform, "percent.mt", sep="_")
volinplot_png(pbmc, c("percent.mt"), "Sample", "", fig_mt, figdir)


fig_tsne = paste(sample, platform, 'umap', "ident", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "ident")
fig_tsne = paste(sample, platform, 'umap', "Sample", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Sample")
fig_tsne = paste(sample, platform, 'umap', "Encode_main_type", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Encode_main_type")
fig_tsne = paste(sample, platform, 'umap', "hpca_main_type", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "hpca_main_type")
fig_tsne = paste(sample, platform, 'umap', "Encode_type_cluster", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Encode_main_cluster")

fig_tsne = paste(sample, platform, 'tsne', "ident", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "ident")
fig_tsne = paste(sample, platform, 'tsne', "Sample", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Sample")
fig_tsne = paste(sample, platform, 'tsne', "Encode_main_type", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Encode_main_type")
fig_tsne = paste(sample, platform, 'tsne', "hpca_main_type", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "hpca_main_type")
fig_tsne = paste(sample, platform, 'tsne', "Encode_type_cluster", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Encode_main_cluster")


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
fig_immune_feature = paste(sample, platform, "Cell_Type_Immunecell_feature", sep="_")
featureplot_umap_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"), fig_immune_feature, figdir)
featureplot_tsne_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"), fig_immune_feature, figdir)
fig_immune_violin = paste(sample, platform, "Cell_Type_Immunecell_violin", sep="_")
volinplot_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"), "seurat_clusters", "Immune", fig_immune_violin, figdir)
#Fibroblast cell
#fig_fibroblast_feature = paste(sample, platform, "Fibroblast_feature", sep="_")
#featureplot_png(pbmc, c("PDPN", "THY1", "CD34", "CDH2", "CDH11", "COL1A1", "S100A4", "VIM", "ACTA2"), fig_fibroblast_feature, figdir)
#fig_fibroblast_violin = paste(sample, platform, "Fibroblast_violin", sep="_")
#volinplot_png(pbmc, c("PDPN", "THY1", "CD34", "CDH2", "CDH11", "COL1A1", "S100A4", "VIM", "ACTA2"), "seurat_clusters", fig_fibroblast_violin, figdir)
#epithelial cell
fig_epithelial_feature = paste(sample, platform, "Cell_Type_Epithelial_feature", sep="_")
featureplot_umap_png(pbmc, c("EPCAM", "ESR1", "ERBB2", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1"), fig_epithelial_feature, figdir)
featureplot_tsne_png(pbmc, c("EPCAM", "ESR1", "ERBB2", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1"), fig_epithelial_feature, figdir)
fig_epithelial_violin = paste(sample, platform, "Cell_Type_Epithelial_violin", sep="_")
volinplot_png(pbmc, c("EPCAM", "ESR1", "ERBB2", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1"), "seurat_clusters", "Epithelial", fig_epithelial_violin, figdir)

#Tcells
fig_tcell_feature = paste(sample, platform, "Cell_Type_Tcell_feature", sep="_")
featureplot_umap_png(pbmc, c("CD3D", "CD3E", "CD3G", "CD2"), fig_tcell_feature, figdir)
featureplot_tsne_png(pbmc, c("CD3D", "CD3E", "CD3G", "CD2"), fig_tcell_feature, figdir)
fig_tcell_violin = paste(sample, platform, "Cell_Type_Tcell_violin", sep="_")
volinplot_png(pbmc, c("CD3D", "CD3E", "CD3G", "CD2"), "seurat_clusters", "T cells", fig_tcell_violin, figdir)
#Bcells
fig_bcell_feature = paste(sample, platform, "Cell_Type_Bcell_feature", sep="_")
featureplot_umap_png(pbmc, c("CD19", "MS4A1", "CD79A", "CD79B"), fig_bcell_feature, figdir)
featureplot_tsne_png(pbmc, c("CD19", "MS4A1", "CD79A", "CD79B"), fig_bcell_feature, figdir)
fig_bcell_violin = paste(sample, platform, "Cell_Type_Bcell_violin", sep="_")
volinplot_png(pbmc, c("CD19", "MS4A1", "CD79A", "CD79B"), "seurat_clusters", "B cells", fig_bcell_violin, figdir)
#Macrophage
fig_macrophage_feature = paste(sample, platform, "Cell_Type_macrophage_feature", sep="_")
featureplot_umap_png(pbmc, c("CD14", "CD68", "CD163", "CSF1R", "FCGR3A"), fig_macrophage_feature, figdir)
featureplot_tsne_png(pbmc, c("CD14", "CD68", "CD163", "CSF1R", "FCGR3A"), fig_macrophage_feature, figdir)
fig_macrophage_violin = paste(sample, platform, "Cell_Type_macrophage_violin", sep="_")
volinplot_png(pbmc, c("CD14", "CD68", "CD163", "CSF1R", "FCGR3A"), "seurat_clusters", "Macrophage", fig_macrophage_violin, figdir)
#Fibroblast
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_feature", sep="_")
featureplot_umap_png(pbmc, c("FAP", "COL1A1", "COL3A1", "THY1", "ACTA2", "VIM", "CDH2", "CDH11", "PDPN"), fig_fibroblast_feature, figdir)
featureplot_tsne_png(pbmc, c("FAP", "COL1A1", "COL3A1", "THY1", "ACTA2", "VIM", "CDH2", "CDH11", "PDPN"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_violin", sep="_")
volinplot_png(pbmc, c("FAP", "COL1A1", "COL3A1", "THY1", "ACTA2", "VIM", "CDH2", "CDH11", "PDPN"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
#EMT
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_EMT_feature", sep="_")
featureplot_umap_png(pbmc, c("CDH1", "CRB3", "DSP", "CDH2", "FN1", "VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2"), fig_fibroblast_feature, figdir)
featureplot_tsne_png(pbmc, c("CDH1", "CRB3", "DSP", "CDH2", "FN1", "VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_EMT_violin", sep="_")
volinplot_png(pbmc, c("CDH1", "CRB3", "DSP", "CDH2", "FN1", "VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2"), "seurat_clusters", "EMT", fig_fibroblast_violin, figdir)
#NK cells
fig_nk_feature = paste(sample, platform, "Cell_Type_NK_feature", sep="_")
featureplot_umap_png(pbmc, c("NK", "GNLY", "CD56", "CD94", "NKP46"), fig_nk_feature, figdir)
featureplot_tsne_png(pbmc, c("NK", "GNLY", "CD56", "CD94", "NKP46"), fig_nk_feature, figdir)
fig_nk_violin = paste(sample, platform, "Cell_Type_NK_violin", sep="_")
volinplot_png(pbmc, c("NK", "GNLY", "CD56", "CD94", "NKP46"), "seurat_clusters", "NK cells", fig_nk_violin, figdir)
#PECAM1
fig_endothelial_feature = paste(sample, platform, "Cell_Type_endothelial_feature", sep="_")
featureplot_umap_png(pbmc, c("PECAM1", "CD34", "CD36", "ENTPD1", "CD44"), fig_endothelial_feature, figdir)
featureplot_tsne_png(pbmc, c("PECAM1", "CD34", "CD36", "ENTPD1", "CD44"), fig_endothelial_feature, figdir)
fig_endothelial_violin = paste(sample, platform, "Cell_Type_endothelial_violin", sep="_")
volinplot_png(pbmc, c("PECAM1", "CD34", "CD36", "ENTPD1", "CD44"), "seurat_clusters", "Endothelial", fig_endothelial_violin, figdir)



if(FALSE){
#find marker genes
cells.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#pdf
#DoHeatmap(object = pbmc, features = top10$gene)
#png
fig_markers = paste(sample, platform, "top10_genes", sep="_")
doheatmap_png(pbmc, top10$gene, fig_markers, figdir)
}

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

