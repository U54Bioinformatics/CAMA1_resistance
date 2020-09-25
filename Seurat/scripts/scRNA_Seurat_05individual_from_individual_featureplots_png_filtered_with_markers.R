library(Seurat)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(MAST)
library(future)
plan(strategy = "multicore", workers = 1)

cell_number_png <- function(cells, filename, figdir, size=4, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    y <- table(cells@meta.data["Sample"])
    print(y)
    df <- data.frame(y)
    print(df)
    names(df) <- c("Sample", "Cell_numbers")
    print(df)
    df <- df[df$Cell_numbers!=0,]
    print(df)
    p <- ggplot(data=df, aes(x = Sample, y = Cell_numbers, width=.5)) + labs(title="Cell numbers", x="", y="Cell numbers") + geom_bar(stat="identity") + geom_text(aes(x = Sample, y = Cell_numbers*1.1, label = Cell_numbers), size= font_size/2) + theme(plot.title=element_text(size=font_size,face="bold", hjust = 0.5), axis.text=element_text(size=font_size, face="bold"), axis.text.x = element_text(angle = 45, hjust = 1), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    print(p)
    dev.off()
}

dimplot_png <- function(cells, method, filename, figdir, group, size=1, legend_nrow=2, width=10, height=10){
    filename <- paste(figdir, '/', filename, '.pdf', sep="")
    pdf(filename, width = width, height = height)
    font_size = 23
    p <- DimPlot(cells, reduction = method, group.by = group, pt.size=size, label.size=12, label=FALSE) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.position = "top", legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(color = guide_legend(nrow = legend_nrow, byrow = TRUE, override.aes = list(size = 4)))
    print(p)
    dev.off()
}

featureplot_png <- function(cells, features, filename, figdir, size=4, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    p <- FeaturePlot(cells, features = features) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 28)))
    print(p)
    dev.off()
}

featureplot_umap_png <- function(cells, features, filename, figdir, size=4, width=1200, height=1200){
    filename_umap <- paste(figdir, '/', filename, '.umap.png', sep="")
    png(filename_umap, width = width, height = height)
    font_size = 46
    p <- FeaturePlot(cells, features = features, reduction='umap', ncol = 2) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 28)))
    print(p)
    dev.off()
}

featureplot_tsne_png <- function(cells, features, filename, figdir, size=4, width=1200, height=1200){
    filename_tsne <- paste(figdir, '/', filename, '.tsne.png', sep="")
    png(filename_tsne, width = width, height = height)
    font_size = 46
    p <- FeaturePlot(cells, features = features, reduction='tsne', ncol = 2) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 28)))
    print(p)
    dev.off()
}


volinplot_png <- function(cells, features, group, filename, figdir, size=2, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    p <- VlnPlot(cells, features = features, group.by = group) + labs(title=features, x="", y=features) +  theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    print(p)
    dev.off()
}


doheatmap_png <- function(cells, features, filename, figdir, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 26
    p <- DoHeatmap(pbmc, features = features) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size-10, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    print(p)
    dev.off()
}


args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
platform=args[2]
#sample="FEL011"
#platform="10x"
prefix=paste(sample, platform, sep="_")
figdir=paste(sample, platform, "seurat_filtered_marker_genes_png", sep="_")

#pdf(paste(prefix, "_Seurat_2kgenes_vst_cc.featureplots.pdf", sep=""), height=7, width=10)

#load cell types for UMAP clusters from singleR
pbmc <- readRDS(paste(prefix, "_Seurat_2kgenes_vst_cc.rds", sep=""))
celltype <- fread(paste(prefix, "_cell_metadata.UMAPcluster.marker_genes.txt", sep=""))
rownames(celltype) <- celltype$Cell.ID
pbmc <- AddMetaData(pbmc, celltype)
#pbmc <- AddMetaData(pbmc, celltype$Marker_Genes, col.name="Marker_genes")
#pbmc <- AddMetaData(pbmc, celltype$Marker_Groups, col.name="Marker_groups")
table(pbmc$orig.ident)
head(colnames(pbmc))
tail(colnames(pbmc))

#remove cells without marker
pbmc <- subset(pbmc, subset = Marker_genes != 'nomarker')
table(pbmc$orig.ident)
table(pbmc$Marker_genes)
table(pbmc$Marker_groups)
#keep cells with 1.5 to 2.5 k gene per cell
#pbmc <- subset(pbmc, subset = nFeature_RNA < 2500 & nFeature_RNA > 1500)

#figdir
fig_cell_num = paste(sample, platform, "Cell_number", sep="_")
cell_number_png(pbmc, fig_cell_num, figdir)

fig_rna = paste(sample, platform, "nCount_RNA", sep="_")
volinplot_png(pbmc, c("nCount_RNA"), "Sample", fig_rna, figdir)
fig_feature = paste(sample, platform, "nFeature_RNA", sep="_")
volinplot_png(pbmc, c("nFeature_RNA"), "Sample", fig_feature, figdir)
fig_mt = paste(sample, platform, "percent.mt", sep="_")
volinplot_png(pbmc, c("percent.mt"), "Sample", fig_mt, figdir)


fig_tsne = paste(sample, platform, 'umap', "ident", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "ident", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'umap', "RNA_snn_res.0.5", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "RNA_snn_res.0.5", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'umap', "RNA_snn_res.0.8", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "RNA_snn_res.0.8", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'umap', "RNA_snn_res.1.2", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "RNA_snn_res.1.2", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'umap', "RNA_snn_res.2", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "RNA_snn_res.2", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'umap', "Sample", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Sample", legend_nrow = 1)
fig_tsne = paste(sample, platform, 'umap', "Marker_genes", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Marker_genes", legend_nrow = 1)
fig_tsne = paste(sample, platform, 'umap', "Marker_groups", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Marker_groups", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'umap', "Cellcycle", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Phase", legend_nrow = 1)


#fig_tsne = paste(sample, platform, 'umap', "Encode_main_type", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Encode_main_type")
#fig_tsne = paste(sample, platform, 'umap', "hpca_main_type", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "hpca_main_type")
#fig_tsne = paste(sample, platform, 'umap', "Encode_type_cluster", sep="_")
#dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Encode_main_cluster")

fig_tsne = paste(sample, platform, 'tsne', "ident", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "ident", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'tsne', "RNA_snn_res.0.5", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "RNA_snn_res.0.5", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'tsne', "RNA_snn_res.0.8", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "RNA_snn_res.0.8", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'tsne', "RNA_snn_res.1.2", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "RNA_snn_res.1.2", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'tsne', "RNA_snn_res.2", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "RNA_snn_res.2", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'tsne', "Sample", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Sample", legend_nrow = 1)
fig_tsne = paste(sample, platform, 'tsne', "Marker_genes", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Marker_genes", legend_nrow = 1)
fig_tsne = paste(sample, platform, 'tsne', "Marker_groups", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Marker_groups", legend_nrow = 2)
fig_tsne = paste(sample, platform, 'tsne', "Cellcycle", sep="_")
dimplot_png(pbmc, 'umap', fig_tsne, figdir, "Phase", legend_nrow = 1)


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
#featureplot_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"), fig_immune_feature, figdir)
#fig_immune_violin = paste(sample, platform, "Immunecell_violin", sep="_")
#volinplot_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"), "seurat_clusters", fig_immune_violin, figdir)
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

#Markers
fig_epithelial_feature = paste(prefix, "marker_feature", sep="_")
featureplot_umap_png(pbmc, c("mVenus", "mCherry"), fig_epithelial_feature, figdir, width=1200, height=600)
featureplot_tsne_png(pbmc, c("mVenus", "mCherry"), fig_epithelial_feature, figdir, width=1200, height=600)
fig_epithelial_violin = paste(prefix, "marker_violin_RNA_snn_res.0.5", sep="_")
volinplot_png(pbmc, c("mVenus", "mCherry"), "RNA_snn_res.0.5", fig_epithelial_violin, figdir, width=1200, height=600)
fig_epithelial_violin = paste(prefix, "marker_violin_RNA_snn_res.0.8", sep="_")
volinplot_png(pbmc, c("mVenus", "mCherry"), "RNA_snn_res.0.8", fig_epithelial_violin, figdir, width=1200, height=600)
fig_epithelial_violin = paste(prefix, "marker_violin_RNA_snn_res.1.2", sep="_")
volinplot_png(pbmc, c("mVenus", "mCherry"), "RNA_snn_res.1.2", fig_epithelial_violin, figdir, width=1200, height=600)
fig_epithelial_violin = paste(prefix, "marker_violin_RNA_snn_res.2", sep="_")
volinplot_png(pbmc, c("mVenus", "mCherry"), "RNA_snn_res.2", fig_epithelial_violin, figdir, width=1200, height=600)

if (FALSE) {
#marker genes
Idents(object = pbmc) <- 'Marker_genes'
cells.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
#output txt file
genes.diff.file <- paste(sample, platform, "top10_genes_Marker_genes.gene_diff.txt", sep="_")
genes.id   <- data.table(gene = rownames(cells.markers))
genes.diff <- cbind(genes.id, cells.markers)
write.table(genes.diff, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#pdf
#DoHeatmap(object = pbmc, features = top10$gene)
#png
fig_markers = paste(sample, platform, "top10_genes_Marker_genes", sep="_")
doheatmap_png(pbmc, top10$gene, fig_markers, figdir)
#dev.off()

#marker groups
Idents(object = pbmc) <- 'Marker_groups'
cells.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
#output txt file
genes.diff.file <- paste(sample, platform, "top10_genes_Marker_groups.gene_diff.txt", sep="_")
genes.id   <- data.table(gene = rownames(cells.markers))
genes.diff <- cbind(genes.id, cells.markers)
write.table(genes.diff, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#pdf
#DoHeatmap(object = pbmc, features = top10$gene)
#png
fig_markers = paste(sample, platform, "top10_genes_Marker_groups", sep="_")
doheatmap_png(pbmc, top10$gene, fig_markers, figdir)
#dev.off()

#marker groups: GV014_mvenus vs. GV013_mvenus
Idents(object = pbmc) <- 'Marker_groups'
cells.markers <- FindMarkers(object = pbmc, ident.1="GV014_mvenus", ident.2="GV013_mvenus", min.pct = 0, logfc.threshold = 0, test.use = "MAST")
#output txt file
genes.diff.file <- paste(sample, platform, "top10_genes_GV014mV_vs_GV013mV.gene_diff.txt", sep="_")
genes.id   <- data.table(gene = rownames(cells.markers))
genes.diff <- cbind(genes.id, cells.markers)
write.table(genes.diff, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)
#filter significant
#genes.diff.sig  <- dplyr::filter(genes.diff, abs(avg_logFC) >= 0.25 & pct.1 >= 0.1 & pct.2 >= 0.1 & p_val_adj <= 0.05)
#genes.diff.file <- paste(sample, platform, "top10_genes_GV014mV_vs_GV013mV.gene_diff.significant.txt", sep="_")
#write.table(genes.diff.sig, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)

#top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#pdf
#DoHeatmap(object = pbmc, features = top10$gene)
#png
#fig_markers = paste(sample, platform, "top10_genes_GV014mV_vs_GV013mV", sep="_")
#doheatmap_png(pbmc, top10$gene, fig_markers, figdir)
#dev.off()
}

if (FALSE){
#marker groups: GV014_mcherry vs. GV015_mcherry
Idents(object = pbmc) <- 'Marker_groups'
cells.markers <- FindMarkers(object = pbmc, ident.1="GV014_mcherry", ident.2="GV015_mcherry", min.pct = 0, logfc.threshold = 0, test.use = "MAST")
#output txt file
genes.diff.file <- paste(sample, platform, "top10_genes_GV014mC_vs_GV015mC.gene_diff.txt", sep="_")
genes.id   <- data.table(gene = rownames(cells.markers))
genes.diff <- cbind(genes.id, cells.markers)
write.table(genes.diff, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)
#top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#pdf
#DoHeatmap(object = pbmc, features = top10$gene)
#png
#fig_markers = paste(sample, platform, "top10_genes_GV014mC_vs_GV015mC", sep="_")
#doheatmap_png(pbmc, top10$gene, fig_markers, figdir)
#dev.off()

#marker groups: GV013_mvenus vs. GV015_mcherry
Idents(object = pbmc) <- 'Marker_groups'
cells.markers <- FindMarkers(object = pbmc, ident.1="GV013_mvenus", ident.2="GV015_mcherry", min.pct = 0, logfc.threshold = 0, test.use = "MAST")
#output txt file
genes.diff.file <- paste(sample, platform, "top10_genes_GV013mV_vs_GV015mC.gene_diff.txt", sep="_")
genes.id   <- data.table(gene = rownames(cells.markers))
genes.diff <- cbind(genes.id, cells.markers)
write.table(genes.diff, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)
#top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#pdf
#DoHeatmap(object = pbmc, features = top10$gene)
#png
#fig_markers = paste(sample, platform, "top10_genes_GV013mV_vs_GV015mC", sep="_")
#doheatmap_png(pbmc, top10$gene, fig_markers, figdir)
#dev.off()

#marker groups: GV014_mvenus vs. GV014_mcherry
Idents(object = pbmc) <- 'Marker_groups'
cells.markers <- FindMarkers(object = pbmc, ident.1="GV014_mvenus", ident.2="GV014_mcherry", min.pct = 0, logfc.threshold = 0, test.use = "MAST")
#output txt file
genes.diff.file <- paste(sample, platform, "top10_genes_GV014mV_vs_GV014mC.gene_diff.txt", sep="_")
genes.id   <- data.table(gene = rownames(cells.markers))
genes.diff <- cbind(genes.id, cells.markers)
write.table(genes.diff, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)
#top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#pdf
#DoHeatmap(object = pbmc, features = top10$gene)
#png
#fig_markers = paste(sample, platform, "top10_genes_GV014mV_vs_GV014mC", sep="_")
#doheatmap_png(pbmc, top10$gene, fig_markers, figdir)
#dev.off()

#marker groups: GV014_mcherry vs. GV013_mvenus
Idents(object = pbmc) <- 'Marker_groups'
cells.markers <- FindMarkers(object = pbmc, ident.1="GV014_mcherry", ident.2="GV013_mvenus", min.pct = 0, logfc.threshold = 0, test.use = "MAST")
#output txt file
genes.diff.file <- paste(sample, platform, "top10_genes_GV014mC_vs_GV013mV.gene_diff.txt", sep="_")
genes.id   <- data.table(gene = rownames(cells.markers))
genes.diff <- cbind(genes.id, cells.markers)
write.table(genes.diff, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)
}

if (FALSE) {
#find marker genes
Idents(object = pbmc) <- 'RNA_snn_res.0.8'
cells.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
#output txt file
genes.diff.file <- paste(sample, platform, "top10_genes_RNA_snn_res.0.8.gene_diff.txt", sep="_")
genes.id   <- data.table(gene = rownames(cells.markers))
genes.diff <- cbind(genes.id, cells.markers)
write.table(genes.diff, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#png
fig_markers = paste(sample, platform, "top10_genes_RNA_snn_res.0.8", sep="_")
doheatmap_png(pbmc, top10$gene, fig_markers, figdir)
#dev.off()

#0.5
Idents(object = pbmc) <- 'RNA_snn_res.0.5'
cells.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
#output txt file
genes.diff.file <- paste(sample, platform, "top10_genes_RNA_snn_res.0.5.gene_diff.txt", sep="_")
genes.id   <- data.table(gene = rownames(cells.markers))
genes.diff <- cbind(genes.id, cells.markers)
write.table(genes.diff, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#png
fig_markers = paste(sample, platform, "top10_genes_RNA_snn_res.0.5", sep="_")
doheatmap_png(pbmc, top10$gene, fig_markers, figdir)
}
