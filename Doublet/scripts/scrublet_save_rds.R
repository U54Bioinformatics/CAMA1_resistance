library(scater)
library(Seurat)
library(scran)
library(BiocSingular)
library(pheatmap)
library(data.table)
library(scds)

args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]

#read seurat obj and convert to scater obj
pbmc <- readRDS(file = paste0(prefix, "_Seurat_2kgenes_vst_cc.raw.rds"))

#read scrublet anno
celltype <- fread(paste(prefix, ".scrublet_out.csv", sep=""))
rownames(celltype) <- celltype$Cell.ID
head(celltype)
pbmc <- AddMetaData(pbmc, celltype)

#remove doublet
pbmc <- subset(pbmc, subset = scrublet_doublet_call2 == 'FALSE')

###Write metadata
out <- data.table(Cell.ID = colnames(x=pbmc))
out_meta <- cbind(out, pbmc@meta.data)
fwrite(out_meta, paste(prefix, "Seurat_2kgenes_vst_cc.doublet_free.cell_metadata.txt", sep="_"), sep="\t", quote=F, col.names=T)

#Write RDS
saveRDS(pbmc, file = paste(prefix, "Seurat_2kgenes_vst_cc.doublet_free.rds", sep="_"))

