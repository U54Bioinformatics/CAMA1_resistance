##this script will read count and meta file from BETSY and create seurat obj for a patient that is specific from agrument
args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
platform=args[2]
#patient=args[3]
prefix=paste(sample, platform, sep="_")

library(dplyr)
library(Seurat)
library(data.table)

# Load the PBMC data from 10x count and meta, Lance
#FEL011_10x_count01.noaggr.counts.txt
#FEL011_10x_cell_metadata.txt
#Load exp and wellList
counts.df <- fread(paste(sample, platform, "count01.noaggr.counts.txt", sep="_"), sep="\t", nThread=1)
counts.df <- subset(counts.df, counts.df$'Gene Symbol' != "")
counts.df <- unique(counts.df, by='Gene Symbol')
counts.df <- as.data.frame(counts.df)
rownames(counts.df) <- counts.df$`Gene Symbol`
counts.df <- counts.df[,4:ncol(counts.df)]
#write.table(counts, paste(prefix, "gene_symbols.counts.txt", sep="_"), sep="\t", quote = F)
#Annotation
Annot.df <- read.table(paste(sample, platform, "cell_metadata.txt", sep="_"), header=TRUE, sep = "\t")
rownames(Annot.df) <- Annot.df$Cell
Annot.df$Cell.orig <- Annot.df$Cell
Annot.df <- Annot.df[,2:ncol(Annot.df)]
Annot.df  <- Annot.df[match(colnames(counts.df), rownames(Annot.df)),]
#Create Seurat Object
pbmc <- CreateSeuratObject(counts.df,  project = sample,  meta.data = Annot.df,  min.cells = 3, min.features = 0)
#tp   <- c("S", "M", "E")
#patient_tp <- paste(sample, tp, sep="_")
#pbmc <- subset(pbmc, subset = Sample == patient_tp[1] |  Sample == patient_tp[2] |  Sample == patient_tp[3])
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
nrow <- length(colnames(x=pbmc))
run_name <- rep(sample, nrow)
pbmc[["Run"]] <- run_name

###Write metadata
out <- data.table(Cell.ID = colnames(x=pbmc))
out_meta <- cbind(out, pbmc@meta.data)
fwrite(out_meta, paste(sample, ".cell_metadata.txt", sep=""), sep="\t", quote=F, col.names=T)

###Write count
counts <- GetAssayData(object = pbmc, slot = "counts")
counts <- as.matrix(counts)
out <- data.table(Gene.ID = rownames(x=pbmc))
out_counts <- cbind(out, counts)
fwrite(out_counts, paste(sample, ".expression_counts.txt", sep=""), sep="\t", quote = F)

saveRDS(pbmc, file = paste(sample, platform, "Seurat_2kgenes_vst_cc.raw.rds", sep="_"))

