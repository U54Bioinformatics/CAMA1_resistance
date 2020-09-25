library(data.table)
args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]   # before any filter
sample=args[2]
#prefix="VG_CAMA1_10x_gene_symbols.filtered.counts" #filtered
#counts.df <- fread("VG_CAMA1_count01.noaggr.txt")
#counts.df <- subset(counts.df, counts.df$'Gene Symbol' != "")
#counts.df <- unique(counts.df, by='Gene Symbol')
counts.df <- read.table(paste(prefix, sample, "txt", sep="."), sep="\t", header=T)
counts.df <- as.data.frame(counts.df)
rownames(counts.df) <- counts.df$`Gene.ID`
counts.df <- counts.df[,2:ncol(counts.df)]

mcherry_cells <- colnames(counts.df)[counts.df["mCherry",]>0]
mvenus_cells <- colnames(counts.df)[counts.df["mVenus",]>0]
write.table(mcherry_cells, paste(prefix, sample, "mcherry_cell.txt", sep="."), sep="\t", quote=F, row.names=F, col.names=F)
write.table(mvenus_cells, paste(prefix, sample, "mvenus_cell.txt", sep="."), sep="\t", quote=F, row.names=F, col.names=F)

shared <- intersect(mcherry_cells, mvenus_cells)
total <- dim(counts.df)[2]
print(paste("total", total, sep=","))
print(paste("mvenus_cells", length(mvenus_cells)-length(shared), (length(mvenus_cells)-length(shared))/total, sep=","))
print(paste("mcherry_cells", length(mcherry_cells)-length(shared), (length(mcherry_cells)-length(shared))/total, sep=","))
print(paste("shared", length(shared), length(shared)/total, sep=","))
print(paste("nomarker", dim(counts.df)[2]-(length(mvenus_cells)+length(mcherry_cells))+length(shared), (dim(counts.df)[2]-(length(mvenus_cells)+length(mcherry_cells))+length(shared))/total, sep=","))


