library(data.table)
args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]
prefix2=args[2]

#prefix="VG_CAMA1_10x_gene_symbols.filtered.counts" #filtered
#counts.df <- fread("VG_CAMA1_count01.noaggr.txt")
#counts.df <- subset(counts.df, counts.df$'Gene Symbol' != "")
#counts.df <- unique(counts.df, by='Gene Symbol')
counts.df <- read.table(paste(prefix, ".GV014.txt", sep=""), sep="\t", header=T)
counts.df <- as.data.frame(counts.df)
rownames(counts.df) <- counts.df$`Gene.ID`
counts.df <- counts.df[,2:ncol(counts.df)]

mcherry_cells <- colnames(counts.df)[counts.df["mCherry",]>0 & counts.df["mVenus",]==0]
mvenus_cells <- colnames(counts.df)[counts.df["mVenus",]>0 & counts.df["mCherry",]==0]
nomarker_cells <- colnames(counts.df)[counts.df["mVenus",]==0 & counts.df["mCherry",]==0]
doublet_cells <- colnames(counts.df)[counts.df["mVenus",]>0 & counts.df["mCherry",]>0]

write.table(mcherry_cells, paste(prefix, ".GV014.mcherry_cell.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
write.table(mvenus_cells, paste(prefix, ".GV014.mvenus_cell.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)

total <- dim(counts.df)[2]
print(paste("total", total, sep=","))
print(paste("mvenus_cells", length(mvenus_cells), length(mvenus_cells)/total, sep=","))
print(paste("mcherry_cells", length(mcherry_cells), length(mcherry_cells)/total, sep=","))
print(paste("shared", length(doublet_cells), length(doublet_cells)/total, sep=","))
print(paste("nomarker", length(nomarker_cells), length(nomarker_cells)/total, sep=","))

#compare gene numebrs, mt.percent, UMI numbers across groups
library(dplyr)
library(magrittr)

meta <- read.table(paste0(prefix2, "_10x_cell_metadata.UMAPcluster.txt"), sep="\t", header=T)
prepare <- function(name, value, xname = x_name, yname = y_name) {
  data_frame(rep(name, length(value)), value) %>%
    set_colnames(c("Cell.Groups", "Cell.ID"))
}

bind_rows(
  prepare("resistant", mcherry_cells),
  prepare("sensitive", mvenus_cells),
  prepare("nomarker", nomarker_cells),
  prepare("both", doublet_cells)
) -> groups.df

groups.df <- as.data.frame(groups.df)
rownames(groups.df) <- groups.df$Cell.ID
meta_groups.df <- merge(meta, groups.df, by="Cell.ID")
crosstable <- xtabs(~Cell.Groups+seurat_clusters, data=meta_groups.df)
write.table(meta_groups.df, paste0(prefix2, "_10x_cell_metadata.UMAPcluster.cell_groups.violin.data.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(crosstable, paste0(prefix2, "_10x_cell_metadata.UMAPcluster.cell_groups.crosstable.txt"), sep="\t", quote=F, col.names=T)

#violin plot
library(ggplot2)
library(stringr)
colors <- c("#E86A10","#56A4AA", "#3A78C4","#F1AB00")
axis_size = 15
ggformat <- theme_classic() + theme(plot.title=element_text(size=24,face="bold"), axis.text=element_text(size=24, face="bold"), axis.title=element_text(size=24,face="bold"), axis.text.x=element_blank(),
 legend.text=element_text(size=24, face="bold"), legend.title=element_text(size=24, face="bold")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
xlab_text = "Cell groups" 
title_text = ""


pdf(paste0(prefix2, "_10x_cell_metadata.UMAPcluster.cell_groups.violin.pdf"), width=8, height=5)
x <- read.table(paste0(prefix2, "_10x_cell_metadata.UMAPcluster.cell_groups.violin.data.txt"), header = T, sep="\t")
#nCount_RNA
p <- ggplot(x, aes(x=x$Cell.Groups, y=x$nCount_RNA, fill=x$Cell.Groups)) + scale_fill_manual(values=colors) + labs(title=title_text, x=xlab_text, y = "nCount_RNA")+ geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1), size=0.5)
p + ggformat
#nFeature_RNA
p <- ggplot(x, aes(x=x$Cell.Groups, y=x$nFeature_RNA, fill=x$Cell.Groups)) + scale_fill_manual(values=colors) + labs(title=title_text, x=xlab_text, y = "nFeature_RNA")+ geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1), size=0.5)
p + ggformat
#Percent.Mitochondria
p <- ggplot(x, aes(x=x$Cell.Groups, y=x$Percent.Mitochondria, fill=x$Cell.Groups)) + scale_fill_manual(values=colors) + labs(title=title_text, x=xlab_text, y = "mt.percent")+ geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1), size=0.5)
p + ggformat
#S
p <- ggplot(x, aes(x=x$Cell.Groups, y=x$S, fill=x$Cell.Groups)) + scale_fill_manual(values=colors) + labs(title=title_text, x=xlab_text, y = "S")+ geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1), size=0.5)
p + ggformat
#G2M
p <- ggplot(x, aes(x=x$Cell.Groups, y=x$G2M, fill=x$Cell.Groups)) + scale_fill_manual(values=colors) + labs(title=title_text, x=xlab_text, y = "G2M")+ geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1), size=0.5)
p + ggformat
dev.off()



