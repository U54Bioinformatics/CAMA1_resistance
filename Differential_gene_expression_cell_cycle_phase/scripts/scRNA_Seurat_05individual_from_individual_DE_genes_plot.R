library(dplyr)
library(ggplot2)
library(ggrepel)

VolcanoPlot <- function(genes.diff, prefix, comparison, logFcThreshold = 0.25, adjPvalThreshold = 0.05) {
#################################################################
#genes.diff is a matrix with three columns: genes, avg_logFC, p_val_adj
#################################################################
pdf(paste0(prefix, ".VolcanoPlot.pdf"), width=8, height=7)
names(genes.diff) <- c("gene", "avg_logFC", "p_val_adj")
print(head(genes.diff))
print(logFcThreshold)
print(adjPvalThreshold)
print(prefix)
#vaconal plot
dataForVolcanoPlot <- genes.diff

dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                     avg_logFC < logFcThreshold | p_val_adj > adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                     avg_logFC >= logFcThreshold & p_val_adj <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                     avg_logFC <= -logFcThreshold & p_val_adj <= adjPvalThreshold)] <- 'DOWN'
fontsize=12
p <- ggplot(dataForVolcanoPlot, aes(x = avg_logFC, y = -log10(p_val_adj))) +
  labs(x=expression('log'[2]*'(Fold Change)'),
       y=(expression('-log'[10]*'(FDR)')),
       title=comparison) +
  geom_point(aes(color=Significance), alpha=1, size=2) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold),
             color='darkgreen',linetype='dashed')+
  #scale_x_continuous(breaks=c(-4,-2,0,2,4,6,8,10)) +
  #scale_y_continuous(expand = c(0.3, 0)) +
  #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
  scale_color_manual(values = c('green3',"black", "red")) +
  geom_text_repel(data = subset(dataForVolcanoPlot, 
                  avg_logFC >= logFcThreshold & -log10(p_val_adj) >= 200),
                  segment.alpha = 0.4, aes(label = gene),
                  size = 3.5, color = 'red', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot, 
                                avg_logFC <= -logFcThreshold & -log10(p_val_adj) >= 200),
                  segment.alpha = 0.4, aes(label = gene),
                  size = 3.5, color = 'green3', segment.color = 'black') +
  labs(title=prefix) +
    theme_classic() +
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1),
                     panel.background = element_blank()) +
    theme(axis.text=element_text(size=fontsize, color='black'),
          axis.title.x =element_text(size=fontsize, face="plain"),
          axis.title.y =element_text(size=fontsize, face="plain"),
          axis.title   =element_text(size=fontsize, face="plain"),
          strip.text = element_text(face = 'plain', size=12)) +
    theme(legend.text = element_text(size=fontsize),
          legend.title = element_blank(),
          legend.position = "none") +
    theme(plot.title=element_text(size=fontsize, face="plain", hjust = 0.95, vjust = 0), axis.text=element_text(size=fontsize, face="bold")) +
    theme(plot.margin =  margin(t = 1.5, r = 4, b = 1.5, l = 1, unit = "cm")) +
    scale_x_discrete(position = "bottom") +
    guides(colour = guide_legend(override.aes = list(size=4), ncol=1))

#by printing the graph we can now close device by dev.off() afterward
print(p)
dev.off()
}#end of function VolcanoPlot

#read DE genes test file
#prefix="VG_CAMA1_D11b_10x_top10_genes_GV014mV_vs_GV013mV"
args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]
comparison=args[2]
fc=as.numeric(args[3])
gene_diff_file=paste0(prefix, ".gene_diff.txt")
genes.diff <- read.table(gene_diff_file, header=T, sep="\t")

#VolcanoPlot
genes.diff.volcano <- genes.diff[c("gene", "avg_logFC", "p_val_adj")]
VolcanoPlot(genes.diff.volcano, prefix, comparison, logFcThreshold = fc, adjPvalThreshold = 0.05)

#filter significant genes
genes.diff.sig  <- dplyr::filter(genes.diff, abs(avg_logFC) >= fc & pct.1 >= 0.1 & pct.2 >= 0.1 & p_val_adj <= 0.05)
genes.diff.file <- paste0(prefix, ".gene_diff.significant.txt")
write.table(genes.diff.sig, genes.diff.file, sep="\t", quote=FALSE, row.names=FALSE)

