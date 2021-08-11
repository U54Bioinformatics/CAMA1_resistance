library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)

#https://yulab-smu.github.io/clusterProfiler-book/chapter3.html#msigdb-analysis
#https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
#export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_enricher/lib/R/library

Enrichment_test <- function (geneList, genes, prefix, m_t2g, pvalue=0.05){
  #enrichment test by hypergeometric test
  em_hyper_h <- enricher(genes, universe=names(geneList), TERM2GENE=m_t2g)
  write.table(em_hyper_h, paste0(prefix, ".hypergeometric.txt"), sep="\t", quote=F, row.names=F)
  
  #gene set enrichment analyses
  em_gsea_h <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = pvalue)
  write.table(em_gsea_h, paste0(prefix, ".GSEA.txt"), sep="\t", quote=F, row.names=F)
  return(em_gsea_h)
}

Enrichment_plot <- function (em_gsea_h, prefix, comparison, reverse = 0){
  em_gsea_h <- as.data.frame(em_gsea_h)
  print(head(em_gsea_h))
  em_gsea_h[c("ID", "NES", "p.adjust")]
  print(head(em_gsea_h))
  em_gsea_h.plot <- em_gsea_h[c("ID", "NES", "p.adjust")]
  if ( reverse ){
     em_gsea_h.plot$NES <- 0 - em_gsea_h.plot$NES
  }
  em_gsea_h.plot$NES_abs <- abs(em_gsea_h.plot$NES)
  em_gsea_h.plot$NES_change[with(em_gsea_h.plot, NES > 0)] <- "Positive"
  em_gsea_h.plot$NES_change[with(em_gsea_h.plot, NES < 0)] <- "Negative"
  em_gsea_h.plot$Comparison <- comparison
  em_gsea_h.plot$GS <- em_gsea_h.plot$ID
  em_gsea_h.plot <- em_gsea_h.plot[order(em_gsea_h.plot$NES),]
  my_level <- em_gsea_h.plot$GS
  em_gsea_h.plot$GS <- factor(em_gsea_h.plot$GS, levels = my_level)
  #plot
  pdf(paste0(prefix, ".Hallmark_NES.pdf"))
  fontsize=12
  p <- ggplot(em_gsea_h.plot, mapping=aes(x=Comparison, y=GS, color=NES_change, size=NES_abs)) +
  geom_point(shape=19) +
    labs(title="Normalized enrichment scores (NES)", x="", y = "") +
    scale_color_manual(values=c("blue", "red")) +
    theme_classic() +
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1),
                     panel.background = element_blank()) +
    theme(axis.text=element_text(size=fontsize, color='black'),
          axis.text.x =element_text(size=fontsize, angle = 0, color='black', hjust = 0.5, vjust = 1),
          axis.title.x =element_blank(),
          axis.title.y =element_text(size=fontsize, face="plain"),
          axis.title   =element_text(size=fontsize, face="plain"),
          strip.text = element_text(face = 'plain', size=12)) +
    theme(legend.text = element_text(size=fontsize),
          legend.title = element_blank(),
          legend.position = "bottom") +
    theme(plot.title=element_text(size=fontsize, face="plain", hjust = 0.95, vjust = -138), axis.text=element_text(size=fontsize, face="bold")) +
    theme(plot.margin =  margin(t = 1.5, r = 4, b = 1.5, l = 1, unit = "cm")) +
    scale_x_discrete(position = "top") +
    guides(colour = guide_legend(override.aes = list(size=4), ncol=1))
    print(p)
    dev.off()
}

#prefix="VG_CAMA1_D11b_10x_top10_genes_GV014mV_vs_GV013mV"
#comparison="GV014mV_vs_GV013mV"
args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]
comparison=args[2]
fc=as.numeric(args[3])
reverse=as.numeric(args[4])
gene_diff_file=paste0(prefix, ".gene_diff.txt")
genes.diff <- read.table(gene_diff_file, header=T, sep="\t")
genes.diff.volcano <- genes.diff[c("gene", "avg_logFC", "p_val_adj")]
head(genes.diff.volcano)

#create ranked gene list
geneList <- genes.diff.volcano$avg_logFC
names(geneList) <- as.character(genes.diff.volcano$gene)
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)

#create DE gene list
genes.diff.sig  <- dplyr::filter(genes.diff, abs(avg_logFC) >= fc & pct.1 >= 0.1 & pct.2 >= 0.1 & p_val_adj <= 0.05)
genes <- genes.diff.sig$gene

#read MSigDB
m_t2g_hallmark <- read.gmt("/home/jichen/Projects/Database/MSigDB/h.all.v6.2.symbols.gmt")
m_t2g_c2 <- read.gmt("/home/jichen/Projects/Database/MSigDB/c2.all.v6.2.symbols.gmt")
m_t2g_c7 <- read.gmt("/home/jichen/Projects/Database/MSigDB/c7.all.v6.2.symbols.gmt")

#perform test
#hallmark
try({
title = paste0(prefix, ".hallmark")
m_t2g = m_t2g_hallmark
em_gsea_h = Enrichment_test(geneList, genes, title, m_t2g, 0.05)
Enrichment_plot(em_gsea_h, prefix, comparison, reverse)
})
#c2
try({
title = paste0(prefix, ".c2")
m_t2g = m_t2g_c2
Enrichment_test(geneList, genes, title, m_t2g, 0.05)
})
#c7
try({
title = paste0(prefix, ".c7")
m_t2g = m_t2g_c7
Enrichment_test(geneList, genes, title, m_t2g, 0.05)
})

