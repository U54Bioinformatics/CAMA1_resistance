library(GSVA)
library(data.table)
library(ggplot2)
library(tidyverse)
library(GO.db)
library(org.Hs.eg.db)
library(msigdbr)

#Definitions
cellline <- "CAMA.1"
condition <- "DMSO"
sample <- "COH069"
type <- "bulk_RNA"
collection <- "C2"
experiment <- "riboR_DMSO.vs.S_DMSO"
prefix <- paste(sample, type, sep = "_")
myTheme <- theme_bw() +
  theme(text = element_text(size=16), axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.placement = "outside")
 
 #Load logCPM
lcpm.filt <- read.table(file = paste(cellline, experiment,
                                     "lcpm.filtered", 
                                     r.lcpm.fil, "genes.txt", 
                                     sep = "_"),
                        sep = "\t", col.names = T, row.names = T, quote = F)

#msigdb pahtways
pathways <- msigdbr(species = "Homo sapiens", category = collection) %>%
  group_by(gs_name)
pathways <- split(x = pathways$gene_symbol, f = pathways$gs_name)

#run GSVA (ssGSEA)
gsva.c2 <- gsva(expr = as.matrix(lcpm.filt), gset.idx.list =  pathways,
                method="ssgsea",  min.sz = 15, max.sz = 500, verbose=TRUE, 
                kcdf="Poisson", parallel.sz=1)

#FIT DATA to model
fit3 <- lmFit(gsva.c2, design = design)
fit3 <- eBayes(fit3)
sigPathways1 <- topTable(fit3, coef=2, number=Inf, p.value=0.25, 
                        adjust.method ="BH", sort.by = "logFC", lfc = 0)
sigPathways.lfc <- sigPathways1
top.gsva <- gsva.c2[rownames(sigPathways.lfc),]
top.gsva2 <- top.gsva[!grepl("DN", rownames(top.gsva)),]

#Set colour for heatmap
require(RColorBrewer)
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
myBreaks <- seq(-1.5, 1.5, length.out=101)
heat1 <- t(scale(t(top.gsva2)))

#Heatmap
heatmap.2(heat1, 
          Rowv = T, 
          Colv = T, 
          dendrogram = "both",
          srtCol = 45,
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          col=myCol, breaks=myBreaks, 
          main=paste(cellline,experiment, "GSVA", sep = " "), 
          key=TRUE, keysize=1.0, key.title="", key.xlab="Enrichment Z-score", 
          scale="none", density.info="none",  trace="none", cexRow=0.6, 
          cexCol=0.8, 
          distfun=function(x) dist(x, method="euclidean"), 
          hclustfun=function(x) hclust(x, method="ward.D2"), margin=c(15,25))
          
          
#Save pathways
write.table(top.gsva2, 
            file=paste(prefix, cellline, experiment,
                       collection, 
                       "collection_lmFit.ebayes_padj0.25_filtered",
                       "pathways_gsva.txt", 
                       sep = "_"), 
            sep = "\t", col.names = T, row.names = T, quote = F)
