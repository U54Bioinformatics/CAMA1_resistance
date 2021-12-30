library(edgeR)
library(data.table)
library(ggplot2)
library(tidyverse)

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
        
#Import Counts
counts <- read.delim(file = paste(sample, "15929R_data.counts.txt", sep = "_"),
                     header = T, quote = "", sep = "\t", row.names = "Gene.ID")

#Insert entrezID, filter NAs, match to counts
counts.df <- read.delim(file = paste(sample, "15929R_data.counts.txt", sep = "_"), 
                        header = T,quote = "", sep = "\t")
g <- counts.df$Gene.ID
g <- g %>% 
  as.matrix() %>% 
  as.data.frame()
annot <- grch38$symbol
annot <- annot %>%
  as.matrix() %>% 
  as.data.frame()
annot[["ENTREZID"]] <- grch38$entrez
colnames(annot) <- c("Symbol","ENTREZID")
annot2 <- annot[match(g$V1, annot$Symbol),]
sum(is.na(annot2))
annot3 <- na.omit(annot2)
sum(is.na(annot3))
counts <- counts[match(annot3$Symbol, rownames(counts)),]

#Generate counts containing only samples to be analyzed
## Here: All CAMA-1 Resistant DMSO vs. All CAMA-1 Sensitive DMSO
counts2 <- counts[,grep(condition, colnames(counts))]
##cellline
counts3 <- counts2[,grep(cellline, colnames(counts2))]
drop <- c("everR")
counts3 <- counts3[,-which(colnames(counts3) %in% drop)]

#Create edgeR object
group <- rep(c(1:2), each = 3)
y <- DGEList(counts = counts3, group = group, genes = annot3)
n <- rownames(y$samples)
n <- n %>%
  as.matrix() %>%
  as.data.frame() %>%
  separate(V1, into = c("CellLine", "R_S", "Treatment", "junk", "Replicate"), 
           sep = "_") %>%
  mutate(Replicate = ifelse(is.na(Replicate), junk, Replicate)) %>%
  mutate(junk = ifelse(junk == Replicate, Treatment, junk)) %>%
  mutate(Treatment = ifelse(Treatment == junk, R_S, Treatment)) %>% 
  mutate(R_S = ifelse(R_S == Treatment, "sensitive", R_S)) %>% 
  dplyr::select(-junk, -Replicate) %>%
  mutate(state = ifelse(R_S == "sensitive", "Sensitive", "Resistant"))
cellLine <- as.factor(n$CellLine)
R_S <- as.factor(n$R_S)
Treatment <- as.factor(n$Treatment)
state <- as.factor(n$state)
data.frame(Sample=colnames(y),cellLine,R_S,Treatment,state)

#Design matrix
design <- model.matrix(~0 + Treat_State, data = y$samples)

#Filer edgeR object based on state
keep <- filterByExpr(y, group = state)
y <- y[keep, , keep.lib.sizes = F]

#Normalize library
y <- calcNormFactors(y, method = "TMM")
nsamples <- ncol(y)
## logCPM counts
lcpm.filt <- cpm(y, log=TRUE)
r.lcpm.fil <- nrow(lcpm.filt)
## save logCPM counts
write.table(lcpm.filt, file = paste(cellline, experiment,
                                    "lcpm.filtered", 
                                    r.lcpm.fil, "genes.txt", 
                                    sep = "_"),
            sep = "\t", col.names = T, row.names = T, quote = F)

#Make comparisons: everR.DMSO.vs.S.DMSO
con <- makeContrasts(
  R.DMSO.vs.S.DMSO = Treat_StateDMSO_Resistant-Treat_StateDMSO_Sensitive,
                     levels = colnames(design))

#Estimate dispersion based on design 
y <- estimateDisp(y, design = design, robust = T)

#Fit and test: QLGLM
fit <- glmQLFit(y, design = design, robust = T)
qlf <- glmQLFTest(fit, contrast = con)

#Export top genes (change n and p.value as needed)
glm.res <- qlf %>% 
  topTags(n = Inf, adjust.method = "BH", sort.by = "logFC", p.value = 1)
glm.res1 <- glm.res$table

#Annotate significant genes
glm.res2 <- glm.res1
colnames(glm.res2) <- c("gene_symbol", "ENTRIZID", "logFC", 
                        "logCPM", "F", "PValue", "FDR")
##0.5 < logFC < -0.5, pval < 0.05
glm.res2 <- glm.res2 %>% 
  mutate(DE.logFC0.5 = ifelse(glm.res2$PValue < 0.05 & glm.res2$logFC > 0.5, 
                     "UP", "Not_Sig"))
glm.res2$DE.logFC0.5[glm.res2$PValue < 0.05 & glm.res2$logFC < -0.5] <- "DOWN" 
mycolors <- c( "blue","purple", "grey")
names(mycolors) <- c("UP", "DOWN", "Not_Sig")
##1.0 < logFC < -1.0, pval < 0.05
glm.res2 <- glm.res2 %>% 
  mutate(DE.logFC1.0 = ifelse(glm.res2$PValue < 0.05 & glm.res2$logFC > 1, 
                     "UP", "Not_Sig"))
glm.res2$DE.logFC1.0[glm.res2$PValue < 0.05 & glm.res2$logFC < -1] <- "DOWN" 

#Arrange and label top genes
glm.res2 <- glm.res2 %>%
  arrange(logFC)
h2 <- head(glm.res2$gene_symbol,30)
t2 <- tail(glm.res2$gene_symbol,40)
glm.res2 <- glm.res2 %>%
  arrange(FDR)
t <- head(glm.res2$gene_symbol, 20)
label <- c(t,h2,t2)

#Look at genes with logFC > 1 and pvalue < 0.05
glm.res2 <- glm.res2 %>% 
  mutate(genelabels = "FALSE") %>%
  mutate(genelabels = ifelse(gene_symbol %in% label, "TRUE", "FALSE"))
glm.res2 <- within(glm.res2,
                  genelabels[genelabels == "TRUE" & 
                               DE.logFC1.0== "Not_Sig"] <- "FALSE")

#Plot all genes: Volcano Plot
ggplot(data = glm.res2,  aes(x = logFC, y = -log10(PValue), 
                         colour=DE.logFC1.0,
                         label=ifelse(genelabels == TRUE, 
                               as.character(gene_symbol),""))) +
  geom_point(aes(size = abs(logCPM))) +
  geom_jitter() +
  geom_label_repel(label.size = 0.6,max.overlaps = Inf) +
  geom_vline(xintercept=c(-1.0,1.0), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  xlab("logFC1.0 (UP = Differentially Expressed in Resistant Cells Treated with
       Everolimus (BLUE))") + 
  theme_minimal() +
  ggtitle(paste(cellline, experiment,
                "DE Genes", sep = " "),) +
  theme(plot.title=element_text(hjust=0.5, size = 14, face = "bold")) +
  scale_color_manual(values=mycolors) +
  myTheme
  
#Plot genes w/ logFC > 1 and pvalue < 0.05: Heatmap
heatmap.2(heat,
   labCol=NULL,
   labRow = NULL,
   lhei=c(2,10),
   Rowv = T, 
   Colv = T,
   srtCol = 45,
   dendrogram = "both",reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
   col=myCol, 
   breaks=myBreaks, 
   main=paste(experiment,"lcpm Sig DEGs",
              sep = " "), 
   key=TRUE, keysize=1.0, key.title="", key.xlab="Enrichment Z-score", 
   scale="none", density.info="none",  trace="none", cexRow=0.8, 
   cexCol=1, distfun=function(x) dist(x, method="euclidean"), 
   hclustfun=function(x) hclust(x, method="ward.D2"), margin=c(10,25),
   revC = F,
   )
   
#Save genes as table (Use for GSEA)
fwrite(glm.res1, file = paste(prefix, cellline, experiment, genes,
                          "DE_genes_glmQLFit_p.value1_edgR.txt", 
                          sep = "_"),
       sep = "\t", col.names = T, row.names = F, quote = F)
