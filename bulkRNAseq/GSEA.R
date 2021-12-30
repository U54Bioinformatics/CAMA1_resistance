library(fgsea)
library(data.table)
library(ggplot2)
library(GO.db)
library(org.Hs.eg.db)
library(tidyverse)
library(org.Sc.sgd.db)

cellline <- "CAMA.1"
condition <- "DMSO"
condition2 <- "ever"
condition3 <- "everR"
sample <- "COH069_15929R"
type <- "bulk_RNA"
#collection <- "H"
collection <- "C2"
###CHOOSE EXPERIMENT HERE
experiment <- "riboR_DMSO.vs.S_DMSO2"
prefix <- paste(sample, type, sep = "_")

#Load gene table from DE analysis
glm.res1 <- read.table(file = paste(prefix, cellline, experiment, genes,
                                              "DE_genes_glmQLFit_p.value1_edgR.txt", 
                                              sep = "_"),
                       sep = "\t", col.names = T, row.names = F, quote = F)

#Generate ranked gene matrix
myranks <- glm.res1 %>%
  #filter(PValue <= 0.05, .preserve = T) %>% # use all genes
  dplyr::select(Symbol, logFC) %>% 
  filter(!is.na(Symbol)) %>% 
  group_by(Symbol) %>% 
  filter(n() == 1) %>% 
  spread(key = Symbol, value = logFC)
ranks <- ncol(myranks)
myranks <- unlist(as.list(myranks))

#Pathways from msigdbr
pathways <- msigdbr(species = "Homo sapiens", category = collection) %>%
  group_by(gs_name)
pathways <- split(x = pathways$gene_symbol, f = pathways$gs_name)

#Run GSEA w/ pathways
coh069_fgsea <- fgsea(pathways = pathways, stats = myranks, minSize = 15, 
                      maxSize = 500, eps = 0.0, nPermSimple = 1000)
sum(is.na(coh069_fgsea))
coh069_fgsea <- coh069_fgsea %>%
  na.omit() %>%
  arrange(desc(NES))

#Isolate statistically significant enriched pathways (p < 0.25 - exploratory)
coh069_fgsea <- coh069_fgsea %>% 
  mutate(Enriched = ifelse(NES > 0, "POSITIVE_ENRICHMENT", 
                           "NEGATIVE_ENRICHMENT"))
cols <- c("POSITIVE_ENRICHMENT" = "blue", 
          "NEGATIVE_ENRICHMENT" = "grey")
summary(coh069_fgsea$padj <= 0.25)
sig <- coh069_fgsea %>%
  filter(padj <= 0.25) %>%
  arrange(desc(NES))

#
sig.up <- sig[! grep("DN", sig$pathway),] #only remove DN
r.sig.up <- dim(sig.up)[1]

#Save pathways
fwrite(sig.up, file=paste(prefix, cellline, experiment,"glmqlfit", 
                          genes, "ranked_genes",
                          collection, "collection_padj0.25",
                          r.sig.up, "pathways_filtered_fgsea.txt",
                          sep = "_"), 
       sep = "\t", col.names = T, row.names = F, quote = F)

#Look at top pathways
summary(sig.up$Enriched == "POSITIVE_ENRICHMENT")
summary(sig.up$Enriched == "NEGATIVE_ENRICHMENT")
top.up <- sig.up %>%
  head(30)
top.down <- sig.up %>% 
  tail(30)
top.all <- bind_rows(top.up, top.down) # use for C2 collection
top.all
## plot
ggplot(top.all, aes(reorder(pathway, NES), NES, fill = Enriched)) +
  geom_col(width = .9) +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 8)) +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.7, size = 14, face = "bold")) +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.6)) 
