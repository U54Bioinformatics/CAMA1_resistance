#
# downstreamAnalysis.R
#
# Analysis to explore mechanisms of facilitation in single-cell RNA seq data.
#
# @author Ben Decato

################################################################################
# Libraries
################################################################################
rm(list=ls())

library(Seurat)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(ggrepel)

options(ggrepel.max.overlaps = Inf)

print(paste("Working directory:", getwd()))

# Set up my ggplot theme for use with all visualizations
myTheme <- theme_bw() +
  theme(text = element_text(size=16), axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.placement = "outside")

################################################################################
# Read and integrate data; high level UMAP plots
################################################################################
seu <- readRDS("~/Desktop/Facilitation/VG_CAMA1_D11_ALL_10x_Seurat_2kgenes_vst_cc.rds")
meta <- read.table("~/Desktop/Facilitation/VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.marker_genes.txt",
                   header = TRUE, sep = "\t")

seu@meta.data <- meta
row.names(seu@meta.data) <- seu@meta.data$Cell.orig

DimPlot(seu, reduction = "umap", 
        group.by = "Marker_groups", 
        label=FALSE, repel=TRUE, pt.size = 0.25, raster = TRUE) + 
  myTheme

FeaturePlot(seu, features = c("mCherry", "mVenus"),
            order=TRUE, pt.size=0.75, raster = TRUE)

################################################################################
## Differential expression between mono and cocultured sensitive cells
################################################################################
facilitation.de.markers <- FindMarkers(seu, ident.1 = "Sensitive (monoculture)",
                                       ident.2 = "Sensitive (coculture)")

# Ran once for every expressed gene; now just read the cached version
facilitation.de.markers <- FindMarkers(seu, ident.1 = "Sensitive (monoculture)",
                                       ident.2 = "Sensitive (coculture)", 
                               min.cells.group = 1, min.cells.feature = 1,
                               min.pct = 0, logfc.threshold = 0, only.pos = FALSE)
facilitation.de.markers <- facilitation.de.markers %>%
  mutate(Significant = ifelse(p_val_adj<0.05, "Yes", "No")) %>%
  rownames_to_column(var = "Gene")

write.table(facilitation.de.markers, file = "facilitation.de.markers.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
facilitation.de.markers <- read.table("~/Desktop/Facilitation/facilitation.de.markers.txt", 
                                      header = TRUE, sep = "\t")

facilitation.de.markers <- facilitation.de.markers %>%
  mutate(TopGenes = ifelse(abs(avg_log2FC)>0.75 & Significant == "Yes", "Yes", "No"))


ggplot(facilitation.de.markers, aes(x=-avg_log2FC, y = -log10(p_val),
                                    color = Significant)) + 
  geom_point(aes(size=abs(pct.2-pct.1))) +
  xlab("log2FC (+ = higher in coculture)") +
  ylab("-log10(p)") +
  geom_label_repel(aes(label = ifelse(TopGenes == "Yes", Gene, "")),
                   box.padding = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "red") +
  geom_vline(xintercept = -0.75, linetype = "dashed", color = "red") +
  myTheme +
  scale_color_npg()

################################################################################
## Exploring mechanistic components of estradiol synthesis
################################################################################

renaGenes <- c("HSD17B1", "HSD17B2", "HSD17B4", "HSD17B5", 
               "HSD17B6", "HSD17B7", "HSD3B1", "CYP17A1", 
               "CYP19A1", "CYP1B1", "CYP17A2")

geneSet <- read.table("~/Desktop/Facilitation/geneset.txt", header = TRUE,
                      sep = "\t")

FeaturePlot(seu, features = geneSet$KEGG_STEROID_HORMONE_BIOSYNTHESIS,
            order=TRUE, pt.size=0.75, raster = TRUE)

seu_subset <- subset(x = seu, subset = Marker_groups != "GV014_nomarker")

Idents(seu_subset) <- "Marker_groups"
levels(seu_subset) <- c("Sensitive (monoculture)", "Sensitive (coculture)",
                        "Resistant (coculture)", "Resistant (monoculture)")
DotPlot(object = seu_subset, 
        features = geneSet$KEGG_STEROID_HORMONE_BIOSYNTHESIS) +
  coord_flip() + 
  xlab("KEGG Steroid Hormone Biosynthesis Pathway Components") + 
  ylab("") +
  myTheme +
  scale_fill_npg()

cutDown <- c("HSD17B1", "HSD17B8", "SULT2B1")

VlnPlot(object = seu_subset, features = cutDown) +
  myTheme






################################################################################
## Assessing proliferation pathway activity across sensitivity and culture conditions
################################################################################

pathways <- readRDS("VG_CAMA1_D11_ALL_10x.zinbwave.normalized.ssGSEA.scores.RDS")
write.table(pathways, file = "full_pathways.txt", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)

pathwayList <- c("BIOCARTA_CELLCYCLE_PATHWAY", "REACTOME_CELL_CYCLE", "KEGG_CELL_CYCLE",
                 "WHITFIELD_CELL_CYCLE_S", "WHITFIELD_CELL_CYCLE_LITERATURE")

prolif_pathways <- pathways %>% 
  filter(`Gene Set` %in% pathwayList) %>% 
  pivot_longer(2:ncol(pathways), 
               names_to = "Cell.orig", values_to = "Score")

prolif_pathways <- left_join(prolif_pathways, seu@meta.data) %>% 
  select(`Gene Set`, Cell.orig, Score, Sample, Marker_genes, Marker_groups,
         Marker_genes1, Marker_groups1)

prolif_pathways$Marker_groups <- factor(prolif_pathways$Marker_groups,
                                        levels = c("Sensitive (monoculture)", 
                                                   "Sensitive (coculture)",
                                                   "Resistant (coculture)", 
                                                   "Resistant (monoculture)"))

ggplot(prolif_pathways, aes(x=`Gene Set`, y = Score, color = Marker_groups)) + 
  geom_boxplot(outlier.shape=NA, alpha = 0.2, color="black", aes(fill = Marker_groups)) +
  geom_point(size = 0.7, position=position_jitterdodge(jitter.width = 0.1)) + 
  stat_compare_means(aes(group = Marker_groups), label = "p.signif") +
  myTheme + 
  scale_color_npg()

results <- prolif_pathways %>%
  group_by(`Gene Set`) %>%
  do(w = kruskal.test(Score~Marker_groups, data=.)) %>%
  summarise(`Gene Set`, KW = w$p.value) %>%
  mutate(FDR = p.adjust(KW))







