#LOAD LIBRARIES
library(data.table)
library(dplyr)
library(ggplot2)
library(Seurat)
library(circlize)
library(fastcluster)
library(RColorBrewer)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(rstatix)
library(gplots)
library(openxlsx)
library(ggridges)
options(ggrepel.max.overlaps = Inf)



#LOAD SEURAT OBJECT AND METADATA
seu <- readRDS(file = paste0("~/facilitation_scRNA/", 
                             "Facilitation/VG_CAMA1_D11_ALL_10x_Seurat_2kgenes_vst_cc.rds"))

meta <- read.table(file = paste0("~/facilitation_scRNA/Facilitation/", 
                                 "VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.", 
                                 "marker_genes.txt"), header = T, quote = "", sep = "\t")

seu@meta.data <- meta
row.names(seu@meta.data) <- seu@meta.data$Cell.orig

seu_subset <- subset(x = seu, 
                     subset = Marker_groups1 != "GV013_nomarker")
seu_subset <- subset(x = seu_subset, 
                     subset = Marker_groups1 != "GV014_nomarker")
seu_subset <- subset(x = seu_subset, 
                     subset = Marker_groups1 != "GV015_nomarker")


#UMAP OF MARKER GROUPS
DimPlot(seu_subset, reduction = "umap", group.by = "Marker_groups",
        label=FALSE, repel=TRUE, raster = FALSE,
        cols = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  theme_classic(base_size = 22)

#UMAP OF CELL CYCLE PHASE
DimPlot(seu_subset, reduction = "umap", group.by = "Phase",
        label=FALSE, repel=TRUE, raster = FALSE,
        cols = c("#D81B60", "#FFC107", "#1E88E5")) +
  theme_classic(base_size = 22)


#GET GENES FOR FIGURE S10
gene <- c("CDK1", "CCNB1", "CDKN2A", "TFF1", "TFF3", "GREB1")

marker.keep <- c("Sensitive (coculture)", "Sensitive (monoculture)")

#check counts in seu object
tt <- seu_subset@assays$RNA@scale.data %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>% 
  dplyr::filter(Gene %in% gene) %>%
  column_to_rownames(var = "Gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  add_column(Marker_groups = seu_subset@meta.data$Marker_groups) %>% 
  dplyr::filter(Marker_groups %in% marker.keep) %>% 
  dplyr::mutate_at(vars(Marker_groups), factor)


tt.tmp <- tt %>% 
  rownames_to_column(var = "IDs") %>% 
  remove_rownames()

setwd("~/facilitation_scRNA/reporting_summary/")
write.table(tt.tmp, file = "scRNA.seq_facil_supp.fig10a.b_Sen.co.mono_scaled_data.txt",
            sep = "\t", col.names = T, row.names = F, quote = F)

tt <- read.table(file = "scRNA.seq_facil_supp.fig10a.b_Sen.co.mono_scaled_data.txt",
                 sep = "\t", header = T, quote = "")

tt <- tt %>% 
  column_to_rownames(var = "IDs")

#MAKE FIGURE FOR EACH INDIVIDUAL GENE

#gene <- "CDK1"
#gene <- "CCNB1"
#gene <- "CDKN2A" 
#gene <- "TFF1"
#gene <- "TFF3" 
gene <- "GREB1"

#PERFORM NON-PARAMETRIC WILCOX RANK SUM
#REMOVE ZERO COUNTS
tt4 <- tt[which(tt$"GREB1" > 0), ]
f <- compare_means(GREB1~Marker_groups, tt4, method = "wilcox.test")
f

my_comparisons <- list(c("Sensitive (coculture)", "Sensitive (monoculture)"))

p <- ggplot(tt[which(tt$GREB1 > 0), ], 
            aes(
              x=Marker_groups,
              y=log(GREB1
                    + 1, base = 2),
              fill=Marker_groups
            )) +
  geom_violin(trim = F, scale = "width") +
  scale_fill_manual(values=c("#FFC107", "#004D40")) + #sensitive co/mono
  labs(title=gene) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  theme_classic() + 
  geom_boxplot(width=0.2, fill="white") +
  geom_jitter(shape=16, position=position_jitter(0.01)) +
  ylab("Expression Level (log2(counts + 1)") +
  theme(axis.title.y = element_text(size = 8, face = "bold")) +
  xlab("") +
  theme(axis.text.x = element_text(size = 19,
                                   angle = 0, 
                                   hjust = 0.5,
                                   face = "bold")) +
  theme(axis.text.y = element_text(size = 14,
                                   angle = 0, 
                                   face = "bold")) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20)) +
  stat_compare_means(comparisons = my_comparisons, 
                     label = "p.adj",
                     method = "wilcox.test", size = 7) +
  NoLegend()
p

setwd("~/facilitation_scRNA/reporting_summary/figures/")
tiff(filename = paste0(gene, "facilitation_figS10.tiff", 
                       sep = "_"), width = 7, height = 6, units = "in", res = 300)
p
dev.off()

#apply(tt.tmp[,2:7], 2, function (x) {print(min(x))
#  print(max(x))
#  print(quantile(x))}
#)


#LOAD SSGSEA SCORES AND MAKE FIGURE FOR PANELS C AND D
setwd(dir = "~/facilitation_scRNA/Facilitation/")
pathways <- readRDS("VG_CAMA1_D11_ALL_10x.zinbwave.normalized.ssGSEA.scores.RDS")


#Put pathway names in data frame
p.names <- pathways %>%
  dplyr::select("Gene Set")

pathwayList <- c("KEGG_CELL_CYCLE", "REACTOME_CELL_CYCLE", 
                 "WILLIAMS_ESR1_TARGETS_UP",
                 "DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP")


#check counts in seu object
dd <- seu_subset@meta.data %>%
  as.data.frame() %>%
  pull(Cell.ID) %>% 
  as.character()

prolif_pathways <- pathways %>%
  dplyr::filter(`Gene Set` %in% pathwayList) %>% 
  pivot_longer(! "Gene Set", 
               names_to = "Cell.orig", values_to = "Score") %>% 
  dplyr::filter(Cell.orig %in% tt) %>% 
  pivot_wider(names_from = "Cell.orig", values_from = "Score") %>% 
  column_to_rownames(var = "Gene Set")

all(colnames(dd) == colnames(prolif_pathways))

prolif_pathways1 <- prolif_pathways[,match(dd, colnames(prolif_pathways))]
all(colnames(prolif_pathways1) == seu_subset@meta.data$Cell.ID)

prolif_pathways2 <- prolif_pathways1 %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  add_column(Marker_groups = seu_subset@meta.data$Marker_groups)


all(rownames(prolif_pathways2) == seu_subset@meta.data$Cell.ID)

prolif_pathways2$Marker_groups <- factor(prolif_pathways2$Marker_groups,
                                         levels = c(
                                           "Resistant (coculture)",
                                           "Resistant (monoculture)",
                                           "Sensitive (coculture)",
                                           "Sensitive (monoculture)"
                                         ))



#ggplot(prolif_pathways2, aes(x=KEGG_CELL_CYCLE, y = Marker_groups, 
#ggplot(prolif_pathways2, aes(x=REACTOME_CELL_CYCLE, y = Marker_groups, 
ggplot(prolif_pathways2, aes(x=WILLIAMS_ESR1_TARGETS_UP, y = Marker_groups, 
                             #ggplot(prolif_pathways2, aes(x=DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP, y = Marker_groups, 
                             fill = Marker_groups)) +
  geom_density_ridges(quantile_lines = T, quantiles = 2)+
  scale_fill_manual(values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  theme_classic()

#SAVE

fig4d2 <- prolif_pathways2 %>% 
  dplyr::select(DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP, Marker_groups)
fig4c1 <- prolif_pathways2%>% 
  dplyr::select(KEGG_CELL_CYCLE, Marker_groups)
fig4c2 <- prolif_pathways2%>% 
  dplyr::select(REACTOME_CELL_CYCLE, Marker_groups)
fig4d1 <- prolif_pathways2%>% 
  dplyr::select(WILLIAMS_ESR1_TARGETS_UP, Marker_groups)


#GET UMAP CORRDINATES
umap_coor <- seu_subset[["umap"]]@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Cell.ID")

all(umap_coor$Cell.ID == seu_subset@meta.data$Cell.ID)
umap_coor <- umap_coor %>% 
  add_column(Marker_groups = seu_subset@meta.data$Marker_groups) %>% 
  add_column(Phase = seu_subset@meta.data$Phase)

ggplot(umap_coor, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = Marker_groups)) +
  scale_color_manual(values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  theme_classic()

ggplot(umap_coor, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = Phase)) +
  scale_color_manual(values = c("#D81B60", "#FFC107", "#1E88E5")) +
  theme_classic()

fig4a <- umap_coor %>% 
  dplyr::select(- Phase) %>% 
  remove_rownames()
fig4b <- umap_coor %>% 
  dplyr::select(- Marker_groups)%>% 
  remove_rownames()

setwd("~/facilitation_scRNA/reporting_summary/")
counts.hsd <- read.table(paste0("facilitation_CAMA.1_R.DMSO_S.DMSO_HSD17B",
                                "_genes_logCPM_normalized_counts.txt"), 
                         sep = "\t", header = T, quote = "")
counts.hsd <- data.frame(t(counts.hsd))
counts.hsd$Iso <- c(rep("Sensitive", 3), rep("Resistant",3))

ggplot(counts.hsd, aes(y=HSD17B1, x=Iso)) +
  geom_point(size=8, aes(color=Iso, alpha=0.85)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="pointrange", color="grey10", size=1, alpha=0.85) +
  scale_color_manual(values = c("#1E88E5", "#004D40")) +
  theme_classic(base_size = 22) +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x="", y="log CPM", title="HSD17B1")

ggplot(counts.hsd, aes(y=HSD17B8, x=Iso)) +
  geom_point(size=8, aes(color=Iso, alpha=0.85)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="pointrange", color="grey10", size=1, alpha=0.85) +
  scale_color_manual(values = c("#1E88E5", "#004D40")) +
  theme_classic(base_size = 22) +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x="", y="log CPM", title="HSD17B8")
dev.off()

t.test(HSD17B1 ~ Iso, data = counts.hsd)
t.test(HSD17B8 ~ Iso, data = counts.hsd)


fig4e1 <- counts.hsd %>% 
  rownames_to_column(var = "Cell.ID") %>% 
  dplyr::select(HSD17B1, Iso)

fig4e2 <- counts.hsd %>% 
  rownames_to_column(var = "Cell.ID") %>% 
  dplyr::select(HSD17B8, Iso)


#MAKE LIST AND SAVE EACH ONTO OWN SHEET

list <- list("fig4a" = fig4a,
             "fig4b" = fig4b,
             "fig4c1" = fig4c1,
             "fig4c2" = fig4c2,
             "fig4d1"= fig4d1,
             "fig4d2" = fig4d2,
             "fig4e1"= fig4e1,
             "fig4e2" = fig4e2,
             "fig10a.b" = tt.tmp)

setwd(dir = "~/facilitation_scRNA/reporting_summary/")
openxlsx::write.xlsx(list, file = "scRNA.seq_facilitation_data_for_Fig4.xlsx")
