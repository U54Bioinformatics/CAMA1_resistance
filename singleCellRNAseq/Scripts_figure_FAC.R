library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(ggridges)

seu <- readRDS("facilitation_seu_subset.rds")  
ssgsea <- readRDS("VG_CAMA1_D11_ALL_10x.zinbwave.normalized.ssGSEA.scores.RDS")

rownames(ssgsea) <- ssgsea$`Gene Set`
ssgsea <- as.matrix(ssgsea[, 2:ncol(ssgsea)])

dat <- data.frame(table(seu$Marker_groups, seu$Phase))


### Stacked bar plot of cell cycle phases ###
pdf("Cell_cycle_proportions.pdf", width=8, height=5)
ggplot(dat, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("darkred", "orange", "blue")) +
  theme_classic(base_size = 22)
dev.off()

### Which resistant co-culture pathways are different from other cells
x1 <- which(seu$Marker_groups == "Resistant (coculture)")
x2 <- which(seu$Marker_groups != "Resistant (coculture)")

y1 <- match(seu$Cell.ID[x1], colnames(ssgsea))
y2 <- na.omit(match(seu$Cell.ID[x2], colnames(ssgsea)))

FC <- array()
FC.p <- array()
for (i in 1:nrow(ssgsea)) {
  FC[i] <- mean(ssgsea[i, y1]) / mean(ssgsea[i, y2])
  FC.p[i] <- wilcox.test(ssgsea[i, y1], ssgsea[i, y2])$p.value
  print(i)
}

FC.fdr <- p.adjust(FC.p, method = 'fdr')

write.csv(data.frame('Gene.Set'=rownames(ssgsea),
                     'FC'=FC,
                     'Wilcox.P'=FC.p,
                     'FDR'=FC.fdr), 
          file="Supp_Table_Wilcox_ssgsea_Resistant_co_vs_all.csv")

### Which sensitive co-culture pathways are different from monoculture
x1 <- which(seu$Marker_groups == "Sensitive (coculture)")
x2 <- which(seu$Marker_groups == "Sensitive (monoculture)")

y1 <- match(seu$Cell.ID[x1], colnames(ssgsea))
y2 <- na.omit(match(seu$Cell.ID[x2], colnames(ssgsea)))

FC2 <- array()
FC2.p <- array()
for (i in 1:nrow(ssgsea)) {
  FC2[i] <- mean(ssgsea[i, y1]) / mean(ssgsea[i, y2])
  FC2.p[i] <- wilcox.test(ssgsea[i, y1], ssgsea[i, y2])$p.value
  print(i)
}

FC2.fdr <- p.adjust(FC2.p, method = 'fdr')

write.csv(data.frame('Gene.Set'=rownames(ssgsea),
                     'FC'=FC2,
                     'Wilcox.P'=FC2.p,
                     'FDR'=FC2.fdr), 
          file="Supp_Table_Wilcox_ssgsea_Sensitive_co_vs_mono.csv")

### 
sen.co <- readRDS("facilitation_sen.co.fc2.rds")
sen.mono <- readRDS("facilitation_sen.mono.fc2.rds")

dim(sen.co); dim(sen.mono)
mat <- cbind(sen.co, sen.mono)


heatmap(mat, Colv=NA, scale = "row", col=colorRampPalette(brewer.pal(8, "PiYG"))(25))


### Violin plots of differentially expressed genes 
marker.meta <- read.delim("VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.marker_genes.txt", 
                          header=T, row.names = 1)

marker.exp <- as.matrix(read.delim("VG_CAMA1_D11_ALL_10x.expression_counts.normalized.txt", 
                          header=T, row.names = 1))

marker.meta <- marker.meta[match(colnames(marker.exp), marker.meta$Cell.orig), ]

marker.data <- cbind(t(marker.exp), marker.meta)
marker.data.sen <- marker.data[grep("Sensitive", marker.meta$Marker_groups), ]



k <- which(marker.data.sen$CDK1 > 0)
ggplot(marker.data.sen[k,], aes(y=CDK1, x=Marker_groups, fill=Marker_groups)) +
  geom_violin() +
  geom_boxplot(width=0.15, fill="white")


k <- which(marker.data.sen$CCNB1 > 0)
ggplot(marker.data.sen[k,], aes(y=CCNB1, x=Marker_groups, fill=Marker_groups)) +
  geom_violin() +
  geom_boxplot(width=0.15, fill="white")

k <- which(marker.data.sen$CDKN2A > 0)
ggplot(marker.data.sen[k,], aes(y=CDKN2A, x=Marker_groups, fill=Marker_groups)) +
  geom_violin() +
  geom_boxplot(width=0.15, fill="white")


ggplot(marker.data.sen[k,], aes(x=CDK1, y=Marker_groups)) +
  geom_density_ridges(scale = 4) + 
  #scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  #scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges()



#### Plot interesting ssGSEA scores ####
y <- match(colnames(ssgsea), seu$Cell.ID)
ssgsea.dat <- data.frame(t(ssgsea[c("WILLIAMS_ESR1_TARGETS_UP", "DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP", 
                "KEGG_CELL_CYCLE", "REACTOME_CELL_CYCLE") ,]), "Group"=seu$Marker_groups[y])

pdf("Interesting_Pathways_1.pdf", height=3, width=9)
ggplot(ssgsea.dat, aes(y=Group, x=KEGG_CELL_CYCLE, fill=Group)) +
  geom_density_ridges() +  
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, size=1, color="white") +
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) + 
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  labs(y="")

ggplot(ssgsea.dat, aes(y=Group, x=REACTOME_CELL_CYCLE, fill=Group)) +
  geom_density_ridges() +  
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, size=1, color="white") +
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) + 
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  labs(y="")

ggplot(ssgsea.dat, aes(y=Group, x=WILLIAMS_ESR1_TARGETS_UP, fill=Group)) +
  geom_density_ridges() +  
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, size=1, color="white") +
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) + 
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  labs(y="")

ggplot(ssgsea.dat, aes(y=Group, x=DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP, fill=Group)) +
  geom_density_ridges() +  
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, size=1, color="white") +
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) + 
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  labs(y="")

dev.off()

###### Make plots from bulk RNAseq data ######
br <- read.delim("facilitation_CAMA.1_R.DMSO_S.DMSO_HSD17B_genes_logCPM_normalized_counts.txt", 
                 header=T)

br <- data.frame(t(br))
br$Iso <- c(rep("Sensitive", 3), rep("Resistant",3))

pdf("Interesting_Genes_1.pdf", height=4, width=9)
ggplot(br, aes(y=HSD17B1, x=Iso)) +
  geom_point(size=8, aes(color=Iso, alpha=0.85)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="grey10", size=1, alpha=0.85) +
  scale_color_manual(values = c("#1E88E5", "#004D40")) +
  theme_classic(base_size = 22) +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x="", y="log CPM", title="HSD17B1")  

ggplot(br, aes(y=HSD17B8, x=Iso)) +
  geom_point(size=8, aes(color=Iso, alpha=0.85)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="grey10", size=1, alpha=0.85) +
  scale_color_manual(values = c("#1E88E5", "#004D40")) +
  theme_classic(base_size = 22) +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x="", y="log CPM", title="HSD17B8")    
dev.off()

t.test(HSD17B1 ~ Iso, data = br)
t.test(HSD17B8 ~ Iso, data = br)
