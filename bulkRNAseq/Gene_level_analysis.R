library(data.table)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(stats)
library(Rcpp)
library(DT)


myTheme <- theme_bw() +
  theme(text = element_text(size=12, face = "bold"),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        strip.background = element_blank(),
        #axis.text.x.bottom = element_blank(),
        axis.title.x.bottom = element_text(size = 16, face = "bold"),
        strip.placement = "outside")

#function to calculate mean and sd for error bars
data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=T),
          sd = sd(x[[col]], na.rm=T),
          se = sd(x[[col]]/sqrt(length(x[[col]]))))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = paste0(varname,"_mean")))
    return(data_sum)

#Definitions
cellline <- "CAMA.1"
condition <- "DMSO"
sample <- "COH069"
type <- "bulk_RNA"
collection <- "C2"
experiment <- "riboR_DMSO.vs.S_DMSO2"
prefix <- paste(sample, type, sep = "_")

#Load logCPM values
m <- read.delim("facilitation_CAMA.1_R.DMSO_S.DMSO_HSD17B_genes_logCPM_normalized_counts.txt", header=T)

#extract gene of interest
m2 <- m[grepl("HSD17B", rownames(m)),]

write.table(m2, paste0("facilitation_CAMA.1_R.DMSO_S.DMSO_HSD17B_genes_", 
                      "logCPM_normalized_counts.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)
m2 <- read.table(file = paste0("facilitation_CAMA.1_R.DMSO_S.DMSO_HSD17B_genes_", 
                      "logCPM_normalized_counts.txt"),
            sep = "\t", header = T, quote = "")

#
g <- setDT(as.data.frame(m2), keep.rownames = T)
g2 <- g %>%
  dplyr::rename(Gene = rn) %>%
  pivot_longer(! Gene, names_to = "Cell", values_to = "logCPM")
n <- g2$Cell
n2 <- n %>%
  as.matrix() %>%
  as.data.frame() %>%
  separate(V1, into = c("CellLine", "R_S", "Treatment", "junk", "Replicate"), 
           sep = "_") %>%
  mutate(Replicate = ifelse(is.na(Replicate), junk, Replicate)) %>%
  mutate(junk = ifelse(junk == Replicate, Treatment, junk)) %>%
  mutate(Treatment = ifelse(Treatment == junk, R_S, Treatment)) %>% 
  mutate(R_S = ifelse(R_S == Treatment, "sensitive", R_S)) %>% 
  dplyr::select(-junk, Replicate) %>%
  mutate(state = ifelse(R_S == "sensitive", "Sensitive", "Resistant"))
g2$CellLine <- n2$CellLine
g2$R_S <- n2$R_S
g2$Treatment <- n2$Treatment
g2$State <- n2$state
g2$Replicate <- n2$Replicate
g2 <- g2 %>%
  unite("Treatment_State", Treatment:State, remove = F) %>% 
  mutate(T_S = case_when(startsWith(Treatment_State, "DMSO_R") ~ "DMSO_R",
                        startsWith(Treatment_State, "DMSO_S") ~ "DMSO_S"))


#SENTITIVE CO VS SENSITIVE MONO
my_comparisons <- list(c("Sensitive (coculture)",
                         "Sensitive (monoculture)"))


#
cc <- as.matrix(m2)
cc <- t(cc)
cc <- as.data.frame(cc)
T_S <- c(rep("DMSO_S", 3), rep("DMSO_R", 3))
cc$T_S <- T_S
tt <- cc

#wilcox test
set.seed(1111)
f <- compare_means(HSD17B1~T_S, tt, method = "wilcox.test")
f

#Filter for Gene
gene <- "HSD17B1"
g3 <- g2 %>%
  filter(Gene == gene)
g3

#plot
ggplot(g3, aes(x= factor(T_S), 
                #x = interaction(T_S,factor(Gene)), 
                y=logCPM,
              color = T_S,
              #group = interaction(Gene, T_S), 
              #shape = T_S
              )) + 
  #geom_point(size = 3) +
  geom_point(position = position_dodge(width = 0.5),
             size = 3) +
  scale_color_manual(values = c("#0000FF", "#00C638")) +
  xlab("") +
  #theme_bw() +
  coord_flip() +
  facet_wrap(~Gene, scales = "free") +
  #facet_grid(space = Gene) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red", size = .3) +
  theme(axis.text.x.bottom = element_text(angle = 0,
                                          vjust = 0.95, 
                                          hjust = 0.95,
                                          size = 9,
                                          face = "bold")) +
  myTheme
