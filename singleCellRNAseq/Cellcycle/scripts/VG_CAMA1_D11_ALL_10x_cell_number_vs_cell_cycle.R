library(data.table)
library(ggplot2)
library(dplyr)
data <- fread("VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.marker_genes.txt", header=T, sep="\t")
data_filtered <- dplyr::filter(data, Marker_genes != "nomarker")
df <- as.data.frame(table(data_filtered$Marker_groups, data_filtered$Phase))
names(df) <- c("Cell", "Cell_cycle", "Cell_Number")

test_table <- dcast(data=df, formula = Cell~Cell_cycle)
rownames(test_table) <- test_table$Cell
test_table$Cell <- NULL
test_table
print("Testing cell cycle between mono-resistant vs. co-resistant")
test_table[1:2, 1:3]
chisq.test(test_table[1:2, 1:3])
print("Testing cell cycle between mono-sensitive vs. co-sensitive")
test_table[3:4, 1:3]
chisq.test(test_table[3:4, 1:3])

fwrite(test_table, "VG_CAMA1_D11_ALL_10x_cell_number_vs_cell_cycle.txt", sep="\t", quote=F, row.names=T)
