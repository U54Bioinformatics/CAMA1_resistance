prefix=VG_CAMA1_D11_ALL
cat $prefix\_10x_gene_symbols.filtered.counts.*.mcherry_cell.txt > $prefix\_10x_gene_symbols.filtered.counts.mcherry_cell.txt
cat $prefix\_10x_gene_symbols.filtered.counts.*.mvenus_cell.txt > $prefix\_10x_gene_symbols.filtered.counts.mvenus_cell.txt
python scRNA_Seurat_05individual_from_individual_add_marker_genes.py --input $prefix\_10x_cell_metadata.UMAPcluster.txt --mcherr $prefix\_10x_gene_symbols.filtered.counts.mcherry_cell.txt --mvenus $prefix\_10x_gene_symbols.filtered.counts.mvenus_cell.txt
