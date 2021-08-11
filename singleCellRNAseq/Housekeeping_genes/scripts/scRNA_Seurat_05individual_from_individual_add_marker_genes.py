#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse

def usage():
    test="name"
    message='''
python Add_cell_type_to_meta.py --cluster_anno FEL012016_10x_SingleR_cluster_table_Blueprint.cluster_anno.txt --meta_cluster FEL012016_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.txt


    '''
    print message


#Cell.ID orig.ident      nCount_RNA	nFeature_RNA Encode_main_cluster     seurat_clusters 
#FEL013_FEL013_E_TTCACCGAGGTATCTC        FEL013  2267    997 Encode_main_cluster     seurat_clusters
def read_meta_with_cluster(infile, mcherry, mvenus):
    clusters     = defaultdict(lambda : str())
    sample_anno  = defaultdict(lambda : defaultdict(lambda : int()))
    linen     = 0
    output = re.sub(r'.txt', r'.marker_genes.txt', infile)
    ofile  = open(output, 'w')
    cluster_index = 0
    anno_index    = 0
    sample_index  = 0 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            if len(line) > 2 and line.startswith(r'Cell'):
                for index, value in enumerate(unit):
                    if value == "seurat_clusters":
                        cluster_index = index
                    elif value == "Encode_main_cluster":
                        anno_index    = index
                    elif value == "Sample":
                        sample_index = index
                unit.append('Marker_genes')
                unit.append('Marker_groups')
                unit.append('Marker_genes1')
                unit.append('Marker_groups1')
                print >> ofile, '\t'.join(unit)
            elif len(line) > 2:
                #print line
                #print cluster_index, unit[cluster_index]
                #unit[anno_index] = cluster_anno[unit[cluster_index]]
                #infercnv_anno    = unit[sample_index]
                #if not unit[anno_index] == 'Epithelial cells':
                #    infercnv_anno = 'Immune'
                marker_gene   = 'nomarker'
                marker_gene1  = 'nomarker'
                if mcherry.has_key(unit[0]):
                    marker_gene  = 'mcherry'
                    marker_gene1 = 'Resistant'
                elif mvenus.has_key(unit[0]):
                    marker_gene = 'mvenus'
                    marker_gene1 = 'Sensitive'
                marker_group  = '%s_%s' %(unit[sample_index], marker_gene)
                marker_group1 = marker_group
                if unit[sample_index] == 'GV013':
                    marker_group1 = 'Sensitive (monoculture)'
                elif unit[sample_index] == 'GV015':
                    marker_group1 = 'Resistant (monoculture)'
                elif unit[sample_index] == 'GV014':
                    if marker_gene1 == 'Resistant':
                        marker_group1 = 'Resistant (coculture)'
                    elif marker_gene1 == 'Sensitive':
                        marker_group1 = 'Sensitive (coculture)'          
                
                unit.append(marker_gene1)
                unit.append(marker_group1)
                unit.append(marker_gene)
                unit.append(marker_group)
                print >> ofile, '\t'.join(unit)
    ofile.close()


#GV014_AAACCCAAGAAGCCAC
#GV014_AAACCCATCGAATGCT
def read_cells(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t',line)
                data[unit[0]] = 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--mcherry')
    parser.add_argument('--mvenus')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.mcherry) > 0 and len(args.mvenus) > 0
    except:
        usage()
        sys.exit(2)

    cell_mcherry = read_cells(args.mcherry)
    cell_mvenus  = read_cells(args.mvenus)
    read_meta_with_cluster(args.input, cell_mcherry, cell_mvenus)
    

if __name__ == '__main__':
    main()

