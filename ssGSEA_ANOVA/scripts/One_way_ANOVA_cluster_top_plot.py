#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import re
import os
import time
import argparse
import glob

def usage():
    test="name"
    message='''
python One_way_ANOVA_cluster_top_plot.py --input FEL011_ANOVA_dir_cluster.test.AdjustP.csv

1. Select significant pathway by FDR < 0.05 and mean difference 0.1.
2. Output top 50 to generate violin plots
    '''
    print message


#File,Pathway,F_value,P_value,FDR,cluster0 mean,cluster1 mean,cluster2 mean,cluster1-cluster0,cluster2-cluster0,cluster2-cluster1
#Line2_ABBUD_LIF_SIGNALING_1_DN.ANOVA_data.test,ABBUD_LIF_SIGNALING_1_DN,4.587188,0.01063858,1,0.0912625878835897,0.0905046891699473,0.102111411183225
def read_pvalues(infile, fdr_cut, diff_cut, top_cut):
    selected   = []
    index_mean = []
    s_mean = re.compile(r'mean')
    header = ''
    outfile_pathway = re.sub(r'.AdjustP.txt', r'.AdjustP.fdr%s_diff%s.txt' %(fdr_cut, diff_cut), infile)
    outfile_list    = re.sub(r'.AdjustP.txt', r'.AdjustP.fdr%s_diff%s.top%s_pathway.txt' %(fdr_cut, diff_cut, top_cut), infile)
    ofile_pathway   = open(outfile_pathway, 'w')
    ofile_list      = open(outfile_list, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            if len(line) > 2 and not line.startswith(r'File'):
                fdr       = float(unit[4])
                index_min = sorted(index_mean)[0]
                index_max = sorted(index_mean)[-1]+1
                means     = unit[index_min:index_max]
                diff      = max(map(float, means)) - min(map(float, means))
                direction = 'NA'
                #print index_min, index_max
                if float(unit[index_max-1]) > float(unit[index_min]):
                    direction = "UP"
                else:
                    direction = "DOWN"
                unit.append(str(direction))
                unit.append(str(diff))
                if fdr < float(fdr_cut) and diff >= float(diff_cut):
                    selected.append(unit)
            elif len(line) > 2:
                header = '%s\tDirection\tDiff' %(line)
                for index, title in enumerate(unit):
                    if s_mean.search(title):
                        index_mean.append(index)
    top_pathway = []            
    selected_sorted = sorted(selected, key=lambda x: float(x[4]))
    linen = 0
    print >> ofile_pathway, header
    print >> ofile_list, header
    for pathway in selected_sorted:
        print >> ofile_pathway, '\t'.join(pathway)
        linen += 1
        if linen <= int(top_cut):
            print >> ofile_list, pathway[1]
            top_pathway.append(pathway[1])
    ofile_pathway.close()
    ofile_list.close()
    convert_txt2xls(outfile_pathway)
    convert_txt2xls(outfile_list)
    return top_pathway

def convert_txt2xls(infile):
    os.system('python ~/software/Project_scripts/Scripts_python/txt2xlsx.py --input %s' %(infile))
 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input')
    parser.add_argument('--fdr')
    parser.add_argument('--diff')
    parser.add_argument('--top')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.fdr:
        args.fdr = 0.05
    if not args.diff:
        args.fdr = 0.1
    if not args.top:
        args.top = 50

    convert_txt2xls(args.input)
    top_pathways = read_pvalues(args.input, args.fdr, args.diff, args.top)
   
    #plot top50 violin plot
    violin_dir   = re.sub(r'.test.AdjustP.txt', r'_violin', args.input)
    violin_dir   = '%s_%s_%s_%s' %(violin_dir, args.fdr, args.diff, args.top)
    violin_pdf   = '%s/pdf' %(violin_dir)
    data_dir     = re.sub(r'.test.AdjustP.txt', r'', args.input)
    if not os.path.exists(violin_dir):
        os.mkdir(violin_dir)
        os.mkdir(violin_pdf)

    ##cp top pathway files into violin folder
    for pathway in top_pathways:
        cmd = 'cp %s/L*%s.* %s/' %(data_dir, pathway, violin_dir)
        #print cmd
        os.system(cmd)

    ##run R to plot in violin folder
    violin_R_files = glob.glob(r'%s/*.violin.R' %(violin_dir))
    for script in violin_R_files:
        cmd = 'cat %s | R --slave > %s.log 2>&1' %(script, script)
        #print cmd 
        os.system(cmd)
        
    ##move pdf into violin folder
    for pathway in top_pathways:
        cmd = 'mv %s/L*%s.*.pdf %s/pdf/' %(data_dir, pathway, violin_dir)  
        #print cmd
        os.system(cmd)
    
    violin_pdf_files = glob.glob(r'%s/pdf/*.violin.pdf' %(violin_dir))
    for pdf in violin_pdf_files:
        png = re.sub(r'.pdf$', r'.png', pdf)
        cmd =  'convert -density 300 %s -quality 100 %s' %(pdf, png)
        os.system(cmd)

#folder=./ANOVA_dir_cluster_violin/
#prefix=ANOVA_dir_cluster_4clusters_wo3
#folder=$prefix\_violin
#mkdir $folder
#awk '{print "cp /home/jichen/Projects/Breast/scRNA/FELINE/FEL013P102/ssGSEA_ANOVA/ANOVA_dir_cluster/L\*"$1 ".* ./ANOVA_dir_cluster_4clusters_wo3_violin/"}' $prefix\_violin.pathway#.list > $prefix\_violin.cp.sh
#bash $prefix\_violin.cp.sh
#ls ./$prefix\_violin/*.violin.R | awk '{print "cat "$1" | R --slave"}' > $prefix\.violin.run.sh
#bash $prefix\.violin.run.sh
#mv ./ANOVA_dir_cluster/*.pdf $prefix\_violin/


if __name__ == '__main__':
    main()
