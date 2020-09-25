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
python One_way_ANOVA.py --input test_score.txt

    '''
    print message

#Gene Set        AB11A_AAACGAACACGGCACT
#ABBUD_LIF_SIGNALING_1_DN        0.0747999808935552
def Prepare_ANOVA(infile, skip, cell_meta, folder):
    test_files = []
    line_num = 0
    header   = []
    ofile = open ('%s.run.sh' %(folder), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            line_num += 1
            if line_num == 1:
                for name in unit[1:]:
                    #name_split = re.split(r'_', name)
                    #name_split[0] = re.sub(r'AB11i', r'AB11I', name_split[0])
                    name = re.sub(r'AB11i', r'AB11I', name)
                    if cell_meta.has_key(name):
                        header.append(cell_meta[name])
                    else:
                        print 'Cluster not found'
                        header.append(['NA', 'NA'])
            else:
                pathway = unit[0]
                test_file, commandline = writefile(folder, pathway, line_num, header, unit[1:], skip)
                test_files.append(test_file)
                print >> ofile, '\n'.join(commandline)
    return test_files


def writefile(outdir, pathway, line_num, header, scores, skip):
    commandline = []
    test_file = ''
    #prepare table
    anova_data_txt = '%s/Line%s_%s.ANOVA_data.txt' %(outdir, line_num, pathway)
    ofile = open(anova_data_txt, 'w')
    print >> ofile, 'Pathway\tSample\tCluster\tScore' 
    for i in range(len(header)):
        if not header[i][0] == 'NA':
            print >> ofile, '%s\t%s\t%s\t%s' %(pathway, header[i][0], header[i][1], scores[i])
    ofile.close()

    #one way anova test
    anova_data_r    = '%s/Line%s_%s.ANOVA_data.R' %(outdir, line_num, pathway)
    anova_data_fp   = '%s/Line%s_%s.ANOVA_data.F_P_values.txt' %(outdir, line_num, pathway)
    anova_data_mean = '%s/Line%s_%s.ANOVA_data.mean_values.txt' %(outdir, line_num, pathway)
    ofile = open(anova_data_r, 'w')
    cmd='''
#https://stats.idre.ucla.edu/r/faq/how-can-i-do-post-hoc-pairwise-comparisons-in-r/
#read in data
x <- read.table("%s", header = T, sep="\\t")
#one way ANOVA
a1 <- aov(x$Score ~ x$Sample)
summary(a1)
write(x=c(summary(a1)[[1]][["F value"]][1], summary(a1)[[1]][["Pr(>F)"]][1]), file="%s")
#mean
x_mean <- aggregate(x$Score, list(x$Sample), mean)
write(x=paste(x_mean$Group.1, x_mean$x), file="%s")
#pairwise test
TukeyHSD(a1)
''' %(os.path.abspath(anova_data_txt), os.path.abspath(anova_data_fp), os.path.abspath(anova_data_mean))
    print >> ofile, cmd
    ofile.close()

    anova_data_test = '%s/Line%s_%s.ANOVA_data.test' %(outdir, line_num, pathway)
    anova_data_log  = '%s/Line%s_%s.ANOVA_data.test.log' %(outdir, line_num, pathway)
    commandtemp     = 'cat %s | R --slave > %s 2> %s' %(os.path.abspath(anova_data_r), os.path.abspath(anova_data_test), os.path.abspath(anova_data_log))
    commandline.append(commandtemp)
    if not skip:
        os.system(commandtemp)
    test_file=os.path.abspath(anova_data_test)
    

    #violin plots
    anova_data_violin_r='%s/Line%s_%s.ANOVA_data.violin.R' %(outdir, line_num, pathway)
    anova_data_violin_pdf='%s/Line%s_%s.ANOVA_data.violin.pdf' %(outdir, line_num, pathway)
    anova_data_violin_png='%s/Line%s_%s.ANOVA_data.violin.png' %(outdir, line_num, pathway)
    cmd='''
pdf("%s", width=8, height=6)
x <- read.table("%s", header = T, sep="\\t")
library(ggplot2)
library(stringr)
x$Cluster <- str_replace_all(x$Cluster, "cluster", "")
names(x) <- c("Pathway","Sample","Subclone","Score")
x$Subclone <- as.factor(x$Subclone)
#with boxplot
#p <- ggplot(x, aes(x=Sample, y=Score, fill=Sample)) +
#geom_violin(trim=TRUE)
#p + geom_boxplot(width=0.1) + theme_minimal() + labs(title="") + theme(plot.title = element_text(hjust = 0.5))

#with jitter
#p <- ggplot(x, aes(x=Sample, y=Score, fill=Sample)) +
#geom_violin(trim=TRUE)
#p + geom_jitter(shape=16, position=position_jitter(0.2))

#lance code
#colors <- c("#E86A10","#56A4AA", "#3A78C4","#F1AB00", "#2F992D", "#70738B")
colors <- c("#E86A10")
fontsize=22
sample_level <- as.vector(sort(unique(factor(x$Sample)), decreasing= T))
p <- ggplot(x, aes(x=factor(Sample, level=sample_level), y=Score)) + labs(title="%s", x="Samples", y = "GSEA_Score")+ geom_violin() + 
geom_jitter(shape=16, position=position_jitter(0.1), size=0.5, aes(color=Subclone)) + scale_color_manual(values=colors) + theme_classic() +
theme_bw()+theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(colour='white'),
                 panel.background = element_blank()) +
theme(axis.text=element_text(size=fontsize, color='black'),
      axis.text.x =element_text(size=fontsize, angle = 25, color='black', hjust = 1),
      axis.title.x =element_blank(),
      axis.title.y =element_text(size=fontsize),
      strip.text = element_text(face = 'bold', size=12)) +
theme(legend.text = element_text(size=fontsize),
      legend.title = element_blank(),
      legend.position = c(0.2, 0.9)) +
theme(plot.title=element_text(size=fontsize,face="bold"), axis.text=element_text(size=fontsize, face="bold")) +
theme(plot.margin =  margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")) +
guides(colour = guide_legend(override.aes = list(size=4)))
#plot
p
dev.off()
''' %(os.path.abspath(anova_data_violin_pdf), os.path.abspath(anova_data_txt), pathway)

    ofile = open(anova_data_violin_r, 'w')
    print >> ofile, cmd
    ofile.close()
    commandtemp1 = 'cat %s | R --slave > pdf.log 2> pdf.log2' %(os.path.abspath(anova_data_violin_r))
    commandtemp2 = 'convert -density 300 %s -quality 100 %s' %(os.path.abspath(anova_data_violin_pdf), os.path.abspath(anova_data_violin_png))
    #commandline.append(commandtemp)
    if not skip:
        os.system(commandtemp1)
        os.system(commandtemp2)
   
    return test_file, commandline

#Line10_ABE_VEGFA_TARGETS_30MIN.ANOVA_data.test
#              Df Sum Sq Mean Sq F value Pr(>F)
#x$Sample       4  0.719 0.17972   78.92 <2e-16 ***
#Residuals   8351 19.018 0.00228

#
#$`x$Sample`
#                     diff          lwr           upr     p adj
#AB11C-AB11A  0.0294157534  0.025267974  3.356353e-02 0.0000000
#AB11F-AB11A  0.0295497537  0.025792737  3.330677e-02 0.0000000
#AB11I-AB11A  0.0251840880  0.021291445  2.907673e-02 0.0000000
#AB11J-AB11A  0.0312143082  0.027146263  3.528235e-02 0.0000000
#AB11F-AB11C  0.0001340003 -0.003873189  4.141190e-03 0.9999844
#AB11I-AB11C -0.0042316654 -0.008366283 -9.704786e-05 0.0418391
#AB11J-AB11C  0.0017985548 -0.002501607  6.098716e-03 0.7846058
#AB11I-AB11F -0.0043656657 -0.008108146 -6.231851e-04 0.0127368
#AB11J-AB11F  0.0016645545 -0.002260046  5.589155e-03 0.7757485
#AB11J-AB11I  0.0060302202  0.001975595  1.008485e-02 0.0004798


#mean
#AB11A 0.0852042788565998
#AB11C 0.0748986828801092
def readanova(test):
    #test_fp   = re.sub(r'ANOVA_data.test', r'ANOVA_data.F_P_values.txt', test)
    #pathway mean
    test_mean = re.sub(r'ANOVA_data.test', r'ANOVA_data.mean_values.txt', test)
    pathway_title, pathway_mean = read_mean(test_mean)
    #pairwise p value
    pairwise_title, pairwise_value  = read_test(test)

    #f and p values
    test_fp   = re.sub(r'ANOVA_data.test', r'ANOVA_data.F_P_values.txt', test)
    name    = os.path.split(test)[1]
    searchpat = re.compile(r'Line\d+_(.*?).ANOVA_data.test')
    pathway   = name
    if searchpat.search(name):
        pathway = searchpat.search(name).groups(0)[0]
    line_num = 0
    anova_stat = ''
    with open (test_fp, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\s+',line)
            line_num += 1
            if line_num == 1:
                fvalue = unit[0]
                pvalue = unit[1]
                #print >> ofile, '%s\t%s\t%s\t%s' %(name, pathway, fvalue, pvalue)
                anova_stat = '%s,%s,%s,%s,%s,%s' %(name, pathway, fvalue, pvalue, ','.join(pathway_mean), ','.join(pairwise_value))
                #print anova_stat
           
    anova_title = 'File,Pathway,F_value,P_value,%s,%s' %(','.join(pathway_title), ','.join(pairwise_title)) 
    return anova_title, anova_stat


#AB11A 0.0120316007731628
def read_mean(infile):
    #data = defaultdict(lambda : str())
    title = []
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r' ',line)
            #print line
            #print unit[1]
            title.append('%s mean' %(unit[0]))
            data.append(unit[1])
    return title, data

#$`x$Sample`
#                     diff          lwr           upr     p adj
#AB11C-AB11A  0.0294157534  0.025267974  3.356353e-02 0.0000000
#AB11F-AB11A  0.0295497537  0.025792737  3.330677e-02 0.0000000
#AB11I-AB11A  0.0251840880  0.021291445  2.907673e-02 0.0000000
#AB11J-AB11A  0.0312143082  0.027146263  3.528235e-02 0.0000000
#AB11F-AB11C  0.0001340003 -0.003873189  4.141190e-03 0.9999844
#AB11I-AB11C -0.0042316654 -0.008366283 -9.704786e-05 0.0418391
#AB11J-AB11C  0.0017985548 -0.002501607  6.098716e-03 0.7846058
#AB11I-AB11F -0.0043656657 -0.008108146 -6.231851e-04 0.0127368
#AB11J-AB11F  0.0016645545 -0.002260046  5.589155e-03 0.7757485
#AB11J-AB11I  0.0060302202  0.001975595  1.008485e-02 0.0004798
def read_test(infile):
    title = []
    value = []
    flag = 0
    line_num = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\s+',line)
            if line.startswith(r'$`x$Sample'):
                flag = 1
            if flag == 1:
                line_num += 1
                if line_num >= 3 and len(line) > 10:
                    print unit[0], unit[3]
                    title.append(unit[0])
                    value.append(unit[4])
    return title, value 

def Sum_ANOVA(test_files, folder):
    anova_stat = []
    test_title = ''
    for test in test_files:
        #print test
        test_title, test_value = readanova(test)
        anova_stat.append(test_value)
    #anova_stat_corrected = multitest_correction(anova_stat)
    ofile = open('%s.test.txt' %(folder), 'w')
    print >> ofile, test_title
    for stat in anova_stat:
        print >> ofile, stat
        print stat
    ofile.close()

    R_code='''
library(data.table)
x <- fread("%s.test.txt", sep=",")
FDR <- p.adjust(x$P_value)
x_new <- cbind(x[,1:4], FDR, x[,5:ncol(x)])
fwrite(x_new, file="%s.test.AdjustP.txt", quote = FALSE, sep = "\t", row.names = FALSE)
''' %(os.path.split(folder)[1], os.path.split(folder)[1])

    R_script='%s.test.AdjustP.R' %(folder)
    ofile = open(R_script, 'w')
    print >> ofile, R_code
    ofile.close

    os.system('cat %s | R --slave > %s.log 2>&1' %(R_script, R_script))


#Cell.ID Sample  seurat_clusters Marker_genes    Marker_groups
#GV014_AAACCCAAGAAGCCAC  GV014   6       mcherry GV014_mcherry
#GV014_AAACCCAAGAGGACTC  GV014   2       mvenus  GV014_mvenus
#GV014_AAACCCAAGGGACACT  GV014   0       nomarker        GV014_nomarker
def read_cluster_file(infile, analysis):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if not line.startswith(r'Cell.ID'):
                unit = re.split(r'\t',line)
                if analysis == '1':
                    #comparing GV014_mvenus vs. GV013_mvenus
                    cell      = unit[0]
                    sample    = unit[4]
                    cluster   = unit[2]
                    if sample in ["GV014_mvenus", "GV013_mvenus"]:
                        data[cell]= [sample, cluster]
                elif analysis == '2':
                    #comparing GV013_mvenus vs. GV015_mcherry
                    cell      = unit[0]
                    sample    = unit[4]
                    cluster   = unit[2]
                    if sample in ["GV013_mvenus", "GV015_mcherry"]:
                        data[cell]= [sample, cluster]
                elif analysis == '3':
                    #comparing mvenus vs. mcherry
                    cell      = unit[0]
                    sample    = unit[3]
                    cluster   = unit[2]
                    if sample in ["mvenus", "mcherry"]:
                        data[cell]= [sample, cluster]
                elif analysis == '4':
                    #comparing GV014_mcherry vs. GV015_mcherry
                    cell      = unit[0]
                    sample    = unit[4]
                    cluster   = unit[2]
                    if sample in ["GV014_mcherry", "GV015_mcherry"]:
                        data[cell]= [sample, cluster]
                elif analysis == '5': 
                    #comparing GV014_mcherry vs. GV014_mvenus
                    cell      = unit[0]
                    sample    = unit[4]
                    cluster   = unit[2]
                    if sample in ["GV014_mcherry", "GV014_mvenus"]:
                        data[cell]= [sample, cluster]

    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--cluster')
    parser.add_argument('--skip')
    parser.add_argument('--summary')
    parser.add_argument('--analysis')
    parser.add_argument('--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0 and len(args.cluster) > 0
    except:
        usage()
        sys.exit(2)


    if not args.analysis:
        args.analysis   = "1"
    if not args.output:
        args.output = "ANOVA_dir_normalcell"
    if not args.summary:
        args.summary = 0
    if not args.skip:
        args.skip = 1

    args.output = '%s_%s' %(args.output, args.analysis)

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    
 
    #read cell -> sample/timepoint or cell -> cluster info into a dict
    cell_meta  = read_cluster_file(args.cluster, args.analysis)
    #prepare ANOVA files
    test_files = Prepare_ANOVA(args.input, args.skip, cell_meta, args.output)
    #summary ANOVA test into a table
    if args.summary:
        Sum_ANOVA(test_files, args.output)


if __name__ == '__main__':
    main()
