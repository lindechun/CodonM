#!/usr/bin/env python
# coding=utf-8
from __future__ import division

__author__ = "Lin Dechun"
__copyright__ = "Copyright 2016, BGI Research."
__credits__ = ["Lin Dechun"]
__version__ = "1.0.0"
__maintainer__ = "Lin Dechun"
__email__ = "lindechun@genomics.cn"

from argparse import ArgumentParser
import re,sys,os

def parseCommand():
    parser = ArgumentParser(description='The flow of Comparative analysis pipeline of Codon Usage',version='1.0.0')
    parser.add_argument('-i', action='store', dest='falist', help='a list file consists of individual Genes, based codon align (megacc)')
    parser.add_argument('-I', action='store', dest='info', help='sample info table')
    parser.add_argument('-t', action='store', dest='set_table', default=1, type=int, help='Genetic Code Table, default is 1')
    parser.add_argument('-j', action='store', dest='orderGenes',default='', help='The order of genes')

    parser.add_argument('-g', action='store', dest='groupColumns',help='the columns of info table, which is the groups of sample.')
    parser.add_argument('-G', action='store', dest='orderGroups',help='the order of -g, which is columns of info table. When -ENC_bar_plot or -CpG, you should set -G')

    parser.add_argument('-CA',action='store_true', dest='CA', default=True,help=' Whether Correspondence Analysis for RSCU matrix of all genes')
    parser.add_argument('-HC', action='store_true', dest='HC', default=True,help='Whether Hierarchical Cluster for RSCU matrix of all genes')
    parser.add_argument('-ATGC_bar_plot', action='store_true', dest='plot_ATGC_percentage', default=True,help='Whether plot ATGC percentage of each gene and groupColumns that you appoint')

    parser.add_argument('-ENC_bar_plot',action='store_true', dest='ENC_bar', default=True,help='Whether plot ATGC percentage of each gene and groupColumns that you appoint')
    parser.add_argument('-ENC_GC3_plot', action='store_true', dest='ENC_GC3', default=True,help='Whether plot ATGC percentage of each gene and groupColumns that you appoint')
    parser.add_argument('-CpG_bar',action='store_true', dest='CpG_bar', default=True,help='Whether plot ATGC percentage of each gene and groupColumns that you appoint')
    parser.add_argument('-CpG_GC_corr', action='store_true', dest='CpG_GC_corr', default=True,help='Whether correlation analysis between GC1s-3s and CpG frequencies')

    parser.add_argument('-o', action='store', dest='Opath',default='results',help='Path of Output File')
    parser.add_argument('-p', action='store', dest='prefix',default='Example', help='Output file prefix')

    return parser.parse_args()

def main():
    para=parseCommand()
    falist=para.falist
    infoTable=para.info
    set_table=para.set_table
    Opath=para.Opath
    prefix=para.prefix
    CA=para.CA
    HC=para.HC
    orderGenes=para.orderGenes
    groupColumns=para.groupColumns
    orderGroups=para.orderGroups
    ENC_bar=para.ENC_bar
    ENC_GC3=para.ENC_GC3
    plot_ATGC_percentage=para.plot_ATGC_percentage
    CpG_bar=para.CpG_bar
    CpG_GC_corr=para.CpG_GC_corr

    if not falist or not groupColumns or not orderGroups or not infoTable or not orderGenes:
        sys.exit("Please check -i, -I, -g, -G, -j")

    ## Step00
    print("# Step00. create results path")
    print('mkdir '+Opath)
    os.system('mkdir '+Opath)

    dirpath=os.path.dirname(os.path.realpath(__file__))

    ## Step01. BasicCodonUsage
    print("\n# Step01. BasicCodonUsage")
    print("mkdir -p "+Opath+"/01.BasicCodonUsage")
    os.system("mkdir -p "+Opath+"/01.BasicCodonUsage")

    for i in open(falist,'r'):
        i=i.strip()
        geneID=os.path.basename(i).split('.')[0]

        variantCalling="python {0}/BasicCodonUsage.py -i {1} -t {2} -o {3}".format(dirpath,i,set_table,Opath+'/'+'/01.BasicCodonUsage/'+geneID)
        print(variantCalling)
        os.system(variantCalling)

    ## Step02. RSCU pattern
    print("\n# Step02. RSCU pattern")
    print("mkdir -p "+Opath+"/02.RSCU_pattern")
    os.system("mkdir -p "+Opath+"/02.RSCU_pattern")

    RSCU_pattern="python {0}/RSCU.pattern.py -i {1} -I {2} -j '{3}' -g {4} -o {5} -p {6}".format(dirpath,Opath+'/01.BasicCodonUsage',infoTable,orderGenes,groupColumns,Opath+'/02.RSCU_pattern',prefix)
    if CA:
        RSCU_pattern += ' -CA'
    if HC:
        RSCU_pattern += ' -HC'
    print(RSCU_pattern)
    os.system(RSCU_pattern)

    ## Step03. ATGC content analysis
    print("\n# Step03. ATGC content analysis")
    print("mkdir -p "+Opath+"/03.ATGC_content")
    os.system("mkdir -p "+Opath+"/03.ATGC_content")

    ATGC_content="python {0}/ATGC.content.py -i {1} -I {2} -g {3} -o {4} -p {5}".format(dirpath,falist,infoTable,groupColumns,Opath+'/03.ATGC_content',prefix)

    if plot_ATGC_percentage:
        ATGC_content += ' -ATGC_bar_plot -G '+"'"+orderGroups+"'"

    print(ATGC_content)
    os.system(ATGC_content)

    ## Step04. ENC and ENC-GC3s analysis
    print("\n# Step04. ENC and ENC-GC3s analysis")
    print("mkdir -p "+Opath+"/results/04.ENC_GC3s")
    os.system("mkdir -p "+Opath+"/results/04.ENC_GC3s")

    ENC_GC3_command="python {0}/ENC_GC3s.py -i {1} -I {2} -g {3} -o {4} -p {5} -G '{6}'".format(dirpath,falist,infoTable,groupColumns,Opath+'/04.ENC_GC3s',prefix,orderGroups)

    if ENC_GC3:
        ENC_GC3_command += ' -ENC_GC3_plot -ATGC '+Opath+'/03.ATGC_content/'+prefix+'.ATGC_Content.overAll.txt'+' -j '+"'"+orderGenes+"'"
    if ENC_bar:
        ENC_GC3_command += ' -ENC_bar_plot'

    print(ENC_GC3_command)
    os.system(ENC_GC3_command)

    ## Step05. Dinucleotide. particularly,CpG
    print("\n# Step05. Dinucleotide. particularly,CpG")
    print("mkdir -p "+Opath+"/05.Dinucleotide")
    os.system("mkdir -p "+Opath+"/05.Dinucleotide")

    dinuc_command="python {0}/Dinucleotide.py -i {1} -I {2} -g {3} -o {4} -p {5}".format(dirpath,falist,infoTable,groupColumns,Opath+'/05.Dinucleotide',prefix)

    if CpG_bar:
        dinuc_command += " -CpG_bar -G "+"'"+orderGroups+"'"
    if CpG_GC_corr:
        dinuc_command +=" -CpG_GC_corr -ATGC "+Opath+'/03.ATGC_content/'+prefix+'.ATGC_Content.overAll.txt'

    print(dinuc_command)
    os.system(dinuc_command)

if __name__  == "__main__":
    main()
