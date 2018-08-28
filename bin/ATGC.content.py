#!/usr/bin/env python
# coding=utf-8

__author__ = "Lin Dechun"
__copyright__ = "Copyright 2017, BGI Research."
__credits__ = ["Lin Dechun"]
__version__ = "0.0.1"
__maintainer__ = "Lin Dechun"
__email__ = "lindechun@genomics.cn"

from argparse import ArgumentParser
import os,sys
import pandas as pd
from glob import glob
from Bio import SeqIO
from collections import defaultdict

def parseCommand():
    parser = ArgumentParser(description='This script use to convert nucleotide to amino acid',version='1.0.0')
    parser.add_argument('-i', action='store', dest='falist', help='a list consists of individual Genes fa, bsed codon align (megacc)')
    parser.add_argument('-ATGC_bar_plot', action='store_true', dest='ATGC_bar_plot', default=False,help='Whether plot ATGC percentage of each gene and groupColumns that you appoint')
    parser.add_argument('-I', action='store', dest='infoTable', help=' sample info table. When you set -R, you should set -I')
    parser.add_argument('-g', action='store', dest='groupColumns',help='the columns of info table, which is the groups of sample. When -R, you should set -g')
    parser.add_argument('-G', action='store', dest='orderGroups',help='the order of -g, which is columns of info table. When -R, you should set -G')
    parser.add_argument('-o', action='store', dest='Opath',default='./', help='Path of Output File')
    parser.add_argument('-p', action='store', dest='prefix', help='prefix of Output File')

    return parser.parse_args()

def ATGC_cal_overall(dnaSequence):
    aCount = 0
    cCount = 0
    tCount = 0
    gCount = 0

    for c in dnaSequence:
        if c.lower() == 'a':
            aCount = aCount + 1
        elif c.lower() == 'c':
            cCount = cCount + 1
        elif c.lower() == 't':
            tCount = tCount + 1
        elif c.lower() == 'g':
            gCount = gCount + 1

    allCount=aCount+cCount+tCount+gCount
    return map(lambda x:round(x/float(allCount),3),[aCount,cCount,tCount,gCount])

def ATGC_cal_anylocus_codon(dnaSequence,n):
    aCount = 0
    cCount = 0
    tCount = 0
    gCount = 0

    for i in range(n-1,len(dnaSequence),3):
        c=dnaSequence[i]
        if c.lower() == 'a':
            aCount = aCount + 1
        elif c.lower() == 'c':
            cCount = cCount + 1
        elif c.lower() == 't':
            tCount = tCount + 1
        elif c.lower() == 'g':
            gCount = gCount + 1

    allCount=aCount+cCount+tCount+gCount
    return map(lambda x:round(x/float(allCount),3),[aCount,cCount,tCount,gCount])


def everySample_ATGC(seq):
    ATGC=defaultdict(dict)

    for fa in SeqIO.parse(seq,'fasta'):

        ### overall ATGC content
        t1,t2,t3,t4=ATGC_cal_overall(fa.seq)
        ATGC[fa.id]['A']=t1
        ATGC[fa.id]['C']=t2
        ATGC[fa.id]['T']=t3
        ATGC[fa.id]['G']=t4

        ### ATGC content of the 1rd of codon
        t1,t2,t3,t4=ATGC_cal_anylocus_codon(fa.seq,1)
        ATGC[fa.id]['A1']=t1
        ATGC[fa.id]['C1']=t2
        ATGC[fa.id]['T1']=t3
        ATGC[fa.id]['G1']=t4
        ATGC[fa.id]['GC1']=round(t2+t4,3)

        ### ATGC content of the 2rd of codon
        t1,t2,t3,t4=ATGC_cal_anylocus_codon(fa.seq,2)
        ATGC[fa.id]['A2']=t1
        ATGC[fa.id]['C2']=t2
        ATGC[fa.id]['T2']=t3
        ATGC[fa.id]['G2']=t4
        ATGC[fa.id]['GC2']=round(t2+t4,3)

        ### ATGC content of the 3rd of codon
        t1,t2,t3,t4=ATGC_cal_anylocus_codon(fa.seq,3)
        ATGC[fa.id]['A3']=t1
        ATGC[fa.id]['C3']=t2
        ATGC[fa.id]['T3']=t3
        ATGC[fa.id]['G3']=t4
        ATGC[fa.id]['GC3']=round(t2+t4,3)
    return ATGC

def plot_ATGC_content(dat,infoTable,groupColumns,Opath,prefix,orderGroups):
    info=pd.read_table(infoTable,sep="\t",index_col=0)

    dat['Groups']=info.loc[dat['ID'],groupColumns].values
    del dat['ID']

    ## ATGC mean(1-3rd)
    dat1=dat.loc[dat.Base.isin(['A','T','G','C']),]
    dat1=pd.melt(dat1,id_vars=['Base','Groups'])
    dat1=dat1.rename(columns={"variable":'Gene','value':'Percentage'})
    dat1=pd.DataFrame(dat1.groupby(['Base','Groups','Gene'])['Percentage'].mean())
    dat1.to_csv(Opath+'/'+prefix+'.ATGC_Content.1-3rd.Codon.withinGroups.txt',sep="\t")

    ## ATGC mean(3rd)
    dat3=dat.loc[dat.Base.isin(['A3','T3','G3','C3']),]
    dat3=pd.melt(dat3,id_vars=['Base','Groups'])
    dat3=dat3.rename(columns={"variable":'Gene','value':'Percentage'})
    dat3=pd.DataFrame(dat3.groupby(['Base','Groups','Gene'])['Percentage'].mean())
    dat3.to_csv(Opath+'/'+prefix+'.ATGC_Content.3rd.Codon.withinGroups.txt',sep="\t")

    os.system("Rscript "+os.path.dirname(os.path.realpath(__file__))+"/ATGC.percentage.plot.R "+Opath+"/"+prefix+".ATGC_Content.1-3rd.Codon.withinGroups.txt "+Opath+"/"+prefix+".ATGC_Content.3rd.Codon.withinGroups.txt "+orderGroups+' '+Opath+'/'+prefix)

def main():
    para=parseCommand()
    falist=para.falist
    ATGC_bar_plot=para.ATGC_bar_plot
    infoTable=para.infoTable
    groupColumns=para.groupColumns
    orderGroups=para.orderGroups
    Opath=para.Opath
    prefix=para.prefix

    all_gene_ATGC={}
    os.system("mkdir -p "+Opath+'/eachGene/')

    if not falist:
        sys.exit("Please set -i")

    for i in open(falist,'r'):
        i=i.strip()

        geneName=os.path.basename(i).split('.')[0]

        ATGC=everySample_ATGC(i)

        ATGC_percent=pd.DataFrame(ATGC).T
        ATGC_percent.index.name='ID'
        ATGC_percent=ATGC_percent.loc[:,['A','T','G','C','A1','T1','G1','C1','GC1','A2','T2','G2','C2','GC2','A3','T3','G3','C3','GC3']]
        ATGC_percent.to_csv(Opath+'/eachGene/'+geneName+'.ATGC.content.txt',sep="\t")

        all_gene_ATGC[geneName]=ATGC

    dat=pd.Panel.from_dict(all_gene_ATGC)
    dat=dat.to_frame()
    dat.index.names=['Base','ID']
    dat.to_csv(Opath+'/'+prefix+'.ATGC_Content.overAll.txt',sep="\t")
    dat=dat.reset_index()

    if ATGC_bar_plot:
        plot_ATGC_content(dat,infoTable,groupColumns,Opath,prefix,orderGroups)


if __name__ == '__main__':
    main()
