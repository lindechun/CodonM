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
    parser = ArgumentParser(description='Calculate ENC and ENC-gc3S combined with GC3s',version='1.0.0')
    parser.add_argument('-i', action='store', dest='falist', help='a list consists of individual Genes fa, based on codon align (megacc)')
    parser.add_argument('-ENC_bar_plot',action='store_true', dest='ENC_bar', default=False,help='Whether plot ATGC percentage of each gene and groupColumns that you appoint')
    parser.add_argument('-ENC_GC3_plot', action='store_true', dest='ENC_GC3', default=False,help='Whether plot ATGC percentage of each gene and groupColumns that you appoint')
    parser.add_argument('-I', action='store', dest='infoTable', help=' sample info table. When you set -ENC_bar_plot or -ENC_GC3_plot, you should set -I')
    parser.add_argument('-g', action='store', dest='groupColumns',help='the columns of info table, which is the groups of sample. When -ENC_bar_plot or -ENC_GC3_plot, you should set -g')
    parser.add_argument('-G', action='store', dest='orderGroups',help='the order of -g, which is columns of info table. When -ENC_bar_plot or -ENC_GC3_plot, you should set -G')
    parser.add_argument('-j', action='store', dest='orderGenes',help='the order of -g, which is columns of info table. When -ENC_GC3_plot, you should set -j')
    parser.add_argument('-ATGC', action='store', dest='ATGC',help='ATGC content overall results. When -ENC_GC3_plot, you should set -gc3s')
    parser.add_argument('-o', action='store', dest='Opath',default='./', help='Path of Output File')
    parser.add_argument('-p', action='store', dest='prefix', help='prefix of Output File')

    return parser.parse_args()


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
    ENC_bar=para.ENC_bar
    ENC_GC3=para.ENC_GC3

    infoTable=para.infoTable
    groupColumns=para.groupColumns
    orderGroups=para.orderGroups
    orderGenes=para.orderGenes
    ATGC=para.ATGC
    Opath=para.Opath
    prefix=para.prefix

    if not falist or not infoTable:
        sys.exit("please set -i and -I.")

    os.system("mkdir -p "+Opath+'/eachGene/')

    ### calculate ENC
    for i in open(falist,'r'):
        i=i.strip()
        geneName=os.path.basename(i).split('.')[0]

        os.system('codonw -enc '+ i+' '+Opath+'/eachGene/'+geneName+'.ENC.txt -silent -nomenu -nowarn -noblk')

    if ENC_bar or ENC_GC3:
        ### gather ENCs
        ENC_list=[]
        for i in glob(Opath+'/eachGene/*.ENC.txt'):
            gene_name=os.path.basename(i)[:-8]
            ENC=pd.read_table(i,sep="\t")
            ENC.columns=['ID_temp','ENC','AA']
            del ENC['AA']
            ENC.ID_temp=ENC.ID_temp.apply(lambda x:x.strip())
            ENC.insert(0,'gene',[gene_name]*ENC.shape[0])
            ENC_list.append(ENC)

        ENCs=pd.concat(ENC_list)

        ENCs.insert(1,'ID',ENCs['gene']+'.'+ENCs['ID_temp'])
        info=pd.read_table(infoTable,sep="\t",index_col=0)
        ENCs['Groups']=info.loc[ENCs.ID_temp,groupColumns].values
        del ENCs['ID_temp']
        ENCs.set_index("ID",inplace=True)

        ENCs.to_csv(Opath+'/'+prefix+".allGene.ENC_GC3s.txt",sep="\t")

        if ENC_bar:
            os.system("Rscript "+os.path.dirname(os.path.realpath(__file__))+"/ENC.barPlot.R "+Opath+"/"+prefix+".allGene.ENC_GC3s.txt "+orderGroups+' '+Opath+'/'+prefix)
        

        if ENC_GC3:
            ## gather GC3s and merge to ENCs
            ATGC_overAll=pd.read_table(ATGC,sep='\t')

            ATGC_gc3s=pd.melt(ATGC_overAll,id_vars=['Base','ID'])
            ATGC_gc3s=ATGC_gc3s.loc[ATGC_gc3s.Base=='GC3']
            ATGC_gc3s['NEW_ID']=ATGC_gc3s.variable+'.'+ATGC_gc3s.ID
            ATGC_gc3s.set_index('NEW_ID',inplace=True)
            ENCs.insert(2,'GC3s',ATGC_gc3s.loc[ENCs.index,'value'].values)

            ENCs.to_csv(Opath+'/'+prefix+".allGene.ENC_GC3s.txt",sep="\t")

            ## individual Genes
            os.system("Rscript "+os.path.dirname(os.path.realpath(__file__))+"/ENC_GC3s.individualGenes.R "+Opath+"/"+prefix+".allGene.ENC_GC3s.txt "+orderGroups+' '+orderGenes+' '+Opath+'/'+prefix)

            ## Group By Group
            os.system("Rscript "+os.path.dirname(os.path.realpath(__file__))+"/ENC_GC3s.GroupByGroup.R "+Opath+"/"+prefix+".allGene.ENC_GC3s.txt "+orderGroups+' '+Opath+'/'+prefix)

if __name__ == '__main__':
    main()
