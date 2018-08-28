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
from pandas import ExcelWriter
import xlsxwriter
from scipy.stats import spearmanr

def parseCommand():
    parser = ArgumentParser(description='This script use to convert nucleotide to amino acid',version='1.0.0')
    parser.add_argument('-i', action='store', dest='falist', help='a list consists of individual Genes fa, based on codon align (megacc)')
    parser.add_argument('-CpG_bar',action='store_true', dest='CpG_bar', default=False,help='Whether plot CpG bar of each gene and groupColumns that you appoint')
    parser.add_argument('-CpG_GC_corr', action='store_true', dest='CpG_GC_corr', default=False,help='Whether correlation analysis between GC1s-3s and CpG frequencies')
    parser.add_argument('-ATGC', action='store', dest='ATGC',help='ATGC content overall results. When -CpG_GC_corr, you should set -gcContent')
    parser.add_argument('-I', action='store', dest='infoTable', help=' sample info table. When you set -ENC_bar_plot or -ENC_GC3_plot, you should set -I')
    parser.add_argument('-g', action='store', dest='groupColumns',help='the columns of info table, which is the groups of sample.')
    parser.add_argument('-G', action='store', dest='orderGroups',help='the order of -g, which is columns of info table. When -ENC_bar_plot or -CpG, you should set -G')
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
    CpG_bar=para.CpG_bar
    CpG_GC_corr=para.CpG_GC_corr
    ATGC=para.ATGC
    infoTable=para.infoTable
    groupColumns=para.groupColumns
    orderGroups=para.orderGroups
    Opath=para.Opath
    prefix=para.prefix

    if not falist or not infoTable:
        sys.exit("{0}\n\tplease set -i and -I.".format(sys.argv[0]))

    os.system("mkdir -p "+Opath+'/eachGene/')

    ### calculate dinucleotide frequency
    for i in open(falist,'r'):
        i=i.strip()
        geneName=os.path.basename(i).split('.')[0]

        os.system('dinuc -enc '+ i+' '+Opath+'/eachGene/'+geneName+'.dinuc.txt -silent -nomenu -nowarn -noblk')

    ### gather all dinuc result
    info=pd.read_table(infoTable,sep="\t",index_col=0)

    dinuc_list=[]
    for i in glob(Opath+"/eachGene/*.dinuc.txt"):
        gene_name=os.path.basename(i)[:-10]

        dinuc=pd.read_table(i,sep="\s+",header=0)
        dinuc=dinuc.loc[dinuc.frame == 'all']
        del dinuc['frame']
        dinuc=dinuc.set_index('title')

        dinuc=dinuc.applymap(lambda x:x*16)
        dinuc=dinuc.reset_index()
        dinuc.insert(1,'gene',[gene_name]*dinuc.shape[0])
        # dinuc.insert(1,'group',dinuc['gene']+'.'+dinuc['title'])
        dinuc_list.append(dinuc)

    dinucs=pd.concat(dinuc_list)
    dinucs.insert(1,'ID',dinucs['gene']+'.'+dinucs['title'])
    dinucs.insert(2,'Groups',info.loc[dinucs.title,groupColumns].values)
    del dinucs['title']
    dinucs.set_index("ID",inplace=True)

    dinucs.to_csv(Opath+'/'+prefix+".allGene.dinuc_frequency.txt",sep="\t")

    #### save mean dinucs of each gene ####
    dinucs_mean=dinucs.groupby('gene').mean()
    dinucs_mean=dinucs_mean.round(2)

    writer=ExcelWriter(Opath+'/'+prefix+".allGene.mean_dinuc_frequency.xlsx",engine='xlsxwriter')
    dinucs_mean.to_excel(writer,sheet_name='sheet1')
    worksheet=writer.sheets["sheet1"]
    worksheet.conditional_format("B2:BH"+str(dinucs_mean.shape[0]+1),{'type':'3_color_scale'})
    writer.close()
    #### save mean dinucs of each gene ####

    if CpG_bar:
        os.system("Rscript "+os.path.dirname(os.path.realpath(__file__))+"/CpG.barPlot.R "+Opath+"/"+prefix+".allGene.dinuc_frequency.txt "+orderGroups+' '+Opath+'/'+prefix)

    ## Correlation analysis between GC1s-3s and CpG frequencies
    if CpG_GC_corr:
        ATGC_overAll=pd.read_table(ATGC,sep='\t')

        ATGC_gc=pd.melt(ATGC_overAll,id_vars=['Base','ID'])
        ATGC_gc=ATGC_gc.loc[ATGC_gc.Base.isin(['GC1','GC2','GC3'])]
        ATGC_gc['NEW_ID']=ATGC_gc.variable+'.'+ATGC_gc.ID
        del ATGC_gc['ID']
        ATGC_gc=ATGC_gc.rename(columns={'NEW_ID':'ID','variable':"Gene","value":"Percentage"})
        ATGC_gc=pd.pivot_table(ATGC_gc,index=['ID','Gene'],columns='Base',values='Percentage')

        ATGC_gc.reset_index(inplace=True)
        ATGC_gc=ATGC_gc.set_index("ID")
        ATGC_gc['CpG']=dinucs.loc[ATGC_gc.index,'CG']

        corr_results=pd.DataFrame()
        pvalue_results=pd.DataFrame()
        for i,j in ATGC_gc.groupby('Gene'):
            corr_gc1=spearmanr(j.GC1,j.CpG)[0]
            corr_gc2=spearmanr(j.GC2,j.CpG)[0]
            corr_gc3=spearmanr(j.GC3,j.CpG)[0]

            pvalue_gc1=spearmanr(j.GC1,j.CpG)[1]
            pvalue_gc2=spearmanr(j.GC2,j.CpG)[1]
            pvalue_gc3=spearmanr(j.GC3,j.CpG)[1]

            corr_results=corr_results.append([[i,corr_gc1,corr_gc2,corr_gc3]])
            pvalue_results=pvalue_results.append([[i,pvalue_gc1,pvalue_gc2,pvalue_gc3]])

        corr_results.columns=['Gene','GC1','GC2','GC3']
        pvalue_results.columns=['Gene','GC1','GC2','GC3']
        corr_results.set_index('Gene',inplace=True)
        pvalue_results.set_index('Gene',inplace=True)


        writer=ExcelWriter(Opath+'/'+prefix+".correlation_coefficient.CpG.GC.xlsx",engine='xlsxwriter')
        corr_results.to_excel(writer,sheet_name='Correlation coefficient')
        pvalue_results.to_excel(writer,sheet_name='P-value')

        worksheet=writer.sheets["Correlation coefficient"]
        worksheet.conditional_format("B2:D"+str(corr_results.shape[0]+1),{'type':'3_color_scale'})
        worksheet1=writer.sheets["P-value"]
        worksheet1.conditional_format("B2:D"+str(pvalue_results.shape[0]+1),{'type':'3_color_scale'})
        writer.close()

if __name__ == '__main__':
    main()
