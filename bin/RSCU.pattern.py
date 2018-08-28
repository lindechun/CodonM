#!/usr/bin/env python
# coding=utf-8

__author__ = "Lin Dechun"
__copyright__ = "Copyright 2016, BGI Research."
__credits__ = ["Lin Dechun"]
__version__ = "0.0.1"
__maintainer__ = "Lin Dechun"
__email__ = "lindechun@genomics.cn"

from argparse import ArgumentParser
import pandas as pd
from glob import glob
import os,sys

def parseCommand():
    parser = ArgumentParser(description='Compare RSCU pattern of all genes',version='1.0.0')
    parser.add_argument('-i', action='store', dest='firstStepPath', help='The path of results of Step01 (01.BasicCodonUsage)')
    parser.add_argument('-CA',action='store_true', dest='CA', default=False,help=' Whether Correspondence Analysis for RSCU matrix of all genes')
    parser.add_argument('-HC', action='store_true', dest='HC', default=False,help='Whether Hierarchical Cluster for RSCU matrix of all genes')
    parser.add_argument('-I', action='store', dest='infoTable', help=' sample info table. When you set -CA or -HC, you should set -I')
    parser.add_argument('-j', action='store', dest='orderGenes',help='the order of -g, which is columns of info table. When -CA or -HC, you should set -j')
    parser.add_argument('-g', action='store', dest='groupColumns',help='the columns of info table, which is the groups of sample. When -CA or -HC, you should set -g')
    parser.add_argument('-o', action='store', dest='Opath',default='./', help='Path of Output File')
    parser.add_argument('-p', action='store', dest='prefix', help='prefix of Output File')
    return parser.parse_args()

def main():
    para=parseCommand()
    firstStepPath=para.firstStepPath
    infoTable=para.infoTable
    CA=para.CA
    HC=para.HC
    orderGenes=para.orderGenes
    groupColumns=para.groupColumns
    Opath=para.Opath
    prefix=para.prefix

    ## Step01. Gather RSCU for each genes
    RSCUs_list=[]

    for RSCUfile in glob(firstStepPath+'/*/*.RSCU.txt'):
        geneName=os.path.basename(RSCUfile)[:-9]
        RSCU=pd.read_table(RSCUfile,sep="\t",index_col=0,header=None)
        RSCU=RSCU.T
        RSCU.insert(0,'gene',[geneName]*RSCU.shape[0])
        RSCUs_list.append(RSCU)

    RSCUs=pd.concat(RSCUs_list)
    RSCUs.insert(1,'Group',RSCUs['gene']+'.'+RSCUs['Codon'])

    info=pd.read_table(infoTable,index_col=0)
    info=info.loc[RSCUs.Codon]
    info.reset_index(inplace=True)
    info.index=list(RSCUs['Group'])
    info.insert(0,'gene',list(RSCUs['gene']))

    del RSCUs['gene'],RSCUs['Codon']
    RSCUs.set_index("Group",inplace=True)
    RSCUs.index.name=''

    RSCUs.to_csv(Opath+'/'+prefix+".allGene.RSCU.txt",sep="\t")
    info.to_csv(Opath+'/'+prefix+".SampleInfo.expand.allGene.txt",sep="\t")

    ## Step02. Correspondence Analysis for RSCU pattern--allGene
    os.system("Rscript "+os.path.dirname(os.path.realpath(__file__))+"/RSCU.CorrespondenceAnalysis.allGene.R "+Opath+'/'+prefix+".SampleInfo.expand.allGene.txt "+Opath+'/'+prefix+".allGene.RSCU.txt "+groupColumns+' '+Opath+'/'+prefix)
    # Rscript $dirpath/../bin/RSCU.CorrespondenceAnalysis.allGene.R results/03.RSCU.pattern/$prefix.SampleInfo.expand.allGene.txt results/03.RSCU.pattern/$prefix.allGene.RSCU.txt SubLineage results/03.RSCU.pattern/$prefix

    ## Step03. Correspondence Analysis for RSCU pattern--individualGenes
    os.system("Rscript "+os.path.dirname(os.path.realpath(__file__))+"/RSCU.CorrespondenceAnalysis.individualGenes.R "+infoTable+' '+firstStepPath+' '+groupColumns+' '+orderGenes+' '+Opath+'/'+prefix)

    # Rscript $dirpath/../bin/RSCU.CorrespondenceAnalysis.individualGenes.R $dirpath/../data/SampleInfo.txt results/01.intermediateResult/ SubLineage "C,prM,E,NS1,NS2A,NS2B,NS3,NS4A,NS4B,NS5" results/03.RSCU.pattern/$prefix

    ## Step04. HierarchicalCluster for RSCU--allGene
    os.system("Rscript "+os.path.dirname(os.path.realpath(__file__))+"/RSCU.HierarchicalCluster.allGene.R "+Opath+'/'+prefix+".SampleInfo.expand.allGene.txt "+Opath+'/'+prefix+".allGene.RSCU.txt "+orderGenes+' '+groupColumns+' '+Opath+'/'+prefix)

    # Rscript $dirpath/../bin/RSCU.HierarchicalCluster.allGene.R results/03.RSCU.pattern/$prefix.SampleInfo.expand.allGene.txt results/03.RSCU.pattern/$prefix.allGene.RSCU.txt "C,prM,E,NS1,NS2A,NS2B,NS3,NS4A,NS4B,NS5" SubLineage results/03.RSCU.pattern/$prefix

if __name__ == '__main__':
    main()