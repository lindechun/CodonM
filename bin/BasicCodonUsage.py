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
import pandas as pd
import numpy as np
from pandas import ExcelWriter
import xlsxwriter
from Bio import SeqIO,AlignIO
from Bio.Seq import Seq,translate,reverse_complement
from Bio.Data import CodonTable
from subprocess import Popen,PIPE
from copy import copy
from multiprocessing import Pool,Queue,Lock,Manager  #使用多进程

def tranTable():
    print '''
    Genetic Code Table

        ID        List of Names
        1    ->    ['Standard', 'SGC0']
        2    ->    ['Vertebrate Mitochondrial', 'SGC1']
        3    ->    ['Yeast Mitochondrial', 'SGC2']
        4    ->    ['Mold Mitochondrial', 'Protozoan Mitochondrial', 'Coelenterate Mitochondrial', 'Mycoplasma', 'Spiroplasma', 'SGC3']
        5    ->    ['Invertebrate Mitochondrial', 'SGC4']
        6    ->    ['Ciliate Nuclear', 'Dasycladacean Nuclear', 'Hexamita Nuclear', 'SGC5']
        9    ->    ['Echinoderm Mitochondrial', 'Flatworm Mitochondrial', 'SGC8']
        10    ->    ['Euplotid Nuclear', 'SGC9']
        11    ->    ['Bacterial', 'Plant Plastid']
        12    ->    ['Alternative Yeast Nuclear']
        13    ->    ['Ascidian Mitochondrial']
        14    ->    ['Alternative Flatworm Mitochondrial']
        15    ->    ['Blepharisma Macronuclear']
        16    ->    ['Chlorophycean Mitochondrial']
        21    ->    ['Trematode Mitochondrial']
        22    ->    ['Scenedesmus obliquus Mitochondrial']
        23    ->    ['Thraustochytrium Mitochondrial']
        24    ->    ['Pterobranchia Mitochondrial']
    '''

    """ ###氨基酸三字母和一字母对应表
    amacid_dict={'F':'Phe','S':'Ser','Y':'Tyr','C':'Cys','L':'Leu','END':'Ter','W':'Trp',
                 'P':'Pro','H':'His','R':'Arg','Q':'Gln','I':'Ile','T':'Thr','N':'Asn',
                 'K':'Lys','M':'Met','V':'Val','A':'Ala','D':'Asp','G':'Gly','E':'Glu'}
    """

def parseCommand():
    parser = ArgumentParser(description='This script use to convert nucleotide to amino acid',version='1.0.0')
    parser.add_argument('-i', action='store', dest='inFa', help='Input is a Align\'s fasta file, based codon way')
    parser.add_argument('-t', action='store', dest='set_table', default=1, type=int, help='Genetic Code Table, default is 1')
    parser.add_argument('-o', action='store', dest='Opath', help='Path of Output File')
    return parser.parse_args()

def read_Queue(q,Opath):
    "处理子进程数据"
    df=pd.DataFrame()
    while not q.empty():
        if df.empty:
            df=q.get()
        else:
            df=df.join(q.get())

    df=df.round(2)

    ### save Codon usage
    writer=ExcelWriter(Opath+'.RSCU.xlsx',engine='xlsxwriter')
    df.T.to_excel(writer,sheet_name='sheet1')
    worksheet=writer.sheets["sheet1"]
    worksheet.conditional_format("B4:BH"+str(df.T.shape[0]+3),{'type':'3_color_scale'})
    writer.close()
    ### save Codon usage

    df.reset_index(inplace=True)
    df['Codon']=df['AmAcid']+'.'+df['Codon']
    del df['AmAcid']
    df.set_index("Codon",inplace=True)
    df.to_csv(Opath+".RSCU.txt",sep='\t')

def condonCount(records,codon,Opath):
    del codon['Number']
    for fasta in records:
        num=0
        fa_id=fasta.id
        codon[fa_id]=''
        while num+3 <= len(fasta.seq):
            temp=fasta.seq[num:num+3]
            if codon.loc[codon[codon.Codon==temp].index,fa_id].values:
                codon.loc[codon[codon.Codon==temp].index,fa_id]+=1
            else:
                codon.loc[codon[codon.Codon==temp].index,fa_id]=1
            num+=3

    codon['ID']=codon['AmAcid']+'.'+codon['Codon']
    del codon['AmAcid']
    del codon['Codon']
    codon.set_index('ID',inplace=True)
    codon=codon.replace('',0)
    codon.index.name='Codon'
    codon.to_csv(Opath+'.codonStats.eachSample.txt',sep='\t')

def codonAmAcidStat(fasta,codon,Opath,q,l):
    '''统计序列的密码子使用情况'''
    num=0
    while num+3 <= len(fasta.seq):
        temp=fasta.seq[num:num+3]
        if codon.loc[codon[codon.Codon==temp].index,'Number'].values:
            codon.loc[codon[codon.Codon==temp].index,'Number']+=1
        else:
            codon.loc[codon[codon.Codon==temp].index,'Number']=1
        num+=3

    codon.Number=codon.Number.replace('',0)
    amacid_filter=['A','C','D','E','F','G','H','I','K','L','N','P','Q','R','S','T','V','Y']

    codon=codon.set_index("AmAcid").sort_index()
    codon=codon.loc[amacid_filter]

    codon.reset_index(inplace=True)
    AmAcid_sum=codon.groupby("AmAcid")[['Number']].sum()
    AmAcid_count=codon.groupby("AmAcid")[['Number']].count()

    codon=codon.set_index("AmAcid").sort_index()
    codon['/1000']=codon.Number*1000/len(fasta.seq)*3
    codon['Fraction']=codon.Number/AmAcid_sum.Number
    codon['RSCU']=codon.Number/AmAcid_sum.Number*AmAcid_count.Number
    codon=codon.fillna(0)

    l.acquire()  ##进程锁
    handle=open(Opath+".codonUsage.eachSample.txt",'a+')
    handle.write('>'+fasta.id+"\n")
    handle.write(codon.index.name+'\t'+'\t'.join(list(codon.columns))+"\n")
    for k,v in codon.iterrows():
        handle.write('{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4:.2f}\n'.format(k,v[0],v[1],v[2],v[3]))
    handle.write("\n")
    handle.close()
    l.release()

    ##单独统计RSCU
    codon.reset_index(inplace=True)
    temp=codon.groupby(['AmAcid','Codon']).sum()
    temp=temp.loc[:,['RSCU']]
    temp.rename(columns={'RSCU':fasta.id},inplace=True)

    q.put(temp)

def main():
    para=parseCommand()
    Opath=para.Opath
    inFa=para.inFa
    set_table=para.set_table

    if not para.inFa:
        print "Please set -i\nsee -h for help"
        tranTable()
        sys.exit(0)

    if not os.path.exists(Opath):
        os.system('mkdir -p '+Opath)

    Opath=Opath+'/'+Opath.split('/')[-1]

    standard_table = CodonTable.unambiguous_dna_by_id[para.set_table]  #密码表
    codon=pd.DataFrame(data=standard_table.forward_table,index=['AmAcid']).T  #建密码子数据框

    codon.index.name='Codon'
    codon.reset_index(inplace=True)
    codon['Number']=0

    writerQueue = Pool(processes = 4)     # 4 process
    m=Manager()
    l=m.Lock()
    q=m.Queue()
    thread=[]

    seq_record=[]
    for fasta in SeqIO.parse(inFa, 'fasta'):
        # fasta.seq=Seq(str(fasta.seq).replace('-',''))

        thread.append(writerQueue.apply_async(codonAmAcidStat,(copy(fasta),codon,Opath,q,l)))
        seq_record.append(fasta)

    writerQueue.close()
    writerQueue.join()

    condonCount(seq_record,copy(codon),Opath)

    ###处理子进程数据队列
    readQueue = Pool()
    readQueue.apply_async(read_Queue,(q,Opath))
    readQueue.close()
    readQueue.join()

if __name__  == "__main__":
    main()
