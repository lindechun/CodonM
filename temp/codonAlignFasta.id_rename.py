#!/usr/bin/env python
# coding=utf-8

__author__ = "Lin Dechun"
__copyright__ = "Copyright 2018, BGI Research."
__credits__ = ["Lin Dechun"]
__version__ = "0.0.1"
__maintainer__ = "Lin Dechun"
__email__ = "lindechun@genomics.cn"

from glob import glob
from Bio import SeqIO
import pandas as pd
import os,sys

os.system('mkdir temp')
info=pd.read_table('/opt/projects/CodonM/data/SampleInfo.txt',sep='\t',index_col=0)

for j in glob('/opt/projects/CodonM/data/codonAlignFasta/*.fa'):
    temp_seq={}
    hd=open('temp/'+os.path.basename(j),'w')
    for i in SeqIO.parse(j,'fasta'):
        i.id=i.id.split('/')[1]    
        temp_seq[i.id]=i.seq
    for k in info.index:
        hd.write('>{0}\n{1}\n'.format(k,temp_seq[k]))
