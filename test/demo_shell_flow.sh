#!/bin/bash
dirpath=$(cd `dirname $0`; pwd)
prefix='Zika'
dat='data.list'

# Step01. BasicCodonUsage
mkdir -p results/01.BasicCodonUsage

while read lines
do
    DataId=`basename $lines .fa`
    python2 $dirpath/../bin/BasicCodonUsage.py -i $lines -t 1 -o results/01.BasicCodonUsage/$DataId
done < $dat

# Step02. RSCU pattern
mkdir -p results/02.RSCU_pattern
python2 $dirpath/../bin/RSCU.pattern.py -i results/01.BasicCodonUsage -I ../data/SampleInfo.txt -CA -HC -j "C,prM,E,NS1,NS2A,NS2B,NS3,NS4A,NS4B,NS5" -g SubLineage -o results/02.RSCU_pattern -p $prefix

# Step03. ATGC content analysis
mkdir -p results/03.ATGC_content
python2 $dirpath/../bin/ATGC.content.py -i $dat -I $dirpath/../data/SampleInfo.txt -g SubLineage -o results/03.ATGC_content -p $prefix -ATGC_bar_plot -G "African,South_Eastern_Asian,Oceania-American"

# Step04. ENC and ENC-GC3s analysis
mkdir -p results/04.ENC_GC3s
python2 $dirpath/../bin/ENC_GC3s.py -i $dat -I $dirpath/../data/SampleInfo.txt -o results/04.ENC_GC3s -g SubLineage -p $prefix -ENC_bar_plot -ENC_GC3_plot -G "African,South_Eastern_Asian,Oceania-American" -j "C,prM,E,NS1,NS2A,NS2B,NS3,NS4A,NS4B,NS5" -ATGC results/03.ATGC_content/${prefix}.ATGC_Content.overAll.txt

## Step05. Dinucleotide. particularly,CpG
mkdir -p results/05.Dinucleotide
python2 $dirpath/../bin/Dinucleotide.py -i $dat -I $dirpath/../data/SampleInfo.txt -o results/05.Dinucleotide -p $prefix -g SubLineage -CpG_bar -G "African,South_Eastern_Asian,Oceania-American" -CpG_GC_corr -ATGC results/03.ATGC.content/${prefix}.ATGC_Content.overAll.txt
