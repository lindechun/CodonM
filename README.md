## CodonM

Comparative analysis pipeline of Population's Codon Usage

## Author

Dechun Lin

## Citation
Please cite the following article when using ggtree:

*Lin, D.*, Li, L., Xie, T., et al. (2018). Codon usage variation of Zika virus: The potential roles of NS2B and NS4A in its global pandemic. *Virus research* 247, 71-83.[Link](https://www.sciencedirect.com/science/article/pii/S016817021730597X?via%3Dihub)

## Introduction

A tool based on `python` to detect Population's Codon Usage Pattern and Bias within homologous gene of Virus complete genome

This lists the basic information for using [`CodonM`](https://github.com/lindechun/CodonM).

## Requirements

* A UNIX based operating system.

* python 2.7

* python packages: biopython, scipy, pandas, statsmodels

* [megacc](https://www.megasoftware.net/)

* [codonW](http://codonw.sourceforge.net)

## Installation

Download CodonM from GitHub. You'll need to add CodonM's bin directory to your $PATH.:

```
git clone https://github.com/lindechun/CodonM.git
/your/path/to/CodonM/CodonM.py -h
```

You will nedd to install all the dependencies packages and `codonw` by `pip` or `conda`.


## Basic test data set

See `/your/path/to/CodonM/data/` for test data set.

## Usage

```
$ CodonM.py -h
usage: CodonM.py [-h] [-v] [-i FALIST] [-I INFO] [-t SET_TABLE]
                 [-j ORDERGENES] [-g GROUPCOLUMNS] [-G ORDERGROUPS] [-CA]
                 [-HC] [-ATGC_bar_plot] [-ENC_bar_plot] [-ENC_GC3_plot]
                 [-CpG_bar] [-CpG_GC_corr] [-o OPATH] [-p PREFIX]

The flow of Comparative analysis pipeline of Codon Usage

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    show program's version number and exit
  -i FALIST        a list file consists of individual Genes, based codon align
                   (megacc)
  -I INFO          sample info table
  -t SET_TABLE     Genetic Code Table, default is 1
  -j ORDERGENES    The order of genes
  -g GROUPCOLUMNS  the columns of info table, which is the groups of sample.
  -G ORDERGROUPS   the order of -g, which is columns of info table. When
                   -ENC_bar_plot or -CpG, you should set -G
  -CA              Whether Correspondence Analysis for RSCU matrix of all
                   genes
  -HC              Whether Hierarchical Cluster for RSCU matrix of all genes
  -ATGC_bar_plot   Whether plot ATGC percentage of each gene and groupColumns
                   that you appoint
  -ENC_bar_plot    Whether plot ATGC percentage of each gene and groupColumns
                   that you appoint
  -ENC_GC3_plot    Whether plot ATGC percentage of each gene and groupColumns
                   that you appoint
  -CpG_bar         Whether plot ATGC percentage of each gene and groupColumns
                   that you appoint
  -CpG_GC_corr     Whether correlation analysis between GC1s-3s and CpG
                   frequencies
  -o OPATH         Path of Output File
  -p PREFIX        Output file prefix
```

See `/your/path/to/CodonM/test/demo_*.sh` for example of Codon usage analysis about 100 Zika sequences.

```
$ cat demo_python_flow.sh
python ../bin/CodonM.py -i data.list -I ../data/SampleInfo.txt -j "C,prM,E,NS1,NS2A,NS2B,NS3,NS4A,NS4B,NS5" -g SubLineage -G "African,South_Eastern_Asian,Oceania-American" -CA -HC -ATGC_bar_plot -ENC_GC3_plot -ENC_bar_plot -CpG_GC_corr -CpG_bar -o results -p Zika
```