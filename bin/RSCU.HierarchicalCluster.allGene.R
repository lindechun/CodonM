#!/usr/bin/env Rscript
argv <- commandArgs(T)
if (length(argv) != 5){
    stop("Rscript allGene.RSCU.HierarchicalCluster.R [SampleInfo.expand.allGene.txt] [allGene.RSCU.txt] [The order of each Gene] [groupColumns] [prefix of output]")
}

library(pheatmap)
info <- read.table(argv[1],sep="\t", header=T, row.names=1)
table <- read.table(argv[2],row.names=1)

geneOrder <- strsplit(argv[3],split = ",")[[1]]
info$gene <- ordered(info$gene,geneOrder)

pdf(paste(argv[5],".allGene.RSCU.HierarchicalCluster.pdf",sep=""),width=10,height=8)
pheatmap(table,show_rownames=F,show_colnames=T,annotation_row=info[c("gene",argv[4])])
dev.off()
