#!/usr/bin/env Rscript
argv <- commandArgs(T)
if (length(argv) != 4){
    stop("Rscript allGene.CorrespondenceAnalysis.R [SampleInfo.expand.allGene.txt] [allGene.RSCU.txt] [groupColumns] [prefix of output]")
}

library(ggplot2)
library(Cairo)
library(ade4)

## read file
info <- read.table(argv[1],row.names=1,sep="\t",header=T)
table <- read.table(argv[2],row.names=1)

rscu.ca<-dudi.coa(table,scannf=F,nf=nrow(table))
attach(rscu.ca)
dat_li<-cbind(li,info)

f1=round(100*eig[1]/sum(eig),1)
f2=round(100*eig[2]/sum(eig),1)
xlab <- paste("f1'(",f1,"%)",sep="")
ylab <- paste("f2'(",f2,'%)',sep="")

shapeLibrary <- c(9,1,2,5,0,4,6,7,8,3,11,10,12,13)

pdf(file=paste(argv[4],".allGene.RSCU.coa.pdf",sep=''),width=9,height=6)

ggplot(dat_li,aes(Axis1,Axis2))+
    geom_point(aes_string(color=argv[3],shape="gene"))+
    geom_hline(yintercept=0,size = 0.5,linetype=2)+
    geom_vline(xintercept=0,size = 0.5,linetype = 2)+
    labs(x=xlab,y=ylab)+
    scale_shape_manual("Gene",values=shapeLibrary[seq(1,length(unique(dat_li$gene)))])+
    theme_bw()+theme(axis.title=element_text(size=15),axis.text=element_text(size=14),panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
    theme(legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text=element_text(size=12))
dev.off()

