#!/usr/bin/env Rscript
argv <- commandArgs(T)
if (length(argv) != 5){
    stop("Rscript eachGene.CorrespondenceAnalysis.R [SampleInfo.txt] [The path of results of Step 01(intermediate Result)] [groupColumns] [The order of each Gene] [prefix of output]")
}

library(ggplot2)
library(Cairo)
library(ade4)
library(cowplot)

info <- read.table(argv[1],sep="\t", header=T, row.names=1)
geneOrder <- strsplit(argv[4],split = ",")[[1]]

plots <- list()
for (i in geneOrder){
    table=read.table(paste(argv[2],'/',i,'/',i,'.RSCU.txt',sep=""),header=T, row.names=1)
    table_info <- info[colnames(table),]
    rscu.pca <- dudi.coa(t(table), scannf=F, nf=ncol(table))
    attach(rscu.pca)
    f1 <- round(100*eig[1]/sum(eig),1)
    f2 <- round(100*eig[2]/sum(eig),1)
    xlab <- paste("f1'(",f1,"%)",sep="")
    ylab <- paste("f2'(",f2,'%)',sep="")
    
    dat_li <- cbind(li[,1:2],table_info)
    
    p <- ggplot(dat_li,aes(Axis1,Axis2))+geom_point(aes_string(color=argv[3]))
    p <- p + geom_hline(yintercept=0, size = 0.5, linetype=2)+ geom_vline(xintercept=0, size = 0.5, linetype = 2) + labs(x=xlab, y=ylab)
    
    legend_b <- get_legend(p + theme(legend.position="bottom"))
    
    plots[[i]] <- p+theme_bw()+theme(legend.position="none")+theme(axis.title=element_text(size=12),axis.text=element_text(size=10),plot.title=element_text(hjust=-.3),panel.grid.minor=element_blank(),panel.grid.major=element_blank())+labs(title=i)+theme(plot.title=element_text(hjust=0.5))
}

dev.off()

plot_row <- ceiling(length(geneOrder)/2)
plot_height <- 4*plot_row + 1

pdf(file=paste(argv[5],".eachGene.RSCU.coa.pdf",sep=''), width=8, height=plot_height)

prow <- plot_grid(plotlist=plots,ncol=2,align=c("hv"),labels=LETTERS[1:length(geneOrder)],hjust=-1)
plot_grid(prow,legend_b, ncol = 1, rel_heights = c(plot_height-1,1))

dev.off()
