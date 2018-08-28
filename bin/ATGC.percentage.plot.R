#!/usr/bin/env Rscript
argv <- commandArgs(T)
if (length(argv) != 4){
    stop("Rscript eachGene.CorrespondenceAnalysis.R [ATGC.content.overall] [ATGC.content.3rd.codon] [The order of the groupColumns][prefix of output]")
}

library(ggplot2)
library(plyr)
library(cowplot)

groupOrder <- strsplit(argv[3],split = ",")[[1]]

#### Nucleotide distribution frequency overall
dat <- read.table(argv[1],sep="\t",header=TRUE)
dat$labels <- round(dat$Percentage*100,1)
dat$Base <- ordered(dat$Base,c("G",'C','T','A'))
dat$Groups <- ordered(dat$Groups,groupOrder)
dat <- arrange(dat,Base)
dat<-ddply(dat,.(Groups,Gene),transform,position=cumsum(Percentage)-0.5*Percentage)
dat$Base <- ordered(dat$Base,c("A",'T','C','G'))

### sort by G
dat_G <- dat[dat$Base %in% c('G'),]
dat_G_sum <- ddply(dat_G,.(Gene,Base),summarise,sum=sum(Percentage))
Gene_sort_by_G <- as.character(arrange(dat_G_sum,desc(sum))$Gene)
dat$Gene <- ordered(dat$Gene,Gene_sort_by_G)

p1 <- ggplot(dat,aes(Groups,Percentage))+geom_bar(aes(fill=Base),stat="identity")+facet_grid(.~Gene)+scale_y_continuous(labels=scales::percent)

p1 <- p1+geom_text(aes(label=labels,y=position),size=2.5)

p1 <- p1+theme_bw()+theme(axis.title=element_text(size=14),axis.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),panel.grid.minor=element_blank(),panel.grid.major=element_blank())+labs(x="Groups",y="Percentage")+scale_fill_discrete(name="Base",labels=paste("%",c('A','T','C','G')))

### Nucleotide distribution frequency calculated only for the 3rd codon base

dat1 <- read.table(argv[2],sep="\t",header=TRUE)
dat1$labels <- round(dat1$Percentage*100,1)
dat1$Base <- ordered(dat1$Base,c("G3",'C3','T3','A3'))
dat1$Groups <- ordered(dat1$Groups,groupOrder)
dat1 <- arrange(dat1,Base)
dat1 <- ddply(dat1,.(Groups,Gene),transform,position=cumsum(Percentage)-0.5*Percentage)
dat1$Base <- ordered(dat1$Base,c("A3",'T3','C3','G3'))

### sort by G
dat1_G <- dat1[dat1$Base %in% c('G3'),]
dat1_G_sum <- ddply(dat1_G,.(Gene,Base),summarise,sum=sum(Percentage))
Gene_sort_by_G <-as.character(arrange(dat1_G_sum,desc(sum))$Gene)
dat1$Gene <- ordered(dat1$Gene,Gene_sort_by_G)

p2 <- ggplot(dat1,aes(Groups,Percentage))+geom_bar(aes(fill=Base),stat="identity")+facet_grid(.~Gene)+scale_y_continuous(labels=scales::percent)

p2 <- p2+geom_text(aes(label=labels,y=position),size=2.5)

p2 <- p2+theme_bw()+theme(axis.title=element_text(size=14),axis.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),panel.grid.minor=element_blank(),panel.grid.major=element_blank())+labs(x="Groups",y="Percentage")+scale_fill_discrete(name="Base",labels=paste("%",c('A3','T3','C3','G3')))

pdf(paste(argv[4],".ATGC.percentage.bar.pdf",sep=''),width=4*length(groupOrder),height=10)
plot_grid(p1,p2,ncol=1,align=c("hv"),labels=LETTERS[1:2],hjust=-1)
dev.off()
