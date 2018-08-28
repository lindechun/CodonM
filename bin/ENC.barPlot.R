#!/usr/bin/env Rscript
argv <- commandArgs(T)
if (length(argv) != 3){
    stop("Rscript ENC.barPlot.R [allGene.ENC_GC3s.txt] [The order of groupColumns][prefix of output]")
}

library(ggplot2)
library(dplyr)

dat=read.table(argv[1],row.names=1,sep="\t",header=T)

groupOrder <- strsplit(argv[2],split = ",")[[1]]

dat$Groups <- ordered(dat$Groups,groupOrder)

temp_mean<-dat %>% group_by(gene,Groups) %>% dplyr::summarise(mean=mean(ENC))
temp_std<-dat %>% group_by(Groups,gene) %>% dplyr::summarise(std=sd(ENC))
temp_mean['std'] <- temp_std['std']

###### order gene ####
gene_enc_mean<-dat %>% group_by(gene) %>% dplyr::summarise(mean=mean(ENC))
order_gene<-gene_enc_mean[order(gene_enc_mean$mean,decreasing = T),]$gene

temp <- dat %>% group_by(gene,Groups) %>% dplyr::summarise(mean=mean(ENC))
temp1 <- temp %>% group_by(gene) %>% mutate(volume = mean-lag(mean, default = mean[1]))
temp2 <- temp1 %>% group_by(gene) %>% dplyr::summarise(mean=mean(volume)) %>% arrange(-mean)

up_gene=c()
lower_gene=c()
for (i in order_gene){
  if (temp2[temp2$gene == i,]$mean>0){
    up_gene=c(up_gene,i)
  }else{
    lower_gene=c(lower_gene,i)
  }
}

order_gene <- c(up_gene,lower_gene)
###### order gene ####

temp_mean$gene<-ordered(temp_mean$gene,order_gene)

pdf(file=paste(argv[3],".ENC.barPlot.pdf",sep=''),width=6,height = 5)
p1 <- ggplot(temp_mean,aes(gene,mean,group=Groups))+geom_bar(aes(fill=Groups),stat="identity",position="dodge",width=0.75)+geom_errorbar(aes(ymax = mean+std, ymin=mean - std),position=position_dodge(0.75),width=0.2,size=0.25)+theme_bw()+theme(axis.title=element_text(size=14),axis.text=element_text(size=12),panel.grid.minor=element_blank(),panel.grid.major=element_blank())+labs(y="Average ENC",x="Gene")+scale_y_continuous(expand=c(0.03,0))+coord_cartesian(ylim=c(0,62),expand = T)

p1 + theme(legend.position = "bottom")+theme(legend.text=element_text(size=12,color='black',family='Times',hjust=0),legend.title=element_text(size = 13,color = 'black',family = 'Times', hjust = 0),legend.background = element_blank())
dev.off()
