#!/usr/bin/env Rscript
argv <- commandArgs(T)
if (length(argv) != 4){
    stop("Rscript ENC_GC3s.individualGenes.R [allGene.ENC_GC3s.txt] [The order of groupColumns] [the order of genes] [prefix of output]")
}

library(ggplot2)
library(cowplot)

rect_point<- function(temp){
  x_min <- min(temp$GC3)
  x_max <- max(temp$GC3)
  y_min <- min(temp$ENC)
  y_max <- max(temp$ENC)
  x_radius=(x_max-x_min)/2
  y_radius=(y_max-y_min)/2
  
  if(x_radius*70 >y_radius){
    y_radius = x_radius*70
  }else{
    x_radius = y_radius/70
  }
  
  x_center=(x_min+x_max)/2
  y_center=(y_min+y_max)/2
  
  x1=round(x_center-x_radius-0.01,2)
  x2=round(x_center+x_radius+0.01,2)
  y1=y_center-y_radius-1.0
  y2=y_center+y_radius+1.0
  return(data.frame(x1=c(x1,x2,x1,x1),x2=c(x2,x2,x2,x1),y1=c(y1,y1,y2,y2),y2=c(y1,y2,y2,y1)))
}

table=read.table(argv[1],row.names=1,sep="\t",header=T)
groupOrder <- strsplit(argv[2],split = ",")[[1]]

table$Groups <- ordered(table$Groups,groupOrder)

x=seq(0,1,0.01)
y=2+x+29/(x^2+(1-x)^2)
enc_gc3s<- data.frame(GC3s=x,ENC=y)

#### each gene
geneOrder <- strsplit(argv[3],split = ",")[[1]]

plots <- list()

for (i in geneOrder){
  temp<-table[table$gene == i,]
  rect_values <-rect_point(temp)

  ## original
  p1 <- ggplot(temp,aes(GC3s,ENC))+geom_point(aes(color=Groups,shape=Groups))+geom_segment(data=rect_values,aes(x=x1,xend=x2,y=y1,yend=y2),color="grey20",linetype="dashed",size=0.3)+
    annotate("segment",x=mean(range(rect_values$x1)),xend=0.5,y=min(rect_values$y1)-2,yend=37,colour="black",size=0.25,arrow=arrow(angle=20,length=unit(0.1,'inches')))+
    scale_x_continuous(labels=c(0,0.2,0.4,0.6,0.8,1),breaks=c(0,0.2,0.4,0.6,0.8,1),limits=c(0,1))+scale_y_continuous(labels=c(0,10,20,30,40,50,60,70),breaks=c(0,10,20,30,40,50,60,70),limits=c(0,70))+coord_cartesian(expand = F)+
    geom_line(data=enc_gc3s,color="#7f7f7f")
  p1 <- p1 + theme_bw() + theme(legend.position="none")+theme(axis.title=element_text(size=10),axis.text=element_text(size=10),plot.title=element_text(hjust=-.3,size=10),panel.grid.minor=element_blank(),panel.grid.major=element_blank())+labs(title=i)+theme(plot.title=element_text(hjust=0.5))

  ## enlarge
  p2 <- ggplot(temp,aes(GC3s,ENC))+geom_point(aes(color=Groups,shape=Groups))+geom_line(data=enc_gc3s,color="#7f7f7f")+
    scale_x_continuous(labels=seq(min(rect_values$x1),max(rect_values$x1),(max(rect_values$x1)-min(rect_values$x1))/2),breaks=seq(min(rect_values$x1),max(rect_values$x1),(max(rect_values$x1)-min(rect_values$x1))/2),limits=range(rect_values$x1))+scale_y_continuous(labels=seq(min(rect_values$y1),max(rect_values$y1),(max(rect_values$y1)-min(rect_values$y1))/2),breaks=seq(min(rect_values$y1),max(rect_values$y1),(max(rect_values$y1)-min(rect_values$y1))/2),limits=range(rect_values$y1))+coord_cartesian(expand = T)
  p2 <- p2+theme_bw()+theme(legend.position="none")+theme(axis.title=element_text(size=8),axis.text=element_text(size=8),plot.title=element_text(hjust=-.3,size=8),panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(size=0.5))+labs(title=i)+theme(plot.title=element_text(hjust=0.5),panel.border=element_blank(),axis.line=element_line(colour='black',size=0.25))+theme(plot.background = element_rect(colour="black",linetype=2,size=0.5))+labs(x=NULL,y=NULL,title=NULL)
  
  legend_b <- get_legend(p1+theme(legend.position="bottom"))
  
  plots[[i]] <- p1 + annotation_custom(ggplotGrob(p2), xmin=0.2, xmax=0.8, ymin=2, ymax=35)
}
dev.off()

plot_row <- ceiling(length(geneOrder)/2)
plot_height <- 3.5*plot_row + 1

pdf(file=paste(argv[4],".ENC-GC3s.individualGenes.pdf",sep=''),width=8,height=plot_height)
prow <- plot_grid(plotlist=plots,ncol=2,align=c("hv"),labels=LETTERS[1:10],hjust=-1)
plot_grid( prow,legend_b, ncol = 1, rel_heights = c(20,1))
dev.off()
