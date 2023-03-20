# load data
sequenza = readRDS("Input/sequenza.rds")
load("Input/SNVS.Rdata")
load("Input/NCG_CRC_cancer_genes.Rdata")

library(plyr)

get.tc.correction.somatic = function( obs, tc, CNt, CNn=2){
  return( min(  obs * ( 1 + (  ( CNn*(1-tc) )/( CNt * tc) ) )  ,  1))
}



for(i in 1:nrow(snvs)){
  snvs$frq.tc[i] = get.tc.correction.somatic(snvs$obs.VAF[i], snvs$Aberrant_Cell_Fraction,
                                                                   snvs$major_cn[i]+snvs$minor_cn[i], snvs$normal_cn[i] )
  snvs$clonality[i] = min(c(snvs$frq.tc[i]*2, 1))
}

save(snvs, file="Results/SNVS.clonality.Rdata")
write.csv(snvs, file="Results/Mutations.csv", row.names=F)

pdf(file="Results/Frequency.pdf", h=5, w=12)
par(mfrow=c(1,2))
hist(snvs$obs.VAF, main='Frequency', xlab ='Frequency', ylab="Counts", xlim=c(0,1), ylim=c(0,700))
hist(snvs$frq.tc, main='Frequency\npost TC correction', xlab ='Frequency', ylab="Counts", xlim=c(0,1), ylim=c(0,700))
dev.off()

pdf(file="Results/Clonality.pdf", h=5, w=12)
par(mfrow=c(1,2))
hist(snvs$frq.tc, main='Frequency\npost TC correction', xlab ='Frequency', ylab="Counts", xlim=c(0,1), ylim=c(0,700))
hist(snvs$clonality, main='Clonality', xlab ='Frequency', ylab="Counts", xlim=c(0,1), ylim=c(0,700))
dev.off()

require(ggplot2)
require(ggrepel)

base_breaks_x = function(x, br, la, reverse=T){
  d <- data.frame(y=-Inf, yend=-Inf, x=0, xend=100)
  if(reverse){
    return(list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_x_reverse(breaks=br, labels=la)))
  }else{
    return(list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_x_continuous(breaks=br, labels=la)))
  }
}

base_breaks_y = function(m, br, la ){
  d <- data.frame(x=Inf, xend=Inf, y=0, yend=m)
  list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_y_continuous(breaks=br, labels=la))
}


p =ggplot(snvs, aes(x=clonality)) + geom_density(fill=rgb(172,152,199, maxColorValue = 255), alpha=.5,aes(y = ..count..)) +
  geom_vline(xintercept=0.80, linetype="dashed", color="grey")+
  geom_vline(xintercept=0.35, linetype="dashed", color="grey")

gg  = ggplot_build(p)$data[[1]]
# gga = unique(ggplot_build(p)$plot[[1]]$Alteration_type) #     ptmp = ggplot(dd, aes(x=cell,fill=Alteration_type)) + geom_histogram() #     gtmp =  ggplot_build(ptmp)$data[[1]]
ymax = max(gg$y)
ymed = ymax/2
lmax = ifelse( ymax>=0.5, round(ymax), 1)
lmed = ifelse(lmax==1, 0.5, ifelse( ymed>0.5, round(ymed), 0.5))

p = p +theme_bw()+scale_x_reverse()+ylab('Expected Alterations')+xlab('Alteration clonality (%)')


drvrs=ddply(
  subset(snvs, Hugo_Symbol%in%crc_cancer_genes$symbol),
  .(clonality), summarise, gname=unique(Hugo_Symbol)
)

drvrs$x = NA
drvrs$y = NA
for (i in 1:nrow(drvrs)){
  x = drvrs$clonality[i]
  drvrs$x[i] = gg$x[which(abs(gg$x-x)==min(abs(gg$x-x)))]
  drvrs$y[i] = gg$y[which(abs(gg$x-x)==min(abs(gg$x-x)))]
}

pdf(file="Results/DensityPlot.pdf",h=6,w=9)
  p +
  geom_point(data=drvrs, aes(x = x, y=y), colour="black", show.legend = F)+
  geom_text_repel(data=drvrs, aes(x=x, y=y, label = gname))
dev.off()


