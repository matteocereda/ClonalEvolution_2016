load("Input/ASCAT.Rdata")
load("Input/NCG_CRC_cancer_genes.Rdata")

write.csv(ascat, file= "Results/ASCAT.results.csv",row.names=F)

# Produce ASCAT plot
png(file="Results/ASCAT.png",w=2000, h=1000, res=200)
y_limit=5
len=sum(ascat$nProbes)
twoColours=T
par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
ticks=seq(0, y_limit, 1)

A_rle = list(lengths=ascat$nProbes, values=ascat$nA)
B_rle=list(lengths=ascat$nProbes, values=ascat$nB)
plot(c(1,len), c(0,y_limit), type = "n", xaxt = "n", yaxt="n", main = NULL, xlab = "", ylab = "")
axis(side = 2, at = ticks)
abline(h=ticks, col="lightgrey", lty=1)

colourMinor="yellow"
colourTotal="blue"
#A_rle<-rle(nAfullPlot)
start=0
#plot minor allele copy number
for(i in 1:length(B_rle$values)){
  val<-B_rle$values[i]
  size<-B_rle$lengths[i]
  rect(start, (val-0.07), (start+size-1), (val+0.07), col=ifelse((twoColours & val>=y_limit), adjustcolor(colourMinor, red.f=0.75,green.f=0.75,blue.f=0.75), colourMinor), border=ifelse((twoColours & val>=y_limit), adjustcolor(colourMinor, red.f=0.75,green.f=0.75,blue.f=0.75), colourMinor))
  start=start+size
}
start=0
#plot total copy number
for(i in 1:length(A_rle$values)){
  val<-A_rle$values[i]
  size<-A_rle$lengths[i]
  rect(start, (val-0.07), (start+size-1), (val+0.07), 
       col=ifelse((twoColours & val>=y_limit), 
                  adjustcolor(colourTotal,red.f=0.75,green.f=0.75,blue.f=0.75), colourTotal), 
       border=ifelse((twoColours & val>=y_limit), adjustcolor(colourTotal,red.f=0.75,green.f=0.75,blue.f=0.75), 
                     colourTotal))
  start=start+size
}

#B_rle<-rle(nBfullPlot)
chr.segs = ddply( ascat, .(Chr), summarise, n=sum(nProbes))
chr.segs$cs = cumsum(chr.segs$n)

chrk_tot_len = 0
abline(v=0,lty=1,col="lightgrey")
for (i in 1:nrow(chr.segs)) {
  abline(v=chr.segs$cs[i],lty=1,col="lightgrey")
  
  chrk_tot_len_prev = chrk_tot_len
  chrk_tot_len = chr.segs$cs[i]
  vpos = chrk_tot_len;
  tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
  text(tpos,y_limit,ifelse(i<23,sprintf("%d",i),ifelse(i==23,"X","Y")), pos = 1, cex = 2)
}
dev.off()

if( install==T) {
  
}else{
  
}



