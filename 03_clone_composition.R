load("Results/SNVS.clonality.Rdata")
library(plyr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
color_clone_composition=c( 'M' = rgb(0,162,205, maxColorValue = 255), 
                           'B' = rgb(243,130,153, maxColorValue = 255), 
                           'P' = rgb(223,207,0, maxColorValue = 255))
theme_cloneR=function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.background = element_blank(),
          panel.border     = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}



upper=0.80
lower=0.35
snvs$upper=upper
snvs$lower=lower

y = ddply(snvs, .(Patient), summarise,
                n = length(clonality),
                n_monoclonal = sum(clonality>=unique(upper)),
                n_biclonal   = sum(clonality<unique(upper) & clonality>=unique(lower)),
                n_polyclonal = sum(clonality<unique(lower)))


y$monoclonal = y$n_monoclonal / y$n
y$biclonal   = y$n_biclonal / y$n
y$polyclonal = y$n_polyclonal / y$n
    
code = c("M","B","P"); names(code)=c('monoclonal','biclonal','polyclonal')
y$composition = code[names(which.max(y[,c('monoclonal','biclonal','polyclonal')]))]

clone.composition.plot= function(x, cl=color_clone_composition){
  names(cl)=c("monoclonal",'biclonal','polyclonal')
  x$composition = factor(x$composition, levels=c("M","B","P"))
  
  if(!is.null(x)){
    mapper = as.list(x$composition)
    names(mapper) = x$Patient
    map_labeller <- function(variable,value){
      return(mapper[value])
    }
    
    m = melt(x[,c('Patient','polyclonal','biclonal','monoclonal')], id.vars = c("Patient")) #1,4:2
    colnames(m)[2] = 'composition'
    p=ggplot(m, aes(x=Patient,y=value,fill=composition))+
      geom_bar(width=0.5, stat="identity")+
      geom_segment(aes(x=-Inf,xend=-Inf,y=0,yend=1),col="black")+
      ylab("Alterations (%)")+xlab("")+
      scale_fill_manual(values=cl,
                        guide = guide_legend(title = NULL),
                        labels=c("Clonality<35%",'35%<Clonality<80%','Clonality>80%'))+
      scale_y_continuous(labels=c("0","25","50","75","100"))+
      theme_cloneR()+theme(panel.background = element_blank(),
                           legend.position  = "top",
                           legend.key=element_rect(size=1.5, color='white'),
                           legend.text      = element_text(color="black",size=10),
                           axis.text        = element_text(color="black",size=10),
                           axis.text.y      = element_blank(),
                           axis.ticks.y     = element_blank(),
                           strip.background = element_rect(fill=NA, colour=NA),
                           strip.text.y     = element_text(angle=0,size=12, colour="black")
      )+
      coord_equal(1/0.1)+
      coord_flip()+
      geom_bar(width=0.5, stat="identity",color="black", show.legend=FALSE)
    
    return(p)
  } else{
    return(NULL)
  }
}

pdf(file="Results/CloneComposition.pdf",h=2, w=6)
clone.composition.plot(y)
dev.off()
