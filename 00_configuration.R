pkgs = c('RColorBrewer','plyr','ggplot2','reshape2','ggrepel','devtools')

lib = installed.packages()
installed=pkgs %in% rownames(lib)
not_installed = which(!installed)
if(length(not_installed)>0){
  for(i in not_installed) install.packages(pkgs[i])
}

source("http://bioconductor.org/biocLite.R")
biocLite("IRanges")
#install devtools if you don't have it already
install.packages("devtools")
library(devtools)
install_github("genome/bmm")
install_github("genome/sciClone")