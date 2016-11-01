load("Input/ASCAT.Rdata")

library(sciClone)

# 1. chromosome 2. position 3. reference-supporting read counts 4. variant-supporting read counts 5. variant allele fraction (between 0-100)

vafs = snvs[,c('Chromosome','Position','ref_counts','var_counts','obs.VAF')]
vafs$obs.VAF = vafs$obs.VAF*100

# 1. chromosome 2. segment start position 3. segment stop position 4. copy number value for that segment. Unrepresented regions are assumed to have a copy number of 2.
cnv = ascat[,c('Chr','Start','End','cn')]

lib = installed.packages()
if("sciClone" %in% rownames(lib)){
sc = sciClone(vafs=vafs,
              copyNumberCalls=cnv,
              sampleNames=snvs$Patient,
              minimumDepth = 50,
              useSexChrs=F,
              maximumClusters = 4)
#save(sc, file="Input/sciClone.Rdata")
}else{
  load("Input/sciClone.Rdata")
}
writeClusterTable(sc, "Results/clones.xls")
tmp = read.delim("Results/clones.xls")


x = paste0(snvs$Chromosome,".",snvs$Position)
y = paste0(tmp[,1],".",tmp[,2])
snvs$cluster =  tmp$cluster[match(x,y)]

stats = ddply(snvs, .(cluster), function(x) summary(x$clonality) )

stats$events=NA
stats$events[1] = "EGFR amplification"
stats$events[2] = "PIK3CA chr3:178952086:T>G c=0.72"
stats$events[3] = "PIK3CA chr3:178916876:G>A c=0.56; KRAS chr12:25398255 G>T c=0.47"

write.csv(stats,"Results/Clones.stats.csv", row.names = F)
