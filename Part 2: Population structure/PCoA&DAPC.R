# This script was modified from Rita Tam's script for PCoA and DAPC analysis

# load required packages
library("adegenet")
library("vcfR")
library("wordcloud")

# specify the path to SNPs' vcf file 
vcf_path = '/g/data/xf3/zl4459/SNPs/pchonmitorepeatrm_covlarger3.merged.freebayes.vcf'

# import SNPs data
vcf = read.vcfR(vcf_path)

# store SNPs into genlight object
genlight = vcfR2genlight(vcf)

# set working directory and checkpoint directory
setwd('/scratch/xf3/zl4459/pca')
checkpoint_path= '/scratch/xf3/zl4459/pca/checkpoint'

# specify the na values to be deleted
toRemove <- is.na(glMean(genlight, alleleAsUnit=FALSE))

# store the data without na values
genlightclean <- genlight[, !toRemove]
save(genlight, file = file.path(checkpoint_path, 'genlight.Rdata'))

# calculate the distance between individuals
dist=dist(tab(genlightclean))

# specify the first twp Principle Coordinates (PCos)
pco1 = dudi.pco(dist, scannf=FALSE, nf=2)
pco2 = dudi.pco(dist, scannf=FALSE, nf=1)

# store and load all the PCos data
save(dist, file=file.path(checkpoint_path, 'fullpco.Rdata'))
load(file.path(checkpoint_path, 'fullpco.Rdata'))

# plotting
s.label(pco1$li*1.2, clab=0, pch="", sub="P. chondrillina SNPs", csub=1.5, cpoint="")
textplot(pco1$li[,1], pco1$li[,2], words=rownames(pco1$li), cex=1.1, new=FALSE, xpd=TRUE)
add.scatter.eig(pco2$eig, nf=2, xax=1, yax=2, posi="topleft")

# calculates the proportion of overall variance occupied by the first two PCos
sum(pco2$eig[1:2])/sum(pco2$eig)

# transfer SNPs data into principle components through dimensional reduction
full.pca <- glPca(genlightclean, scale=FALSE, nf=40, parallel=require("parallel"), n.cores=48, useC=FALSE)
save(full.pca, file=file.path(checkpoint_path, "full.pca.Rdata"))
load(file.path(checkpoint_path, "genlightclean.Rdata"))
load(file.path(checkpoint_path, "full.pca.Rdata"))

# find optimal cluster number depending on BIC values
find.clusters.genlight(genlightclean, max.n.clust=14, glPca=full.pca, n.pca=2)

# create a function for DAPC analysis
DAPC <- function(gl_clean, input_pca, retained_pca, retained_da, n_clust) {
  grp <- find.clusters.genlight(gl_clean, glPca=input_pca, n.pca=retained_pca, n.clust=n_clust)
  dapc <- dapc(gl_clean, pop=grp$grp, glPca=input_pca, n.pca=retained_pca, n.da=retained_da)
  print(grp$grp)
  typeof(grp$grp)
  return(dapc)
}

# create a function to plot the cumulative variance (%) against the PC number
myInset <- function(dapc){
  temp <- dapc$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated variance (%)",
       cex=1, pch=20, type="h", lwd=2)
}

# detect the classification of isolates in the subpopulations
full.dapc.clust3 <- DAPC(genlightclean, full.pca, retained_pca=3, retained_da=2, n_clust=3)
full.dapc.clust3

# plotting of DAPC results with three clusters
myCol <- c('brown1', 'darkblue', 'chocolate1')

scatter(full.dapc.clust3, posi.da="topright", bg="white",
        col=myCol, scree.pca=FALSE)

add.scatter(myInset(full.dapc.clust3), posi="topleft", inset=c(-0.00,-0.008), ratio=0.15, bg=transp("white"))