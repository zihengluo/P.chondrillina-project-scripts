# Before using this script, please divide the genome-wide SNPs into contigs,
# and place each subset SNPs vcf file into a directory. The required command lines
# can be found in the file 'Divide SNPs into contigs'.
# Please prepare a csv with two columns containing contig id and contig length.
# Also please prepare a txt file containing the population information,
# including two column with isolates names and population names. 

# load package
library(tidyverse)
library(PopGenome)

# import a data frame containing contig id and length
contigs_length = read.csv("./length854.csv")
contigs = data.frame(contigs_length)

# import the population information
pop_info <- read_delim("./pch_pops.txt", delim = "\t")
populations <- split(pop_info$ind, pop_info$pop)

#create an empty dataframe to contain results
genomeDataDF = as_tibble(data.frame(matrix(ncol=16,nrow=0)))
column.names = c("contig", "start", "stop", "mid", "CIN_pi","HCN_pi", "HisI_pi","CIN_HCN_fst", "CIN_HisI_fst","HCN_HisI_fst", "CIN_HCN_dxy", "CIN_HisI_dxy", "HCN_HisI_dxy", "CIN_Taj", "HCN_Taj", "HisI_Taj" ) 
colnames(genomeDataDF) <- column.names
snpsum = 0

for (i in 1:nrow(contigs)) {
  
  # obtain length and id of each contig
  l=contigs$length[i]
  id=contigs$Contigs[i]
  id = toString(id)
  print(id)
  
  # filter out all the contigs shorter than 30kb
  if (l < 30000) { next }
  # import the SNPs of each contig
  genomeclass = readData(id, format = "VCF", include.unknown = TRUE, FAST = TRUE)
  get.sum.data(genomeclass)
  #sum up SNPs
  snpsum = snpsum + genomeclass@n.biallelic.sites + genomeclass@n.polyallelic.sites
  
  #setup population
  genomeclass = set.populations(genomeclass, populations)
  
  #set up sliding windows, size 20kb, step 10kb.
  snpsites =  as_tibble(data.frame(get.sum.data(genomeclass)))
  contigEnd =snpsites$n.sites
  if (contigEnd < 20000) {
    contigEnd = 20001
  }
  window_size <- 20000
  window_jump <- 10000
  
  
  # specify the start, stop and mid-point of windows, 
  window_start <- seq(from = 1, to = contigEnd-window_size+1, by = window_jump)
  window_stop <- window_start + window_size
  windowsDf <- data.frame(contig = id,start = window_start, stop = window_stop, 
                          mid = window_start + (window_stop-window_start)/2)
  
  # make a sliding window dataset
  genome_sw <- sliding.window.transform(genomeclass, width = window_size, jump = window_jump, type = 2)
  
  # calculate statistics
  genome_sw = F_ST.stats(genome_sw, mode = "nucleotide")
  genome_sw = neutrality.stats(genome_sw, FAST=TRUE)
  nucDiv = genome_sw@nuc.diversity.within/window_size
  thetas = genome_sw@theta_Watterson/window_size
  Dxy = t(genome_sw@nuc.diversity.between/window_size)
  fst = t(genome_sw@nuc.F_ST.pairwise)
  Taj = genome_sw@Tajima.D
  pops <- sort(unique(pop_info$pop))
  colnames(nucDiv) <- paste0(pops, "_pi")
  colnames(Taj) = paste0(pops, '_Taj')
  colnames(thetas) = paste0(pops, '_thetas')
  
  # merge all the statistic results into a data frame
  x <- colnames(fst)
  x <- sub("pop1", pops[1], x)
  x <- sub("pop2", pops[2], x)
  x <- sub("pop3", pops[3], x)
  x <- sub("/", "_", x)
  paste0(x, "_fst")
  paste0(x, "_dxy")
  colnames(fst) <- paste0(x, "_fst")
  colnames(Dxy) <- paste0(x, "_dxy")
  contigDataDf <- as_tibble(data.frame(windowsDf, nucDiv, thetas, fst, Dxy,Taj))
  contigDataDf['contig'] = id
  genomeDataDF = rbind(genomeDataDF, contigDataDf)
  
}

#save statistic results into a csv file
write.csv(genomeDataDF,"./genomeDataDF_30kcnt.20kwd.thetas.csv", row.names = FALSE)


