# P.chondrillina-project-scripts

## The command lines and programming script used in this project are divided into three parts (folders):
### Part 1: Genome assembly & SNPs calling:
#### The file 'Genome assembly' contains the command lines used in Genome assembly, polishing and quality control.
#### The file 'SNPs calling' contains the command lines used in mapping, PCR duplicates removal, SNPs calling, SNPs quality filtering, SNPs counting.

### Part 2: population structure 
#### The file 'ADMIXTURE for admixture analysis' contains the command lines used in admixture analysis with software ADMIXTURE and result visualization.
#### The file 'IQTREE2 for phylogenomic tree' contains the command lines used in phylogenomic analysis with software IQTREE2.
#### The R script 'PCoA&DAPC.R' is for Principal Coordinates Analysis and Discriminant analysis of principal components with R package adegenet.

### Part 3: genome scans
#### The file 'Fst&Tajima'sD' contains command lines used to calculate Fst and Tajima'sD using VCFtools
#### The file 'Divide SNPs into contigs' contains command lines used to divide genome-wide SNPs into contigs in order to calculate pi and Dxy with the R package Popgenome.
#### The R script 'pi&Dxy.R' is for calculating pi and Dxy statistics with R package Popgenome.
#### The R script 'visualization.R' is for merging all the statistic results (pi, Dxy, Fst, Tajima'sD) into a data frame and visualizing the results.
