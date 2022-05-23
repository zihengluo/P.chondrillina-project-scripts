library(tidyverse)
library(dplyr)
library(magrittr)
library(ggpubr)

# import statistics results obtained from Popgenome and add a column of window id
genomeDataDF =data.frame(read.csv("./genomeDataDF_30kcnt.20kwd.csv"))
id = seq(1,8023)
genomeDataDF = cbind(genomeDataDF, id)

# import Fst obtained from VCFtools
fst12 <- data.frame(read.delim("pop12_fst_20kwd.tab", header=F))[,c(2,3,4,5)]
names(fst12) <- c("contig", "start", "stop", "pop12avg_fst" )
fst12[,2]=fst12[,2]+1
fst12[,3]=fst12[,3]+1
id1 = seq(9397)
fst12 = cbind(fst12, id1)
fst12['pop12avg_fst'] = as.numeric(unlist(fst12['pop12avg_fst']))
#view(fst)

fst13 <- data.frame(read.delim("pop13_fst_20kwd.tab", header=F))[,c(2,3,4,5)]
names(fst13) <- c("contig", "start", "stop", "pop13avg_fst" )
fst13[,2]=fst13[,2]+1
fst13[,3]=fst13[,3]+1
id1 = seq(9406)
fst13 = cbind(fst13, id1)
fst13['pop13avg_fst'] = as.numeric(unlist(fst13['pop13avg_fst']))

fst23 <- data.frame(read.delim("pop23_fst_20kwd.tab", header=F))[,c(2,3,4,5)]
names(fst23) <- c("contig", "start", "stop", "pop23avg_fst" )
fst23[,2]=fst23[,2]+1
fst23[,3]=fst23[,3]+1
id1 = seq(9407)
fst23 = cbind(fst23, id1)
fst23[fst23 == '.'] = 0
fst23['pop23avg_fst'] = as.numeric(unlist(fst23['pop23avg_fst']))

# merge the three data frames of Fst
df = merge(fst12,genomeDataDF, by = c('contig','start','stop'))
df1 = merge(df, fst13, by = c('contig','start','stop'))
df2 = merge(df1, fst23, by = c('contig','start','stop'))

# import Tajima's D obtained from VCFtools
tajima1 <- data.frame(read.delim("pop1.Tajima.D.tab", header=T))[,c(1,2,4)]
names(tajima1) <- c("contig", "start",  "pop1_tajimasd" )
tajima1['pop1_tajimasd'] = as.numeric(unlist(tajima1['pop1_tajimasd']))
tajima1[,2]=tajima1[,2]+1


tajima2 <- data.frame(read.delim("pop2.Tajima.D.tab", header=T))[,c(1,2,4)]
names(tajima2) <- c("contig", "start",  "pop2_tajimasd" )
tajima2['pop2_tajimasd'] = as.numeric(unlist(tajima2['pop2_tajimasd']))
tajima2[,2]=tajima2[,2]+1


tajima3 <- data.frame(read.delim("pop3.Tajima.D.tab", header=T))[,c(1,2,4)]
names(tajima3) <- c("contig", "start",  "pop3_tajimasd" )
tajima3['pop3_tajimasd'] = as.numeric(unlist(tajima3['pop3_tajimasd']))
tajima3[,2]=tajima3[,2]+1

# merge all the statistics (pi, dxy, Fst, Tajima's D) together into a data frame
tajima12 = merge(tajima1,tajima2, by = c('contig','start'))
d123 = merge(tajima12,tajima3, by = c('contig','start'))
df3 = df2[,c(1,2,3,19,4,6,7,8,9,13,14,15,25,27)]
df4 = merge(df3, d123, by = c('contig','start'), all.x = TRUE)
df4 = df4[,c(4,1,2,3,6,9,8,7,12,11,10,5,13,14,15,16,17)]
colnames(df4)[seq(6,17)] = c('pop1_pi', 'pop2_pi', 'pop3_pi', 'pop12_dxy', 'pop13_dxy', 
                             'pop23_dxy','pop12_fst', 'pop13_fst','pop23_fst',
                             'pop1_tajimasd','pop2_tajimasd','pop3_tajimasd')

# Plot Tajima's D
tajplot = df4%>% select(id, pop1_tajimasd,pop2_tajimasd, pop3_tajimasd)
# use gather to rearrange everything
alltaj = gather(tajplot,-id, key = "stat", value = "value")

#let's rename the columns to make it easier for axis labelling the in the facet plot.
# first make a factor
x = factor(alltaj$stat)
# then reorder the levels
x = factor(x, levels(x)[c(1,2,3)])
# add to data.frame
alltaj$stat <- x
# construct a plot with facets
HVsSplot = ggplot(alltaj, aes(id, value, colour = stat)) + geom_point(size =0.3) 
HVsSplot = HVsSplot + facet_grid(stat~., scales = "free_y")
HVsSplot = HVsSplot + xlab("window_id")
HVsSplot + theme_light() + theme(legend.position = "none")
ggplot(alltaj, aes(x = stat, y = value, colour = stat)) +
  geom_boxplot() + theme_light() + xlab('')
box = ggplot(alltaj, aes(x = stat, y = value, colour = stat)) +
  geom_boxplot() + theme_light()
figure <- ggarrange(HVsSplot, box, 
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)
figure

# Plot Fst
x =ggplot(df4, aes(x= pop12_fst)) +
  geom_histogram(fill ='red',binwidth = 0.01)+theme_light() + xlab('pop1vspop2_Fst')  +ylab('Count')+
  geom_vline(aes(xintercept = mean(pop12_fst)),col='black',size=0.5)+
  annotate("text",x = 0.5 ,y = 200 ,label = paste("Mean =",0.37 ),
           col = "black",size = 5)
y = ggplot(df4, aes(x= pop13_fst)) +
  geom_histogram(fill ='red',binwidth = 0.01)+theme_light()+ xlab('pop1vspop3_Fst')+
  geom_vline(aes(xintercept = mean(pop13_fst)),col='black',size=0.5)+
  annotate("text",x = 0.4 ,y = 200 ,label = paste("Mean =",0.25 ),
           col = "black",size = 5)+ylab('Count')
z = ggplot(df4, aes(x= pop23_fst)) +
  geom_histogram(fill ='red',binwidth = 0.01)+theme_light()+ xlab('pop2vspop3_Fst')+
  geom_vline(aes(xintercept = mean(pop23_fst)),col='black',size=0.5)+
  annotate("text",x = 0.5 ,y = 200 ,label = paste("Mean =",0.37),
           col = "black",size = 5)+ylab('Count')
xyz <- ggarrange(x,y,z,
                 labels = c("", "", ''),
                 ncol = 1, nrow = 3)
xyz
