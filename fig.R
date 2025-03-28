library(dplyr)
library(tidyverse)
library(patchwork)
library(devtools)
library(circlize)

setwd("~/Documents/IARC_HPV group/Figs/")

####


###
circos.info()

cytoband.file = system.file(package = "circlize", "extdata", "cytoBandIdeo.txt")
cytoband.file = system.file("cytoBandIdeo.txt")

setwd("~/Documents/")
cytoband.df = read.table("cytoBandIdeo.txt", colClasses = c("character", "numeric",
                                                            "numeric", "character", "character"), sep = "\t")


cytoband.df = subset(cytoband.df,cytoband.df$V1 %in% c("chr1" , "chr2",  "chr3" , "chr4",  "chr5",  "chr6" , "chr7",  "chr8",
                                                       "chr9"  ,"chr10" ,"chr11","chr12" ,"chr13" ,"chr14", "chr15", "chr16" ,
                                                       "chr17" ,"chr18" ,"chr19" ,"chr20", "chr21" ,"chr22","chrX"))
cytoband.df[852,] = c("HPV",1,100000000,"HPV","gneg")

chromosome <- c("chr1" , "chr2",  "chr3" , "chr4",  "chr5",  "chr6" , "chr7",  "chr8",
                "chr9"  ,"chr10" ,"chr11","chr12" ,"chr13" ,"chr14", "chr15", "chr16",
                "chr17" ,"chr18" ,"chr19" ,"chr20", "chr21" ,"chr22","chrX","HPV")
cytoband.df[,2:3] = apply(cytoband.df[,2:3],2,as.numeric)
### HIVID info plot
###
circos.clear()
circos.par("start.degree" = 90) ### 
circos.initializeWithIdeogram(cytoband.df)
options(repr.plot.width = 30,repr.plot.height = 8)
###
library(rtracklayer)

plot_data = read.csv("pHR_HPV_integration_breakpoints.csv")
plot_data = plot_data[,c(2,4,5,6,8)]
plot_data_CG = read.csv("pHR_HPV_integration_Cell_genomics.csv")
plot_data[110:116,] = plot_data_CG[,c(2,3,5,6,4)]
 
plot_data$Viral_Type = gsub("^HpV","HPV",plot_data$Viral_Type)
plot_data$HPV_Position = plot_data$HPV_Position*10000

#####HPV34
HPV34 = subset(plot_data,Viral_Type %in% c("HPV34"))

HPV_mid = HPV34[,c(2,5)];HPV_mid$Viral_Type = "HPV"
Human_mid = HPV34[,c(3,4)]

circos.genomicLink(Human_mid, HPV_mid,col = "#588157")

#####
HPV66 = subset(plot_data,Viral_Type %in% c("HPV66"))

HPV_mid = HPV66[,c(2,5)];HPV_mid$Viral_Type = "HPV"
Human_mid = HPV66[,c(3,4)]

circos.genomicLink(Human_mid, HPV_mid,col = "#4361ee")

#####
HPV73 = subset(plot_data,Viral_Type %in% c("HPV73"))

HPV_mid = HPV73[,c(2,5)];HPV_mid$Viral_Type = "HPV"
Human_mid = HPV73[,c(3,4)]

circos.genomicLink(Human_mid, HPV_mid,col = "#ffc857")

#####
HPV82 = subset(plot_data,Viral_Type %in% c("HPV82"))

HPV_mid = HPV82[,c(2,5)];HPV_mid$Viral_Type = "HPV"
Human_mid = HPV82[,c(3,4)]

circos.genomicLink(Human_mid, HPV_mid,col = "#fb6f92")

###HPV genes comparison
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggsignif)
library(reshape2)

p <-ggplot(data_plot, aes(x=variable, y=value),fill = Type) +
  geom_boxplot(alpha = 1, outlier.shape = NA) + 
 geom_point(mapping = aes(x=variable, y=value,
                         color = Type),size = 3, shape = 16) +
  scale_color_manual(values = c("#123f81","#cf4728","#81b29a"),
                     name = "Spe") +
#geom_hline(aes(yintercept = 2.532410), colour="#7a6a56", linetype="dashed") +
theme(legend.position = 'right',
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 30, hjust = 1)) +
  theme(panel.grid.major = element_blank(), #????????
        panel.grid.minor = element_blank(), #????????
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"), #?߿?
       
        panel.background = element_rect(fill = 'white'), #????ɫ
        plot.background= element_rect(fill="white"))+
  geom_signif(comparisons = list(c("E1", "E6"),c("E1", "E7"),c("E2", "E7")), 
              map_signif_level = TRUE)

options(repr.plot.width = 12, repr.plot.height = 6.5)
p
