#!/usr/bin/env Rscript
#.libPaths(c(.libPaths(), "/net/eichler/vol2/eee_shared/modules/anaconda3/envs/py3.20161130/lib/R/library"))
library(ggplot2)
library(plyr)
require(gridExtra)
require(scales)
library(RColorBrewer)
library(reshape2)
#install.packages("tidyr")
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
#install.packages("argparse")
library(argparse)
library(splitstackshape)
#suppressPackageStartupMessages(library("argparse"))
theme_set(theme_classic(base_size = 18))

plotextra = NULL
for( myf in Sys.glob("~/Desktop/chm13_hifi_assembly/recalled_canu_asm/results/*per.tbl")){
  print(myf)
  plotextra = rbind(plotextra, fread(myf) ) 
}
colnames(plotextra) = c("Extra", "PercentResolved", "Assembly")
plotextra = plotextra[plotextra$Extra%%500 == 0]

old = c( 
  "GRCh38" ,
  "AK1"    ,
  "HX1",
  "NA12878_ONT",  
  "CHM1" ,
  "recalled_canu1.7_2xracon",
  "chm13_clr_quiver_pilon"   
  )


ID=c(
  "GRCh38 (BAC, Hierarchical)",
  "AK1 (CLR, FALCON)",
  "HX1 (CLR, FALCON)",
  "NA12878 (ONT, Canu)",
  "CHM1 (CLR, FALCON)",
  "CHM13 (HiFi, Canu)", 
  "CHM13 (CLR, FALCON)"
)
i = 1
for(pre in old){
  new = ID[i]
  plotextra$Assembly[plotextra$Assembly == pre ] = new
  i = i + 1
}

color = c(
  "#000000",
  "#808080",
  "#999999", 
  "#b2b2b2",
  "#cccccc",
  "#b20000", 
  "#b20000"
)
names(color) = ID

myxlab = "Minimal extension from duplication boundary (kbp)"

p = ggplot(plotextra) + 
  geom_point(aes(Extra/1000, PercentResolved, group=Assembly, color=Assembly, shape=Assembly), size=3, alpha =1) +
  geom_line(aes(Extra/1000, PercentResolved, group=Assembly, color=Assembly), size = .75) +
  scale_color_manual(values=color) +
  ylab("% of segmental duplications resolved\nin each genome assembly") +
  theme(legend.position="top") +
  xlab(myxlab) +  scale_x_continuous(labels = comma);p


ggsave("~/Desktop/chm13_hifi_assembly/recalled_canu_asm/results/SDvsRef.pdf", plot = p, width = 16, height = 9)



hifi= "CHM13 (HiFi, Canu)"
clr =  "CHM13 (CLR, FALCON)"
ont =  "NA12878 (ONT, Canu)"

diff = (plotextra[plotextra$Assembly==hifi]$PercentResolved - plotextra[plotextra$Assembly==clr]$PercentResolved) /plotextra[plotextra$Assembly==clr]$PercentResolved
mean(diff)
mean(plotextra[plotextra$Assembly==hifi]$PercentResolved)
mean(plotextra[plotextra$Assembly==clr ]$PercentResolved)
mean(plotextra[plotextra$Assembly==ont ]$PercentResolved)

diff = (plotextra[plotextra$Assembly==hifi]$PercentResolved - plotextra[plotextra$Assembly==ont]$PercentResolved) /plotextra[plotextra$Assembly==ont]$PercentResolved
mean(diff)
