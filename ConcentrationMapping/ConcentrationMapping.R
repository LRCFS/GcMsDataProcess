getwd()
rm(list=ls())
file.choose()
setwd("/Users/victoriamarland/Documents/Cardquant")
library(readxl)
library(xlsx)
library(plyr)

#============
# HEAT MAP
#============

#Import Data as Matrix
Map <- as.matrix(read_xlsx("ConcentrationMapping_Metadata.xlsx", col_names=FALSE)) ; Map

Map<- as.data.frame(Map[-c(1), ])

colnames(Map) <- c("A","B","C","D","E","F","G","H","I","J","K","L")
rownames(Map) <- c(1:9)

Map <- mutate_all(Map, function(x) as.numeric(as.character(x)))
Map2 <- data.matrix(Map)

#Select Colour Palette
library(RColorBrewer)
coul = colorRampPalette(brewer.pal(9,"YlOrRd"))(6)

#Make Heat Map
library("gplots")
png(filename="Concentration_Mapping.png",units="in",width=12,height=9,res=300)
heatmap.2(x=Map2,col=c(coul),Rowv=NULL,Colv=NULL, trace="none",density.info="none")
dev.off()

#If you want to explore possible colour palettes
display.brewer.all()
