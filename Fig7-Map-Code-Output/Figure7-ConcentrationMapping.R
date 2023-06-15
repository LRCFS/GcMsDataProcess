###########################################################################
#
# Detection and quantitation of novel benzodiazepines in tablets, powders,
# blotters and infused materials from Scottish prisons  - Copyright (C) 2021
#
# Leverhulme Research Centre for Forensic Science

# Centre for Forensic Science, Department of Pure and Applied Chemistry,
# University of Strathclyde, Royal College, 204 George Street, Glasgow

# Victoria Marland, Robert Reid, Andrew Brandon, Kevin Hill, Fiona Cruickshanks,
# Craig McKenzie, Caitlyn Norman, Niamh Nic Daéid, Herve Menard

# Website: https://github.com/LRCFS/GcMsDataProcess
# Contact: lrc@dundee.ac.uk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
###########################################################################
#
# This general code to run first
#
###########################################################################

#Clear all lists from memory to avoid unintentional errors

rm(list=ls()) 

#############################################################
#####                 File requirement                  #####
#############################################################
# load the libraries
library(readxl)
library(xlsx)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(gplots)

##################
# HEAT MAP
##################

#Import Data as Matrix
Map <- as.matrix(read_xlsx("Fig7-ConcentrationMapping/ConcentrationMapping_Metadata.xlsx", col_names=FALSE))
Map

Map<- as.data.frame(Map[-c(1), ])

colnames(Map) <- c("A","B","C","D","E","F","G","H","I","J","K","L")
rownames(Map) <- c(1:9)

Map <- mutate_all(Map, function(x) as.numeric(as.character(x)))
Map2 <- data.matrix(Map)

#Select Colour Palette

coul = colorRampPalette(brewer.pal(9,"YlOrRd"))(6)

#Make Heat Map

png(filename="Fig7-ConcentrationMapping/Fig7-Concentration_Mapping.png",units="in",width=12,height=9,res=300)
heatmap.2(x=Map2,col=c(coul),Rowv=NULL,Colv=NULL, trace="none",density.info="none")
dev.off()