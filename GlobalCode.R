
# To clean the Global environment
rm(list=ls()) 

#############################################################
#####                 File requirement                  #####
#############################################################
# The files to be imported are generated from Scopus and Web Of Science databases
# The columns will need to contain:
# Year; Title; Source.title; Authors; AuthorID; DE; DES; EID; SO; DT
library(readxl)
library(xlsx)
library(plyr)
library(dplyr)
library(stringr)
# library(stringr)
library(tidyr)
library(hablar)
library(ggplot2)
# library(maps)
# library(countrycode)
# library(RColorBrewer)
# library(bibliometrix)
# library(tidyverse)
# library(plotly)
# library(gridExtra)
# library(extrafont)
library("scales")
library(zoo)
library(hyperSpec)


#############################################################
#####                      Function                     #####
#############################################################

# source("Functions/SearchAndReplace.R")
# source("Functions/Diacritics.R")
# source("Functions/PartialMatch.R")
# source("Functions/RemoveDuplicate.R")


#############################################################
#####                Folder & Files                     #####
#############################################################

# set extension and Citation
# extensionTXT <- ".TXT"
# extensionBIB <- ".bib"
extensionCSV <- ".csv"
extensionXLSX <- ".xlsx"
extensionMS1 <- ".ms1"

# where the generated figures are saved, create folder if not existing
Metadata.dir <- "Metadata/"
GcData.dir <- "GcData/"
HeatMapData <- "HeatMap/"
TempData.dir <- "TempData/"
GcDataConverterMs.dir <- "GcDataConverterMs/"

Results.dir <- "Results/"
dir.create(file.path(Results.dir),recursive = TRUE) # will create folder if not already there.
#

HeatResults.dir <- "HeatMapResults/"
dir.create(file.path(Results.dir),recursive = TRUE) # will create folder if not already there.

#############################################################
#####           Constants and thresholds                #####
#############################################################

threshold <- 0.2 # % peak height of max intensity TIC

#############################################################
#####                       Codes                       #####
#############################################################

# This codes can be run subsequently or independently as each create necessary outputs for the next codes.

# This script changes the converted *.D Agilent data using Proteo Wizard - MsConvert
# The converted data file format is *.ms1 and placed in the GcDataConverterMs folder
# The Msconverter can handle bath process, the R code will run one file after another and save the peak area results in GcData folder 
source("Code/AgilentDataHeatmap.R")

# This script use the results of the previous code (or self entered) in GcData and combined it to the metadata
#source("Code/Metadata.R")

# This script produce a heatmap of concentration for a given sample
#source("Code/SampleHeatmap.R")


