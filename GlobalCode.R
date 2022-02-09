
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
library(tidyverse)
# library(plotly)
# library(gridExtra)
# library(extrafont)
library("scales")
library(zoo)
library(hyperSpec)
library(reshape2)
library(baseline)
library(smooth)
library(IDPmisc)

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
MetadataOutput.dir <- "Metadata/FiguresOutput/"

GcData.dir <- "GcData/"
HeatMapData <- "HeatMap/"
TempData.dir <- "TempData/"
# This is where the MS1 files converted using ProteoWizard - MSConvert are placed  
GcDataConverterMs.dir <- "GcDataConverterMs/"
# This is where the MS1 files are saved after their content is reordered in specific columns and rows  
GcDataConvertedRcode.dir <- "GcDataConvertedRcode/"
dir.create(file.path(GcDataConvertedRcode.dir),recursive = TRUE) # will create folder if not already there.

Results.dir <- "Results/"
dir.create(file.path(Results.dir),recursive = TRUE) # will create folder if not already there.

Backup.dir <- "Results/Backup/"

Library.dir <- "Results/Library/"
dir.create(file.path(Library.dir),recursive = TRUE) # will create folder if not already there.

LibraryFound.dir <- "Results/LibraryFound/"
dir.create(file.path(LibraryFound.dir),recursive = TRUE) # will create folder if not already there.

LibraryNotFound.dir <- "Results/LibraryNotFound/"
dir.create(file.path(LibraryNotFound.dir),recursive = TRUE) # will create folder if not already there.

LibraryComparison.dir <- "Results/LibraryComparison/"
dir.create(file.path(LibraryComparison.dir),recursive = TRUE) # will create folder if not already there.


HeatResults.dir <- "HeatMapResults/"
dir.create(file.path(Results.dir),recursive = TRUE) # will create folder if not already there.

#############################################################
#####           Constants and thresholds                #####
#############################################################

threshold <- 0.2 # % peak height of max intensity TIC
angular.vertor <- 0.7
# m/z accepted precision
Mass.precision <- 0.2

#smooth.loop <- 10 # number of smmothing repeat applied to the data
rolling.average <- 3 # value for rolling average
IS.Exp.Range <- 390  # the retention time Internal Standard is expected at
Etizolam.Exp.Range <- 612  # the retention time Etizolam is expected at

#############################################################
#####                       Codes                       #####
#############################################################

# This codes can be run subsequently or independently as each create necessary outputs for the next codes.

# This script changes the converted *.D Agilent data using Proteo Wizard - MsConvert
# The converted data file format is *.ms1 and placed in the GcDataConverterMs folder
# The Msconverter can handle bath process, the R code will run one file after another and save the peak area results in GcData folder

# This only need to be done once per ms1 files and export is saved to a new folder: GcDataConvertedRcode
source("Code/MsFilesReorganiser.R")

# This code takes the files reordered from the previous code and calculate the areas for the relevant peaks 
source("Code/AgilentDataHeatmap.R")

# This script use the results of the previous code (or self entered) in GcData and combined it to the metadata
source("Code/Metadata.R")

# This script produce a heatmap of concentration for a given sample
#source("Code/SampleHeatmap.R")


