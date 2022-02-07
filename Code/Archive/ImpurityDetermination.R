#
# Prior running, GC-MS data must be converted using Proteo Wizard/ MsConvert
# 
#####################################################################
#####              Run Global Code in first instance            #####
#####################################################################
# Load Metadata files
filenameGcDataConverterMs <- list.files(GcDataConverterMs.dir, extensionMS1, full.names=TRUE)

for (file in filenameGcDataConverterMs) {
DataName <- gsub(extensionMS1, "", file)

# for testing code on a single file (first on list) in filenameGcDataConverterMs
# DataName <- gsub(extensionMS1, "", filenameGcDataConverterMs[2])

DataName <- gsub(".*/", "", DataName)

#####################################################################
#####  re-organise converted data from  Proteo Wizard MsConvert #####
#####################################################################
# for testing code on a single file (first on list) in filenameGcDataConverterMs
# GcDataConverter <- read.csv2(filenameGcDataConverterMs[2], sep = "\t", header = FALSE)

GcDataConverter <- read.csv2(file, sep = "\t", header = FALSE)

temp_data_0 <- strsplit(as.character(GcDataConverter$V1),' ') 
#do.call(rbind, temp_data_0)

temp_data_1 <-data.frame(GcDataConverter,do.call(rbind,temp_data_0))
# remove rows containing "H"
temp_data_1 <- temp_data_1 %>% 
  filter(!str_detect(V1, 'H'))

# remove rows containing "NativeID"  
temp_data_1 <- temp_data_1 %>% 
  filter(!str_detect(V2, 'NativeID'))

# remove rows containing "BPI"
temp_data_1 <- temp_data_1 %>% 
  filter(!str_detect(V2, 'BPI'))

# remove rows containing "BPM"
temp_data_1 <- temp_data_1 %>% 
  filter(!str_detect(V2, 'BPM'))

# replace blank with "NA"
temp_data_1$V3 <- as.numeric(temp_data_1$V3)

# fill cells with previous values in the same column
temp_data_1 <- temp_data_1 %>%
  fill(V3)
# remove rows containing "TIC"
temp_data_1 <- temp_data_1 %>% 
  filter(!str_detect(temp_data_1$V2, "TIC"))

# replace blanck with "NA"
temp_data_1$V2[temp_data_1$V2==""]<-0

n <- nrow(temp_data_1)

# fill cells with values taken in different column
for(i in 1:nrow(temp_data_1)) { # for-loop over rows
  m <- temp_data_1[i,2]
  if (m =="RTime") {
    Rtime <- temp_data_1[i,3]
    i <- i+1
    q <- temp_data_1[i,2]
    
    while (q == "0" & i <= n ) {
      temp_data_1[i,2] <- Rtime
      i <- i+1
      q <- temp_data_1[i,2]
      }
  }
  i <- i+1
}

# remove rows containing "S"
temp_data_1 <- temp_data_1 %>% 
  filter(!str_detect(V1, 'S'))

# remove rows containing "S"
temp_data_1 <- temp_data_1 %>% 
  filter(!str_detect(V1, 'I'))

DataHeatMap <- temp_data_1 %>%
  select(-V1)

names(DataHeatMap)[1] <- "RetentionTime"
names(DataHeatMap)[2] <- "TIC"
names(DataHeatMap)[3] <- "Mass"
names(DataHeatMap)[4] <- "Intensity"

# convert all columns to numeric
DataHeatMap <- mutate_all(DataHeatMap, function(x) as.numeric(as.character(x)))

# convert the retention time to sec
DataHeatMap$RetentionTime <- DataHeatMap$RetentionTime *60

# plot figure including all data
p <- ggplot(DataHeatMap, aes(RetentionTime, Mass))+
   geom_tile(aes(fill = Intensity))
 show(p)

#########################################################
####### filter data according to threshold value  #######
#########################################################
# # determine mass peak intensity per group
# DataHeatMap <- DataHeatMap %>%
#   group_by(RTime) %>%
#   mutate(freq = Intensity / sum(Intensity))
# 
# # scale TIC intensity 0 to 1 
# DataHeatMap$TICfreqScaled <- rescale(DataHeatMap$TIC)
# 
# # filter values higher than threshold
# DataHeatMapThreshold <- DataHeatMap %>%
#   filter(TICfreqScaled >= threshold)
# 
# #plot figure including all data
# p <- ggplot(DataHeatMapThreshold, aes(RTime, Mass))+
#   geom_tile(aes(fill = Intensity))
# show(p)
# 
#########################################################
#######           calculate TIC peak areas        #######
#########################################################

DataTicRentention <- DataHeatMap %>%
  select(RetentionTime, TIC) %>%
  distinct()

# plot figure including all data
 p <- ggplot(DataTicRentention, aes(RetentionTime, TIC))+
   geom_line()+
   ylim(0,250000)
 show(p)

########################################################
#######    baseline corrections using hyperSpec  #######
########################################################

# baseline corrections adapted from and using the library package hyperSpec
# for more information: https://cran.r-project.org/web/packages/hyperSpec/index.html
# 
# the data must be in hyperspectral format (this could be extended to creating matrix of data)

# convert data to hyperSpec format
fileTIC <- as.data.frame(DataTicRentention)

 
tempFileDataTIC <- fileTIC$TIC

tempFileDataTIC <- as.data.frame(t(as.data.frame(tempFileDataTIC)))

tempFileDataTICmatrix <-data.matrix(tempFileDataTIC)

# bc.modpolyfit <- baseline(tempFileDataTICmatrix[1,, drop=FALSE], method='modpolyfit', deg=4)
# 
# plot(bc.modpolyfit)
# 
# 
# # extract the corrected spectra
# tempFileCorrected <- bc.modpolyfit@corrected
# # convert to dataframe and transpose to 1 column
# tempFileCorrected <- as.data.frame(t(as.data.frame(tempFileCorrected)))
# # rename column to Corrected
# names(tempFileCorrected)[1] <- "CorrectedTrend"
# 
# 
# # extract the baseline
# tempFileBaseline <- bc.modpolyfit@baseline
# # convert to dataframe and transpose to 1 column
# tempFileBaseline <- as.data.frame(t(as.data.frame(tempFileBaseline)))
# # rename column to Corrected
# names(tempFileBaseline)[1] <- "BaselineTrend"
# 
# Combined <- cbind(fileTIC,tempFileCorrected, tempFileBaseline)
# names(Combined)[1] <- "RetentionTime"
# 
# Combined$Subtracted <- Combined$TIC-Combined$BaselineTrend
# 
# CombinedNarrow <- Combined %>%
#   select(RetentionTime,Subtracted) %>%
#   mutate(srate_ma01 = rollmean(Subtracted, k = 5, fill = NA))
# 
# # Plot it
# p <- ggplot(CombinedNarrow, aes(x=RetentionTime)) +
#   geom_line(aes(y = Subtracted, colour = "subtracted")) + 
#   geom_line(aes(y = srate_ma01, colour = "average")) +
#   ylim(0,150000) +
#   xlim(400,600)
#   
# show(p)
#   


# Description
# A translation from Kevin R. Coombes et al.’s MATLAB code for detecting peaks and removing
# baselines
# https://cran.r-project.org/web/packages/baseline/baseline.pdf

baseline.peakDetection <- baseline(tempFileDataTICmatrix[1,, drop=FALSE], method='peakDetection',
                             left=200, right=200, lwin=5, rwin=5)

# extract the corrected spectra
tempFileCorrected <- baseline.peakDetection@corrected

# convert to dataframe and transpose to 1 column
tempFileCorrected <- as.data.frame(t(as.data.frame(tempFileCorrected)))

# rename column to Corrected
names(tempFileCorrected)[1] <- "CorrectedTrend"

# extract the baseline
tempFileBaseline <- baseline.peakDetection@baseline

# convert to dataframe and transpose to 1 column
tempFileBaseline <- as.data.frame(t(as.data.frame(tempFileBaseline)))

# rename column to BaselineTrend
names(tempFileBaseline)[1] <- "BaselineTrend"

Combined.bc.fillPeak <- cbind(fileTIC,tempFileCorrected, tempFileBaseline)
names(Combined.bc.fillPeak)[1] <- "RetentionTime"

Combined.bc.fillPeak$Subtracted <- Combined.bc.fillPeak$TIC-Combined.bc.fillPeak$BaselineTrend

# Plot it
p <- ggplot(Combined.bc.fillPeak, aes(x=RetentionTime)) +
  geom_line(aes(y = Subtracted, colour = "subtracted")) +
  geom_line(aes(y = BaselineTrend, colour = "baseline")) +
  ylim(0,25000) #+
#  xlim(200,800)

show(p)


ggsave(
  sprintf("%s_BL.tiff",DataName),
  plot = p,
  device = NULL,
  path = file.path(LibraryComparison.dir),
  scale = 1,
  width = 7.5,
  height = 5.0,
  units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = TRUE
)


# Finding Peaks in Raw Data
# Description
# Returns position, signal height and approximate width at half maximum peak height.
# Adapted from https://www.imsbio.co.jp/RGM/R_rdfile?f=IDPmisc/man/peaks.Rd&d=R_CC

baseline.peakDetection.Narrow <- Combined.bc.fillPeak %>%
  select(RetentionTime,Subtracted) %>%
  mutate(srate_ma01 = rollmean(Subtracted, k = 3, fill = NA))  #http://uc-r.github.io/ts_moving_averages

baseline.peakDetection.smooth <- baseline.peakDetection.Narrow %>%
  select(RetentionTime,srate_ma01)

Q <- peaks(baseline.peakDetection.smooth,minPH=9000, minPW=0.3)

# save the plot
tiff(paste0(LibraryComparison.dir,sprintf("%s_PK.tiff",DataName)), units="in", width=5, height=5, res=300)
plot(baseline.peakDetection.smooth, type="l", xlab="m/z", ylab="Intensity /count", ylim = c(0, 500000))
points(Q,col="blue",cex=1.6)

dev.off()

# extract a list of peaks in Q and the equivalent MS data
List.of.peaks <- DataHeatMap %>%
  filter(RetentionTime %in% Q$x)

# calculate the relative MS peak intensity per Retention Time, normaile to 100 Max peak
List.of.peaks <- List.of.peaks %>%
  group_by(RetentionTime) %>%
  mutate(freq = (Intensity)/(max(Intensity))*100)

List.of.peaks <- List.of.peaks %>%
  select(RetentionTime,Mass,freq)

# ###############################################
# ##### Creation of a library list of peaks #####
# ###############################################
#                                  # this is temporary. Need a library to combine so take first sample for now
# Peak.library <- List.of.peaks    # this will be created by importing a library instead
# names(Peak.library)[3] <-"freqLib"
# 
# Peak.library <- as.data.frame(Peak.library)
# 
# # reshape the dataframe from wide to long
# Peak.library.wide <- reshape(Peak.library, idvar = "Mass", timevar = "RetentionTime", direction = "wide")
# 
# # Rename the peaks of the library
# for (i in 2:ncol(Peak.library.wide)) {
#   # start peak series names for library
#   j <- i-1
#   names(Peak.library.wide)[i] <- paste0("L",j)
# }
# 
# #Export to text file 
# write.csv(Peak.library.wide, file=paste0(Library.dir,"ImpurityLibrary.csv"), row.names = F)
#
#
#
# # create a list of library peaks
# List.Peak.Library <- group_split(Peak.library)
# 
# # remove Retention time to allow full join
# for (i in 1:length(List.Peak.Library)) {
#   #  test[[i]] <- test[[i]] %>%
#   #  select(Mass,freq)
#   List.Peak.Library[[i]]$RetentionTime <- NA
#   i <- i+1
# }

###############################################
#####      read library list of peaks     #####
###############################################

Peak.library.wide <- read.csv(file=paste0(Library.dir,"ImpurityLibrary.csv"), check.names=FALSE)

Peak.library.wide$Mass <- round(Peak.library.wide$Mass/Mass.precision)*Mass.precision

###############################################
#####          sample peaks               #####
###############################################

# peaks from the sample of interest
Peak.sample <- List.of.peaks

Peak.sample <- as.data.frame(Peak.sample)

# reshape the dataframe from wide to long
Peak.sample.wide <- reshape(Peak.sample, idvar = "Mass", timevar = "RetentionTime", direction = "wide")

# Rename the peaks of the library
for (i in 2:ncol(Peak.sample.wide)) {
  # start peak series names for library
  j <- i-1
  names(Peak.sample.wide)[i] <- paste0("S",j)
}

Peak.sample.wide$Mass <- round(Peak.sample.wide$Mass/Mass.precision)*Mass.precision

# # create a list of peaks for sample
# List.Peak.Sample <- group_split(List.of.peaks)
# 
# # remove Retention time to allow full join
# for (i in 1:length(List.Peak.Sample)) {
#   #  test[[i]] <- test[[i]] %>%
#   #  select(Mass,freq)
#   List.Peak.Sample[[i]]$RetentionTime <- NA
#   i <- i+1
# }
# 
# temp.merge.lists <- Map(c,List.Peak.Sample,List.Peak.Libray)
# 
# for (i in length(test)) {
#   for (j in length(test[[i]]$freq)){
# test$Ia2 <- test[[i]]$freq[j]*test[[i]]$freq[j]
# }
# }
# 
# temp2 <- test[[2]]$freq*test[[2]]$freq

###############################################
#####   join the library to the sample    #####
###############################################

# # # test demo library peaks example
# Peak.library.wide.test <- data.frame(Mass = c(100, 125.3, 125.7, 200.6),
#                   L1 = c(100, 50, 75, 0),
#                   L2 = c(100, 75, 50, 0),
#                   L3 = c(100, 0, 25, 25))
# 
# # test demo sample peaks example
# Peak.sample.wide.test <- data.frame(Mass = c(100, 125, 150, 200),
#                   S1 = c(100, 50, 75, 0),
#                   S2 = c(50, 100, 75, 0),
#                   S3 = c(25, 50, 0, 100),
#                   S4 = c(100, 0, 25, 10))
# 

Joined.sample.lib <- full_join(Peak.library.wide,Peak.sample.wide)

# replace NA wth 0 
Joined.sample.lib[is.na(Joined.sample.lib)] <- 0

# set value "l" for the number of columns associated to the number of peaks in the library
# set value "s" for the number of columns associated to the number of peaks in the sample

m <- ncol(Peak.library.wide) +1   # all cols included to set a floating limit in the joined dataframe  
l <- ncol(Peak.library.wide)-1   # to exclude the first col m/z and part of
s <- ncol(Peak.sample.wide)-1   # to exclude the first col m/z and part of 
t <- 2 # first column to store products of sample and library

rm(numerator)
numerator <- Joined.sample.lib %>%
  select(Mass)

# Each column in the library will be processed against the series of peak for the sample
# the column name will be take to the number of the library peak position and increment it by the sample peak number -1
# example: the first column for the library peak (i.e. 1) processed with the first sample peak (i.e. 1)
# will be stored in a col P1.0 but showing as P1 in the dataframe 

# This is for the nominator part of the equation, products of the two lists
# followed by their sum: SUM(A.i*B.i)
for(l in 2:ncol(Peak.library.wide)) {
  for(s in m:ncol(Joined.sample.lib)) {
    numerator[,t] <- Joined.sample.lib[l] * Joined.sample.lib[s]
    
    # give name to the columns associated with the product
    lib.peak <- names(Joined.sample.lib)[l]
    sam.peak <- names(Joined.sample.lib)[s]
    name.var <-paste0(lib.peak,".",sam.peak)
    names(numerator)[t] <- name.var
    
    t <- t+1
    }
}

SumNumerator <- as.data.frame(colSums(numerator[-1]))
names(SumNumerator)[1] <- c("Numerator")

# This is for the denominator part of the equation, the sum of the square of each peak,
# followed by their products, then their square root,
# followed by their sum: SUM(A.i*B.i)

# reset t for index in new dataframe 
t <- 2 # first column to square

# remove any existing database called Denominator.Library
rm(Denominator.Library)

Peak.library.wide[is.na(Peak.library.wide)] <- 0

# create one dataframe using the library.wide and Mass
Denominator.Library <- Peak.library.wide %>%
  select(Mass)

# Square all the peak intensity columns
for(l in 2:ncol(Peak.library.wide)) {
    Denominator.Library[,t] <- Peak.library.wide[l] * Peak.library.wide[l]
    
    # give name to the columns associated with the product
    lib.peak <- names(Peak.library.wide)[l]
    names(Denominator.Library)[t] <- lib.peak
    
    t <- t+1
  }

# fill first row with 0 prior summing all cols
Denominator.Library$Mass <- 0

SumDenominator.Library <- as.data.frame(t(colSums(Denominator.Library)))

# reset t for index in new dataframe 
t <- 2 # first column to square

# remove any existing database called Denominator.Sample
rm(Denominator.Sample)

Peak.sample.wide[is.na(Peak.sample.wide)] <- 0

# create one using the library.wide and Mass
Denominator.Sample <- Peak.sample.wide %>%
  select(Mass)

# Square all the peak intensity columns
for(l in 2:ncol(Peak.sample.wide)) {
  Denominator.Sample[,t] <- Peak.sample.wide[l] * Peak.sample.wide[l]
  
  # give name to the columns associated with the product
  sam.peak <- names(Peak.sample.wide)[l]
  names(Denominator.Sample)[t] <- sam.peak
  
  
  t <- t+1
}

# fill first row with 0 prior summing all cols
Denominator.Sample$Mass <- 0

SumDenominator.Sample <- as.data.frame(t(colSums(Denominator.Sample)))

# Join by mass the two denominators into one dataframe 
Join.Denominator <- full_join(SumDenominator.Library,SumDenominator.Sample)

##### calculate the products of the to denominators ##### 
# reset t for index in new dataframe 
t <- 2 # first column to square

# remove any existing database called Sqrt.Denominator
rm(Sq.Denominator)
# create one using the library.wide and Mass
Sq.Denominator <- Join.Denominator %>%
  select(Mass)

# Square all the peak intensity columns
for(l in 2:ncol(Join.Denominator)) {
  Sq.Denominator[,t] <- Join.Denominator[l] * Join.Denominator[l]
  
  # give name to the columns associated with the product
  name.peak <- names(Sq.Denominator)[l]
  names(Sq.Denominator)[t] <- name.peak
  
  t <- t+1
}

# remove any existing database called Denominator
rm(Sqrt.Denominator)
# create one using the library.wide and Mass
Sqrt.Denominator <- Sq.Denominator %>%
  select(Mass)

m <- ncol(Denominator.Library) + 1   # all cols included to set a floating limit in the joined dataframe  
t <- 2 # first column to store products of sample and library

for(l in 2:ncol(Denominator.Library)) {
  for(s in m:ncol(Join.Denominator)) {
    Sqrt.Denominator[,t] <- sqrt(Join.Denominator[l] * Join.Denominator[s])
    
    # give name to the columns associated with the product
    lib.peak <- names(Join.Denominator)[l]
    sam.peak <- names(Join.Denominator)[s]
    name.var <-paste0(lib.peak,".",sam.peak)
    names(Sqrt.Denominator)[t] <- name.var
    
    t <- t+1
  }
}

Denominator <- as.data.frame(t(Sqrt.Denominator[-1]))

#############################################
##### Combine numerator and denominator #####
#############################################

Ratio <- merge(SumNumerator, Denominator, by=0, all=TRUE)

# calculate angular vertor ratios
Ratio$ratio <- Ratio$Numerator/Ratio$V1

# split names for library and sample peaks attribution
Ratio$L.peak <- sub('\\..*', '', as.character(Ratio$Row.names))

Ratio$P.peak <- sub('.*\\.', '', as.character(Ratio$Row.names))

p <- ggplot(Ratio, aes(x=P.peak, y= L.peak, fill = ratio)) + 
  geom_tile()  +
  geom_text(aes(label = round(ratio, 3)), size=3)+
  scale_fill_gradient2(high="red",mid="white",low="blue", 
                       na.value="yellow", midpoint=angular.vertor)+
  xlab("Sample peaks") +
  ylab("Library peaks")

 show(p)

ggsave(
  sprintf("%s.tiff",DataName),
  plot = p,
  device = NULL,
  path = file.path(LibraryComparison.dir),
  scale = 1,
  width = 7.5,
  height = 5.0,
  units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = TRUE
)



# create two lists for
# peaks matching (or partially) matching the reference library
# peaks not matching the reference library

Sample.Peak.list <- Ratio %>%
  select(P.peak) %>%
  distinct()

Sample.Peak.list <- Sample.Peak.list$P.peak

rm(temp.peak.sample.not.list)
rm(temp.peak.sample.list)

for (series in Sample.Peak.list) {
  temp.peak.sample <- Ratio[Ratio$P.peak == series,]
  temp.peak.sample.max <- max(temp.peak.sample$ratio)
  temp.peak.sample$name <- DataName
  
  if (temp.peak.sample.max <= angular.vertor) {
    #Export to text file 
    # if the dataset does exist, append to it
    if (exists("temp.peak.sample.notlist")){
      temp.peak.sample.not.list <- rbind(temp.peak.sample.not.list, temp.peak.sample)
      write.csv(temp.peak.sample.not.list, file=paste0(LibraryNotFound.dir,sprintf("%s.csv",DataName)), row.names = F)
    }

    # if the dataset doesn't exist, create it
    if (!exists("temp.peak.sample.not.list")){
      write.csv(temp.peak.sample, file=paste0(LibraryNotFound.dir,sprintf("%s.csv",DataName)), row.names = F)
      temp.peak.sample.not.list <- temp.peak.sample
    }
  }
  
  else  {
    # if the merged dataset does exist, append to it
    if (exists("temp.peak.sample.list")){
      temp.peak.sample.list <- rbind(temp.peak.sample.list, temp.peak.sample)
      write.csv(temp.peak.sample.list, file=paste0(LibraryFound.dir,sprintf("%s.csv",DataName)), row.names = F)
    }

    # if the dataset doesn't exist, create it
    if (!exists("temp.peak.sample.list")){
      write.csv(temp.peak.sample, file=paste0(LibraryFound.dir,sprintf("%s.csv",DataName)), row.names = F)
      temp.peak.sample.list <- temp.peak.sample
    }
  }
}
  

}
