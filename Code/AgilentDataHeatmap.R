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
DataName <- gsub(".*/", "", DataName)

#####################################################################
#####  re-organise converted data from  Proteo Wizard MsConvert #####
#####################################################################
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


# replace blanck with "NA"
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

names(DataHeatMap)[1] <- "RTime"
names(DataHeatMap)[2] <- "TIC"
names(DataHeatMap)[3] <- "Mass"
names(DataHeatMap)[4] <- "Intensity"

# convert all columns to numeric
DataHeatMap <- mutate_all(DataHeatMap, function(x) as.numeric(as.character(x)))

# convert the retention time to sec
DataHeatMap$RTime <- DataHeatMap$RTime *60

#plot figure including all data
# p <- ggplot(DataHeatMap, aes(RTime, Mass))+
#   geom_tile(aes(fill = Intensity))
# show(p)

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
  select(RTime, TIC) %>%
  distinct()

# plot figure including all data
# p <- ggplot(DataTicRentention, aes(RTime, TIC))+
#   geom_line()
# show(p)

########################################################
#######    baseline corrections using hyperSpec  #######
########################################################

# baseline corrections adapted from and using the library package hyperSpec
# for more information: https://cran.r-project.org/web/packages/hyperSpec/index.html
# 
# the data must be in hyperspectral format (this could be extended to creating matrix of data)

# convert data to hyperSpec format
file <- as.data.frame(DataTicRentention)
temp_hyper_0 <- new ("hyperSpec", wavelength = file [,1], spc = t (file [, -1]))

# apply the background
bl <- spc.rubberband (temp_hyper_0 [,, 185 ~ 3000], noise = 3000, df = 20)

# plot the background to check fit to data
# plot(bl)
# plot(temp_hyper_0,add=TRUE)

# subtract the background to the original data and plot to check
temp_hyper_1 <- temp_hyper_0-bl
# plot(temp_hyper_1, plot.args = list (ylim = c(0, 4500)))

# export the data in temp folder
write.txt.long(temp_hyper_1, file = paste0(TempData.dir, "temp.txt"))

# open the data from temp folder
SignalMinusBg <- read.table(file = paste0(TempData.dir, "temp.txt"), sep = "\t", header = TRUE)

names(SignalMinusBg)[1] <- "RTime"
names(SignalMinusBg)[2] <- "TIC"

# apply limits to RTime: peak integration boundaries
################################################################
#####             Internal standard Peak Area             ######
################################################################
PeakIS <- SignalMinusBg %>%
  filter( between(RTime, 400, 405) )

# plot figure including all data
# p <- ggplot(PeakIS, aes(RTime, TIC))+
#   geom_line()
# show(p)

idIS <- order(PeakIS$RTime)
I.S.PA <- sum(diff(PeakIS$RTime[idIS])*rollmean(PeakIS$TIC[idIS],2))

################################################################
#####                 Etizolam Peak Area                  ######
################################################################
PeakEtizolam <- SignalMinusBg %>%
  filter( between(RTime, 645, 670) )

# plot figure including all data
# p <- ggplot(PeakEtizolam, aes(RTime, TIC))+
#   geom_line()
# show(p)

idEti <- order(PeakEtizolam$RTime)
PA <- sum(diff(PeakEtizolam$RTime[idEti])*rollmean(PeakEtizolam$TIC[idEti],2))

Results <-data.frame(DataName,PA,I.S.PA)

write.table(Results, file=paste0(GcData.dir,sprintf("%s.csv",DataName)), sep = ",", row.names = FALSE)

}

