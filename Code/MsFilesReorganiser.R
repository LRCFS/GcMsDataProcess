#
# Prior running, GC-MS data must be converted using Proteo Wizard/ MsConvert
# This code will convert the MS1 files to a readable format and save them in a folder called in the subsequent code 
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

# write the files in the converted folder
write.csv(DataHeatMap, file=paste0(GcDataConvertedRcode.dir,DataName,".csv"), row.names = F)
}
