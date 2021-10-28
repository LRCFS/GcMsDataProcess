#####################################################################
#####              Run Global Code in first instance            #####
#####################################################################
# Load Metadata files
filenameMetadata <- list.files(Metadata.dir, pattern=extensionXLSX, full.names=TRUE)

Metadata <- read_xlsx(filenameMetadata)

# create the date from DataName
Metadata$Date <- str_sub(Metadata$DataName,1,6)

# Load GC results
filenameGcData <- list.files(GcData.dir, pattern=extensionCSV, full.names=TRUE)

GcResults <- read.csv(filenameGcData, header = TRUE, encoding = "UTF-8")
names(GcResults)[1] <- "DataName"

# join GC results to metadata
CombinedResults <- full_join(Metadata,GcResults)

# peak area ratio
CombinedResults$ratio <- CombinedResults$PA/CombinedResults$I.S.PA

# calibration range
calibration <- c("25","50","75","100","150","200","250")
#convert to a dataframe
calibration <- as.data.frame(calibration)
calibration$calibration <- as.numeric(calibration$calibration)

# filter the calibration out of data
calibrationValues <- CombinedResults %>%
  filter(CombinedResults$Type=="Cal")

# bind calibration range and GC output together
CalibrationData <- cbind(calibration,calibrationValues)

# determine the quadratic fit 
model <- lm(CalibrationData$ratio ~ poly(CalibrationData$calibration, degree = 2, raw = T))

# extract the values from list
QuadraticValues <- model[[1]]

# convert to dataframe
QuadraticValues <- as.data.frame(QuadraticValues)

# solving for x the quadractic equation: y = c*x^2 + b*x + a
# x = (±(b^2+4c(y−a))^(1/2)−b)/2c
# The order from the list [1], [2] and [3] are for a, b and c
# only the positive part of the equation is considered
CombinedResults$ValuesPositive <- (((QuadraticValues$QuadraticValues[2])^2 + 4*QuadraticValues$QuadraticValues[3]*(CombinedResults$ratio-QuadraticValues$QuadraticValues[1]))^(1/2)-QuadraticValues$QuadraticValues[2])/(2*QuadraticValues$QuadraticValues[3])

# to run only for the first time an export needs to be created, to be commented afterward or all saved data will be overwritten.
# write.table(CombinedResults,file = paste0(Results.dir,"GCMSResults.csv"),  sep = ",", row.names = F)

# Load already processed data
filenameData <- list.files(Results.dir, pattern=extensionCSV, full.names=TRUE)
ProcessedData <- read.csv(filenameData, sep=",", header=TRUE)

# Combined the existing results to the new one
CombinedData <- rbind(ProcessedData,CombinedResults)

write.table(CombinedData,file = paste0(Results.dir,"GCMSResults.csv"),  sep = ",", row.names = F)

