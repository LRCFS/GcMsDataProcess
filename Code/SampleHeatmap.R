#####################################################################
#####              Run Global Code in first instance            #####
#####################################################################

# Load heatmap results
filenameHeatMap <- list.files(HeatMapData, pattern=extensionCSV, full.names=TRUE)

for (file in filenameHeatMap) {
HeatResults <- read.csv(file, header = TRUE, encoding = "UTF-8")
names(HeatResults)[1] <- "DataName"

# peak area ratio
HeatResults$ratio <- HeatResults$PA/HeatResults$I.S.PA

# filter the calibration out of data
calibrationValues <- HeatResults %>%
  filter(HeatResults$Type=="Cal")

# determine the quadratic fit 
model <- lm(calibrationValues$ratio ~ poly(calibrationValues$Concentration, degree = 2, raw = T))

# extract the values from list
QuadraticValues <- model[[1]]

# convert to dataframe
QuadraticValues <- as.data.frame(QuadraticValues)

# solving for x the quadractic equation: y = c*x^2 + b*x + a
# x = (±(b^2+4c(y−a))^(1/2)−b)/2c
# The order from the list [1], [2] and [3] are for a, b and c
# only the positive part of the equation is considered
HeatResults$ValuesPositive <- (((QuadraticValues$QuadraticValues[2])^2 + 4*QuadraticValues$QuadraticValues[3]*(HeatResults$ratio-QuadraticValues$QuadraticValues[1]))^(1/2)-QuadraticValues$QuadraticValues[2])/(2*QuadraticValues$QuadraticValues[3])

# Include the dilution and convertion to mg/ml: multiply by dilution dived by 1000
HeatResults$Conc <- HeatResults$ValuesPositive * HeatResults$Dilution /1000

# to run only for the first time an export needs to be created, to be commented afterward or all saved data will be overwritten.
# write.table(HeatResults,file = paste0(HeatResults.dir,"HeatResults.csv"),  sep = ",", row.names = F)

# Load already processed data
filenameData <- list.files(HeatResults.dir, pattern=extensionCSV, full.names=TRUE)
ProcessedData <- read.csv(filenameData, sep=",", header=TRUE)

# Combined the existing results to the new one
CombinedData <- rbind(ProcessedData,HeatResults)

write.table(CombinedData,file = paste0(HeatResults.dir,"HeatResults.csv"),  sep = ",", row.names = F)
}

#remove any duplicates and filter to sample only
CombinedData <- CombinedData %>%
  distinct()

#filter to sample only
CombinedDataNarrow <- CombinedData %>%
  filter(CombinedData$Type == "Sample")

CombinedDataNarrow$PositionY <- as.character(CombinedDataNarrow$PositionY)


# a plot with numbers
p <- ggplot(CombinedDataNarrow, aes(PositionX, PositionY)) +
  geom_tile(aes(fill = Conc)) + 
  geom_text(aes(label = round(Conc, 2))) +
  xlab("") + ylab("")+
  scale_fill_gradient(low = "white", high = "red")

# reserve y axis
# p + scale_y_reverse()

show(p)
