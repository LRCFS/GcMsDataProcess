#
# Need to run GlobalCode first. 
# 

# Load Data Set GCMS Results
filenameData <- list.files(Results.dir, pattern=extensionCSV, full.names=TRUE)
ProcessedData <- read.csv(filenameData, sep=",", header=TRUE)
# filter the calibration out of data
calibrationValues <- ProcessedData %>%
  filter(ProcessedData$Type=="Cal") %>%
  select(ratio, ValuesPositive, Date)

calibrationValues$Date <- as.character(calibrationValues$Date)

# Plot the results

p<- ggplot(calibrationValues, aes(x=ValuesPositive, y=ratio)) + 
  geom_point(aes(colour = factor(Date)), size = 2) +
  # geom_line(data=calibrationValues, aes(x=ValuesPositive, y=ratio), colour="green") +
  theme_bw() +
  ylab("Etizolam - Internal Standard peak area ratio / arb. unit") +
  xlab(bquote("Concentration Etizolam /" ~mu~g%.%mL^{-1})) +
  theme(text = element_text(size = 12))


show(p)

# output figure can be saved if required.