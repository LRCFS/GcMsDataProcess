
CardData <- read.csv("CardData.csv", sep=",", header=TRUE)
names(CardData)[1] <-"SampleName"

CardData <- CardData[!is.na(CardData$Dilution), ]


TempDataCard <- full_join(CardData,CombinedData)
