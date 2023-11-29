### Algorithm
# Import Data.
# Compute median of sel_coeffs.
# Make new data frame.
# Make pair plots using ggplot2.

################################################################################

## Import data from files.
drug <- c("anidulafungin", "caspofungin", "micafungin")

File <- paste(getwd(), "processed_data", paste("BY4741_FKS1-HS1_single_ortho_", drug, sep = ""), "selcoeff_all_libraries.csv", sep = "/")

Readfile <- list()

for (i in File)
{
    Readfile <- append(Readfile, list(read.csv(i)))
}

DataAni <- subset(Readfile[[1]], Readfile[[1]]$seq_type == "ortho")
DataCaspo <- subset(Readfile[[2]], Readfile[[2]]$seq_type == "ortho")
DataMica <- subset(Readfile[[3]], Readfile[[3]]$seq_type == "ortho")

## Compute median sel_coeffs.
MedianSelectionCoeffAni <- c()
MedianSelectionCoeffCaspo <- c()
MedianSelectionCoeffMica <- c()

Mutation <- c()

for (i in DataAni$X)
{
    MedianSelectionCoeffAni <- append(MedianSelectionCoeffAni,
                                      median(c(DataAni$selcoeff_1[DataAni$X ==  i],
                                               DataAni$selcoeff_2[DataAni$X == i])))
    
    MedianSelectionCoeffCaspo <- append(MedianSelectionCoeffCaspo,
                                        median(c(DataCaspo$selcoeff_1[DataCaspo$X ==  i],
                                                 DataCaspo$selcoeff_2[DataCaspo$X == i])))
    
    MedianSelectionCoeffMica <- append(MedianSelectionCoeffMica,
                                        median(c(DataMica$selcoeff_1[DataMica$X ==  i],
                                                 DataMica$selcoeff_2[DataMica$X == i])))
    
    Mutation <- append(Mutation, DataAni$Nmut_codons[DataAni$X == i])
}

Mutation <- as.factor(Mutation)

## Make new data frame.
NewData <- data.frame(Anidulafungin = MedianSelectionCoeffAni,
                      Caspofungin = MedianSelectionCoeffCaspo,
                      Micafungin = MedianSelectionCoeffMica,
                      NMutation = Mutation)

## Make pair plot.
pairs(NewData)

library(ggplot2)
library(GGally)

ggpairs(NewData, columns = 1:3) + ggtitle("Median Selection coefficient")

# mapping = ggplot2::aes(colour = NMutation)
