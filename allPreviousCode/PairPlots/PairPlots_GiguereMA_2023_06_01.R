### Algorithm
# Import Data.
# Compute median of sel_coeffs.
# Make new data frame.
# Make pair plots using ggplot2.

################################################################################

# Import data from files.
drug <- c("anidulafungin", "caspofungin", "micafungin")

File <- paste(getwd(), "processed_data", paste("BY4741_FKS1-HS1_single_ortho_", drug, sep = ""), "selcoeff_all_libraries.csv", sep = "/")

Readfile <- list()

# Readfile <- list(list(read.csv(File[1])), list(read.csv(File[2])), list(File[3]))

for (i in File)
{
    # What and how do I want to append/rbind/cbind ?
    Readfile <- append(Readfile, list(read.csv(i)))
}

# wow this works
DataAni <- subset(Readfile[[1]], Readfile[[1]]$seq_type == "ortho")
DataCaspo <- subset(Readfile[[2]], Readfile[[2]]$seq_type == "ortho")
DataMica <- subset(Readfile[[3]], Readfile[[3]]$seq_type == "ortho")


# BigDataAni <- BigDataAni[, c("X", "seq_type", "nt_seq", "aa_seq", "Nham_nt", "Nham_aa", "Nmut_codons", "selcoeff_1", "selcoeff_2")]

MedianSelectionCoeffAni <- c()
MedianSelectionCoeffCaspo <- c()
MedianSelectionCoeffMica <- c()

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
}

NewData <- data.frame(Anidulafungin = MedianSelectionCoeffAni,
                      Caspofungin = MedianSelectionCoeffCaspo,
                      Micafungin = MedianSelectionCoeffMica)

# Make pair plot
pairs(NewData)

library(ggplot2)
library(GGally)

ggpairs(NewData) + ggtitle("Median Selection coefficient")
