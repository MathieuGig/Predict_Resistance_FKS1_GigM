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

MutationCodon <- c()
MutationAA <- c()
MutationNt <- c()

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
    
    MutationCodon <- append(MutationCodon, DataAni$Nmut_codons[DataAni$X == i])
    MutationAA <- append(MutationAA, DataAni$Nham_aa[DataAni$X == i])
    MutationNt <- append(MutationNt, DataAni$Nham_nt[DataAni$X == i])
}

MutationCodon <- as.factor(MutationCodon)
MutationAA <- as.factor(MutationAA)
MutationNt <- as.factor(MutationNt)

## Make new data frame.
NewData <- data.frame(Anidulafungin = MedianSelectionCoeffAni,
                      Caspofungin = MedianSelectionCoeffCaspo,
                      Micafungin = MedianSelectionCoeffMica,
                      NMutationCodon = MutationCodon,
                      NMutationAA = MutationAA,
                      NMutationNt = MutationNt)

## Make pair plot.
#pairs(NewData[, 1:3], main = "Median Selection Coefficient",
#      col = c("pink", "red", "brown", "orange", "yellow", "green", "cyan", "blue", "purple")[NewData$NMutation])

library(ggplot2)
library(GGally)

## L'erreur était dû au fait que seulement 2 valeurs avaient un NMutation == 2.
# Solution 1: filtrer les données pour enlever les valeurs problématiques.
Graph <- ggpairs(NewData[c(-108, -109),], columns = 1:3, aes(color = NMutationCodon))
Graph

# Solution 2: Ne pas calculer la corrélation.
Graph <- ggpairs(NewData, columns = 1:3, upper = "blank", aes(color = NMutationCodon)) + ggtitle("Median Selection coefficient") 
Graph

# Solution 3: Colorier sur chaque graphique individuellement.
Graph <- ggpairs(NewData, columns = 1:3) 
wow <- Graph[2,1]
wow <- wow + aes(color = NMutationCodon) + ggtitle("Median Selection coefficient")
wow
