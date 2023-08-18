makeComparisonPlot <- function(sample, gene, hotspot)
{
    drug <- c("anidulafungin", "caspofungin", "micafungin")
    
    nameOfFile <- paste(sample, "_", gene, "-", hotspot, "_single_ortho_", sep = "")
    
    File <- paste(getwd(), "processed_data", paste(nameOfFile, drug, sep = ""),
                  "selcoeff_all_libraries.csv", sep = "/")
    
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
    
    MutationAA <- c()
    AASequence <- c()
    
    for (i in DataAni$X)
    {
        MedianSelectionCoeffAni <- append(MedianSelectionCoeffAni, DataAni$median_s[DataAni$X == i])
        
        MedianSelectionCoeffCaspo <- append(MedianSelectionCoeffCaspo, DataCaspo$median_s[DataCaspo$X == i])
        
        MedianSelectionCoeffMica <- append(MedianSelectionCoeffMica, DataMica$median_s[DataMica$X == i])
        
        AASequence <- append(AASequence, DataAni$aa_seq[DataAni$X == i])
        MutationAA <- append(MutationAA, DataAni$Nham_aa[DataAni$X == i])
    }
    
    MutationAA <- as.factor(MutationAA)
    
    ## Make new data frame.
    NewData <- data.frame(Anidulafungin = MedianSelectionCoeffAni,
                          Caspofungin = MedianSelectionCoeffCaspo,
                          Micafungin = MedianSelectionCoeffMica,
                          MutatedAminoAcids = MutationAA,
                          Sequence = AASequence)
    
    BestData <- unique(NewData)
    
    library(ggplot2)
    library(GGally)
    
    Graph1 <- ggplot(BestData, aes(Anidulafungin, Caspofungin,
                                   color = MutatedAminoAcids,
                                   label = ifelse((Anidulafungin > 1 & Caspofungin < 1) | (Anidulafungin < 1 & Caspofungin > 1), Sequence, ""))) +
        geom_text(hjust = -0.1, vjust = 0) + geom_jitter() +
        ggtitle(paste(nameOfFile, "Median Selection coefficient", sep = " ")) + xlim(-0.5, 2.5) + ylim(-0.5, 2.5)
    
    Graph2 <- ggplot(BestData, aes(Anidulafungin, Micafungin,
                                   color = MutatedAminoAcids,
                                   label = ifelse((Anidulafungin > 1 & Micafungin < 1) | (Anidulafungin < 1 & Micafungin > 1), Sequence, ""))) +
        geom_text(hjust = -0.1, vjust = 0) + geom_jitter() +
        ggtitle(paste(nameOfFile, "Median Selection coefficient", sep = " ")) + xlim(-0.5, 2.5) + ylim(-0.5, 2.5)
    
    Graph3 <- ggplot(BestData, aes(Caspofungin, Micafungin,
                                   color = MutatedAminoAcids,
                                   label = ifelse((Caspofungin > 1 & Micafungin < 1) | (Caspofungin < 1 & Micafungin > 1), Sequence, ""))) +
        geom_text(hjust = -0.1, vjust = 0) + geom_jitter() +
        ggtitle(paste(nameOfFile, "Median_selection_coefficient", sep = " ")) + xlim(-0.5, 2.5) + ylim(-0.5, 2.5)
    
    GraphList <- list(Graph1, Graph2, Graph3)
    GraphList
}

makeComparisonPlot("BY4741", "FKS1", "HS1")
#makeComparisonPlot("BY4741", "FKS1", "HS2")

#makeComparisonPlot("R1158", "FKS1", "HS1")
#makeComparisonPlot("R1158", "FKS1", "HS2")

#makeComparisonPlot("R1158", "FKS2", "HS1")
#makeComparisonPlot("R1158", "FKS2", "HS2")

# Parametres axes.
# cow plot
