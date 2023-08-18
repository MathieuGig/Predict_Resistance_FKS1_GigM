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
    library(cowplot)
    
    origin <- 0
    LabelCutOff <- 1
    
    Graph1 <- ggplot(BestData, aes(Anidulafungin, Caspofungin,
                                   color = MutatedAminoAcids,
                                   label = ifelse((Anidulafungin > LabelCutOff & Caspofungin < origin) | (Anidulafungin < origin & Caspofungin > LabelCutOff), Sequence, ""))) +
        geom_text(hjust = -0.1, vjust = -0.1) + geom_jitter()
    
    Graph2 <- ggplot(BestData, aes(Anidulafungin, Micafungin,
                                   color = MutatedAminoAcids,
                                   label = ifelse((Anidulafungin > LabelCutOff & Micafungin < origin) | (Anidulafungin < origin & Micafungin > LabelCutOff), Sequence, ""))) +
        geom_text(hjust = -0.1, vjust = -0.1) + geom_jitter()
        
    Graph3 <- ggplot(BestData, aes(Micafungin, Caspofungin,
                                   color = MutatedAminoAcids,
                                   label = ifelse((Micafungin > LabelCutOff & Caspofungin < origin) | (Micafungin < origin & Caspofungin > LabelCutOff), Sequence, ""))) +
        geom_text(hjust = -0.1, vjust = -0.1) + geom_jitter()

    title <- ggdraw() + 
        draw_label(
            paste(sample, "_", gene, "-", hotspot,sep = ""),
            fontface = 'bold',
            x = 0,
            hjust = 0
        ) +
        theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 97)
        )
    
    #plot_grid(title, Graph1, Graph3, Graph2)
    #ggsave(paste("pairPlot_", nameOfFile, ".png", sep = ""), width = 10, height = 10)
    write.csv(BestData, paste("Data_pairPlot_", nameOfFile, ".csv", sep = ""), row.names=FALSE)
}

makeComparisonPlot("BY4741", "FKS1", "HS1")
makeComparisonPlot("BY4741", "FKS1", "HS2")

makeComparisonPlot("R1158", "FKS1", "HS1")
makeComparisonPlot("R1158", "FKS1", "HS2")

makeComparisonPlot("R1158", "FKS2", "HS1")
makeComparisonPlot("R1158", "FKS2", "HS2")