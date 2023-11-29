makeBoxplot <- function(sample, drug)
{
    File <- paste(getwd(), "processed_data", paste(sample, "_FKS1-HS1_single_ortho_", drug, sep = ""), "selcoeff_all_libraries.csv", sep = "/")
    Readfile <- read.csv(File)
    
    BigData <- subset(Readfile, Readfile$seq_type == "ortho")
    
    NCodon <- c()
    Naa <- c()
    Nnt <- c()
    SelectionCoeff <- c()
    
    for (i in BigData$X)
    {
        NCodon <- append(NCodon, BigData$Nmut_codons[BigData$X == i])
        Naa <- append(Naa, BigData$Nham_aa[BigData$X == i])
        Nnt <- append(Nnt, BigData$Nham_nt[BigData$X == i]) 
        SelectionCoeff <- append(SelectionCoeff, median(c(BigData$selcoeff_1[BigData$X ==  i], BigData$selcoeff_2[BigData$X == i])))
    }
    
    NCodon <- as.factor(NCodon)
    Naa <- as.factor(Naa)
    Nnt <- as.factor(Nnt)
    
    #return(NCodon)
    
    NewData <- data.frame(NombreCodonMuté = NCodon, NombreAaMuté = Naa, NombreNtMuté = Nnt, SelectionCoeff)
    
    #return(NewData)
    
    library(ggplot2)
    
    # Graph x == NCodon
    GraphNCodon <- ggplot(data = NewData, aes(x = NombreCodonMuté, y = SelectionCoeff)) + geom_boxplot() + geom_jitter(position = position_jitter(0.2)) + labs(title = drug) + ylim(-1, 3)
    
    # Graph x == Naa
    GraphNaa <- ggplot(data = NewData, aes(x = NombreAaMuté, y = SelectionCoeff)) + geom_boxplot() + geom_jitter(position = position_jitter(0.2)) + labs(title = drug) + ylim(-1, 3)
    
    # Graph x == Nnt
    GraphNnt <- ggplot(data = NewData, aes(x = NombreNtMuté, y = SelectionCoeff)) + geom_boxplot() + geom_jitter(position = position_jitter(0.2)) + labs(title = drug) + ylim(-1, 3)
    
    graphList <- list(GraphNCodon, GraphNaa, GraphNnt)
    graphList
    
}

makeBoxplot("BY4741", "anidulafungin")
#makeBoxplot("BY4741", "caspofungin")



# Readfile <- read.csv("C:/Users/User/Documents/stageEte2023/processed_data/BY4741_FKS1-HS1_single_ortho_anidulafungin/selcoeff_all_libraries.csv")


