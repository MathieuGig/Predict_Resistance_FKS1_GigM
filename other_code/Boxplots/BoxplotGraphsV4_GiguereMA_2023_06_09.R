###
### makeBoxplot(sample, drug)
###
##
## Generates 3 boxplot graphs using ggplot and the DMS data.
##
##
## Arguments
##
## sample: character string, the name of the analysed sample.
## drug: character string, the name of the antifungal drug.
##
##
## Value
##
## 3 boxplot graphs representing the median selection coefficient as a function
## of the number of mutated codons, the number of mutated amino acids,
## the number of mutated nucleotides respectively.
##
##
## Example
##
## makeBoxplot("BY4741", "anidulafungin")
##

makeBoxplot <- function(sample, drug)
{
    # Validation of arguments
    stopifnot("Invalid sample" = sample %in% c("BY4741", "R1158"))
    stopifnot("Invalid drug" = drug %in% c("anidulafungin", "caspofungin", "micafungin"))

    # Import data from file
    File <- paste(getwd(), "processed_data", paste(sample, "_FKS1-HS1_single_ortho_", drug, sep = ""), "selcoeff_all_libraries.csv", sep = "/")
    Readfile <- read.csv(File)
    BigData <- subset(Readfile, Readfile$seq_type == "ortho")

    NCodon <- c()
    Naa <- c()
    Nnt <- c()
    SelectionCoeff <- c()

    # Make simpler data frame from the imported data.
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

    NewData <- data.frame(MutatedCodons = NCodon, MutatedAA = Naa, MutatedNucleotides = Nnt, SelectionCoeff)
    
    # Load library to make graphs.
    library(ggplot2)
    
    # Custom graph title
    graphTitle <- paste(sample, "_FKS1-HS1_single_ortho_", drug, sep = "")
    
    # Make Graph where x axis = NCodon.
    GraphNCodon <- ggplot(data = NewData, aes(x = MutatedCodons, y = SelectionCoeff)) +
        geom_boxplot(fill = "cyan") + geom_jitter(position = position_jitter(0.2)) +
        labs(title = graphTitle) + ylim(-1, 3)
    
    # Make Graph where x axis = Naa.
    GraphNaa <- ggplot(data = NewData, aes(x = MutatedAA, y = SelectionCoeff)) +
        geom_boxplot(fill = "springgreen") + geom_jitter(position = position_jitter(0.2)) +
        labs(title = graphTitle) + ylim(-1, 3)
    
    # Make Graph where x axis = Nnt.
    GraphNnt <- ggplot(data = NewData, aes(x = MutatedNucleotides, y = SelectionCoeff)) +
        geom_boxplot(fill = "moccasin") + geom_jitter(position = position_jitter(0.2)) +
        labs(title = graphTitle) + ylim(-1, 3)

    # Return all 3 graphs.
    graphList <- list(GraphNCodon, GraphNaa, GraphNnt)
    graphList
    
}

### Testing
makeBoxplot("BY4741", "anidulafungin")
#makeBoxplot("BY4741", "caspofungin")
#makeBoxplot("BY4741", "micafungin")

#makeBoxplot("R1158", "anidulafungin")
#makeBoxplot("R1158", "caspofungin")
#makeBoxplot("R1158", "micafungin")
