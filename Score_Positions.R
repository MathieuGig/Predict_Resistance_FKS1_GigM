###
### makeBoxplot(sample, gene, hotspot, drug)
###
##
## Generates a boxplot graph using ggplot and the DMS data.
##
##
## Arguments
##
## sample: character string, the name of the analysed sample.
## gene: character string, the name of the analysed gene.
## hotspot: character string, HS1 or HS2.
## drug: character string, the name of the antifungal drug.
##
##
## Value
##
## A boxplot graph representing the median selection coefficient as a function
## of the number of mutated amino acids.
##
##
## Example
##
## makeBoxplot("BY4741", "FKS1", "HS1" "anidulafungin")
##

makeBoxplot <- function(sample, gene, hotspot, drug)
{
    # Validation of arguments
    stopifnot("Invalid sample" = sample %in% c("BY4741", "R1158"))
    stopifnot("Invalid gene" = gene %in% c("FKS1", "FKS2"))
    stopifnot("Invalid hotspot" = hotspot %in% c("HS1", "HS2"))
    stopifnot("Invalid drug" = drug %in% c("anidulafungin", "caspofungin", "micafungin"))
    
    # Import data from file
    File <- paste(getwd(), "processed_data", paste(sample, "_", gene, "-", hotspot, "_single_ortho_", drug, sep = ""), "selcoeff_all_libraries.csv", sep = "/")
    Readfile <- read.csv(File)
    BigData <- subset(Readfile, Readfile$seq_type == "ortho")
    
    AApos <- c()
    SelectionCoeff <- c()
    
    # Make simpler data frame from the imported data.
    for (i in BigData$X)
    {
        AApos <- append(AApos, BigData$aa_pos[BigData$X == i])
        SelectionCoeff <- append(SelectionCoeff, BigData$median_s[BigData$X == i])
    }
    
    AApos <- as.factor(AApos)
    
    Data <- data.frame(Mutated_Position = AApos, SelectionCoeff)
    
    # Load library to make graphs.
    library(ggplot2)
    
    # Custom graph title
    graphTitle <- paste(sample, "_", gene, "-", hotspot, "_single_ortho_", drug, sep = "")
    
    # Make Graph where x axis = Naa.
    Graph_Mut_Pos <- ggplot(data = Data, aes(x = Mutated_Position, y = SelectionCoeff)) +
        geom_boxplot(fill = "springgreen")+ geom_jitter(position = position_jitter(0.2)) +
        labs(title = graphTitle)
    
    # Return graph
    Graph_Mut_Pos
}

### Testing
makeBoxplot("BY4741", "FKS1", "HS1", "anidulafungin")
makeBoxplot("BY4741", "FKS1", "HS1", "caspofungin")
makeBoxplot("BY4741", "FKS1", "HS1", "micafungin")

makeBoxplot("BY4741", "FKS1", "HS2", "anidulafungin")
makeBoxplot("BY4741", "FKS1", "HS2", "caspofungin")
makeBoxplot("BY4741", "FKS1", "HS2", "micafungin")


makeBoxplot("R1158", "FKS1", "HS1", "anidulafungin")
makeBoxplot("R1158", "FKS1", "HS1", "caspofungin")
makeBoxplot("R1158", "FKS1", "HS1", "micafungin")

makeBoxplot("R1158", "FKS1", "HS2", "anidulafungin")
makeBoxplot("R1158", "FKS1", "HS2", "caspofungin")
makeBoxplot("R1158", "FKS1", "HS2", "micafungin")


makeBoxplot("R1158", "FKS2", "HS1", "anidulafungin")
makeBoxplot("R1158", "FKS2", "HS1", "caspofungin")
makeBoxplot("R1158", "FKS2", "HS1", "micafungin")

makeBoxplot("R1158", "FKS2", "HS2", "anidulafungin")
makeBoxplot("R1158", "FKS2", "HS2", "caspofungin")
makeBoxplot("R1158", "FKS2", "HS2", "micafungin")