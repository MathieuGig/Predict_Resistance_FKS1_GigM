# Import libraries
library(ggplot2)
library(UpSetR)

# For FKS1 Hotspot1
ggplot(BeauTableauFKS1, aes(x = Hotspot1)) + geom_bar(fill = "orange2") + 
    theme(axis.text.x = element_text(angle = 90))
# For FKS1 Hotspot2
# For FKS2 Hotspot1
# For FKS2 Hotspot2

## Premier Upset graph. Observe Hotspot1 per Phylum
listInput <- list(Chlorophyta = BeauTableauFKS1$Hotspot1[BeauTableauFKS1$Phylum == "Chlorophyta"],
                  Streptophyta = BeauTableauFKS1$Hotspot1[BeauTableauFKS1$Phylum == "Streptophyta"],
                  Mucoromycota = BeauTableauFKS1$Hotspot1[BeauTableauFKS1$Phylum == "Mucoromycota"],
                  Ascomycota = BeauTableauFKS1$Hotspot1[BeauTableauFKS1$Phylum == "Ascomycota"],
                  Basidiomycota = BeauTableauFKS1$Hotspot1[BeauTableauFKS1$Phylum == "Basidiomycota"])

upset(fromList(listInput), order.by = "freq",
      mainbar.y.label = "Hotspot1 sequence intersection",
      sets.x.label = "Unique hotspot1 sequence per phylum")


## Graph for Hotspot1 of Ascomycota specie having 2 or more hotspot1 sequences
AscomycotaDoubles <- unique(BeauTableauFKS1$Species[duplicated(BeauTableauFKS1$Species[BeauTableauFKS1$Phylum == "Ascomycota"])])

Input <- list()
for (i in AscomycotaDoubles)
{
    Input <- append(Input, list(BeauTableauFKS1$Hotspot1[BeauTableauFKS1$Species == as.character(i)]))
}
names(Input) <- AscomycotaDoubles
Input
upset(fromList(Input), order.by = "freq",
      mainbar.y.label = "Hotspot1 sequence intersection",
      sets.x.label = "Unique hotspot1 sequence per Ascomycota specie")



################################################################################
# Attempting to optimize
makeUpsetGraph <- function(ObserveY, PerX)
{
    Input <- list()
    for (i in PerX)
    {
        Input <- append(Input, list(ObserveY[PerX == as.character(i)]))
    }
    names(Input) <- PerX
    graph <- upset(fromList(Input), order.by = "freq")
    graph
}
makeUpsetGraph(BeauTableauFKS1$Hotspot1, unique(BeauTableauFKS1$Species[BeauTableauFKS1$Phylum == "Ascomycota"]))


################################################################################

## Graph pour les ascomycota sans filtre
Ascomycota <- unique(BeauTableauFKS1$Species[BeauTableauFKS1$Phylum == "Ascomycota"])

Input <- list()
for (i in Ascomycota)
{
    Input <- append(Input, list(BeauTableauFKS1$Hotspot1[BeauTableauFKS1$Species == as.character(i)]))
}
names(Input) <- Ascomycota
upset(fromList(Input), order.by = "freq",
      mainbar.y.label = "Hotspot1 sequence intersection",
      sets.x.label = "Unique hotspot1 sequence per Ascomycota specie")

################################################################################

## Graph FKS1 Hotspot1 for all species
Input <- list()
for (i in unique(BeauTableauFKS1$Species))
{
    Input <- append(Input, list(BeauTableauFKS1$Hotspot1[BeauTableauFKS1$Species == as.character(i)]))   
}
names(Input) <- unique(BeauTableauFKS1$Species)

upset(fromList(Input), nsets = 10 , order.by = "freq",
      mainbar.y.label = "Hotspot1 sequence intersection",
      sets.x.label = "hotspot1 sequence per specie")
