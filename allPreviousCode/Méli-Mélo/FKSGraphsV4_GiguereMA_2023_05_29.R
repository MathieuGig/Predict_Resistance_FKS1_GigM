## Import libraries
library(ggplot2)
library(UpSetR)

################################################################################

## Bar Plots
# For FKS1 Hotspot1
ggplot(TableauFKS1, aes(x = Hotspot1)) + geom_bar(fill = "orange2") + 
    theme(axis.text.x = element_text(angle = 90)) + ggtitle("Count FKS1 Hotspot 1 sequences")

# For FKS1 Hotspot2
ggplot(TableauFKS1, aes(x = Hotspot2)) + geom_bar(fill = "orange2") +
    theme(axis.text.x = element_text(angle = 90)) + ggtitle("Count FKS1 Hotspot 2 sequences")

# For FKS2 Hotspot1
ggplot(TableauFKS2, aes(x = Hotspot1)) + geom_bar(fill = "green3") +
    theme(axis.text.x = element_text(angle = 90)) + ggtitle("Count FKS2 Hotspot 1 sequences")

# For FKS2 Hotspot2
ggplot(TableauFKS2, aes(x = Hotspot2)) + geom_bar(fill = "green3") +
    theme(axis.text.x = element_text(angle = 90)) + ggtitle("Count FKS2 Hotspot 2 sequences")

################################################################################

## Premier Upset graph. Observe FKS1 Hotspot1 per Phylum.
listInput <- list(Chlorophyta = TableauFKS1$Hotspot1[TableauFKS1$Phylum == "Chlorophyta"],
                  Streptophyta = TableauFKS1$Hotspot1[TableauFKS1$Phylum == "Streptophyta"],
                  Mucoromycota = TableauFKS1$Hotspot1[TableauFKS1$Phylum == "Mucoromycota"],
                  Ascomycota = TableauFKS1$Hotspot1[TableauFKS1$Phylum == "Ascomycota"],
                  Basidiomycota = TableauFKS1$Hotspot1[TableauFKS1$Phylum == "Basidiomycota"])

upset(fromList(listInput), order.by = "freq",
      mainbar.y.label = "Hotspot1 sequence intersection",
      sets.x.label = "Unique hotspot1 sequence per phylum")

################################################################################

## Graph for FKS1 Hotspot1 of Ascomycota species having 2 or more hotspot1 sequences.
AscomycotaDoubles <- unique(TableauFKS1$Species[duplicated(TableauFKS1$Species[TableauFKS1$Phylum == "Ascomycota"])])

Input <- list()
for (i in AscomycotaDoubles)
{
    Input <- append(Input, list(TableauFKS1$Hotspot1[TableauFKS1$Species == as.character(i)]))
}
names(Input) <- AscomycotaDoubles
Input
upset(fromList(Input), order.by = "freq",
      mainbar.y.label = "Hotspot1 sequence intersection",
      sets.x.label = "Unique hotspot1 sequence per Ascomycota specie")

################################################################################

## Graph for FKS1 Hotspot 1 of all Ascomycota species.
Ascomycota <- unique(TableauFKS1$Species[TableauFKS1$Phylum == "Ascomycota"])

Input <- list()
for (i in Ascomycota)
{
    Input <- append(Input, list(TableauFKS1$Hotspot1[TableauFKS1$Species == as.character(i)]))
}
names(Input) <- Ascomycota
upset(fromList(Input), order.by = "freq",
      mainbar.y.label = "Hotspot1 sequence intersection",
      sets.x.label = "Unique hotspot1 sequence per Ascomycota specie")

################################################################################

## Graph FKS1 Hotspot1 for all species.
Input <- list()
for (i in unique(TableauFKS1$Species))
{
    Input <- append(Input, list(TableauFKS1$Hotspot1[TableauFKS1$Species == as.character(i)]))   
}
names(Input) <- unique(TableauFKS1$Species)

upset(fromList(Input), nsets = 10 , order.by = "freq",
      mainbar.y.label = "FKS1 Hotspot1 sequence intersection",
      sets.x.label = "hotspot1 sequence per specie")

################################################################################

## Graph FKS1 Hotspot2 for all species.
Input <- list()
for (i in unique(TableauFKS1$Species))
{
    Input <- append(Input, list(TableauFKS1$Hotspot2[TableauFKS1$Species == as.character(i)]))   
}
names(Input) <- unique(TableauFKS1$Species)

upset(fromList(Input), nsets = 10 , order.by = "freq",
      mainbar.y.label = "FKS1 Hotspot2 sequence intersection",
      sets.x.label = "hotspot2 sequence per specie")

################################################################################

## Graph FKS2 Hotspot1 for all species.
Input <- list()
for (i in unique(TableauFKS2$Species))
{
    Input <- append(Input, list(TableauFKS2$Hotspot1[TableauFKS2$Species == as.character(i)]))   
}
names(Input) <- unique(TableauFKS2$Species)

upset(fromList(Input), nsets = 10 , order.by = "freq",
      mainbar.y.label = "FKS2 Hotspot1 sequence intersection",
      sets.x.label = "hotspot1 sequence per specie")

################################################################################

## Graph FKS2 Hotspot2 for all species.
Input <- list()
for (i in unique(TableauFKS2$Species))
{
    Input <- append(Input, list(TableauFKS2$Hotspot2[TableauFKS2$Species == as.character(i)]))   
}
names(Input) <- unique(TableauFKS2$Species)

upset(fromList(Input), nsets = 10 , order.by = "freq",
      mainbar.y.label = "FKS2 Hotspot2 sequence intersection",
      sets.x.label = "hotspot2 sequence per specie")

################################################################################