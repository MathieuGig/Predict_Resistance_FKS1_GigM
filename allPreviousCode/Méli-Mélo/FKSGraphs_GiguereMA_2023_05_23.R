library(ggplot2)
library(UpSetR)

# For FKS1 Hotspot1
ggplot(BeauTableauFKS1, aes(x = Hotspot1)) + geom_bar(fill = "orange2") + 
    theme(axis.text.x = element_text(angle = 90))

# For FKS1 Hotspot2
# For FKS2 Hotspot1
# For FKS2 Hotspot2

# fonction "unique" pour avoir par espèce unique...

# Use summary ?
# Look into factors to sort data ...

# voir UpSetR
# Mettre chaque phylum dans les bulles en dessous.
# Nombre de hotspot uniques
# combinaison

listInput <- list(EspecesUniques = unique(BeauTableauFKS1$Species), Isolats = unique())