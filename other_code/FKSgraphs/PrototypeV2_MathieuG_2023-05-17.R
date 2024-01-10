# Importe le fichiers des sequences filtrees
SequencesFiltreesFKS1 <- readLines("FKS1_filtered.aln")

x <- data.frame(Identifiant = substr(grep(">", SequencesFiltreesFKS1, value = TRUE), start = 2, stop = 12))

lignesSpecies <- grep("\\[*\\]", SequencesFiltreesFKS1, value = TRUE)


EspecesNonFormate <- gsub("^.*\\[", "", lignesSpecies)
Especes <- gsub("\\]", "", EspecesNonFormate)
Especes

x <- cbind(x, Especes)

################################################################################

## Trouve hotspots1 et hotspots 2

ListeHotspots1 <- list()
ListeHotspots2 <- list()

for (i in grep(">", SequencesFiltreesFKS1))
{
    ListeHotspots1 <- append(ListeHotspots1, substr(SequencesFiltreesFKS1[i+17], start = 52, stop = 60))
    ListeHotspots2 <- append(ListeHotspots2, substr(SequencesFiltreesFKS1[i+38], start = 1, stop = 8))
}

Hotspot1 <- unlist(ListeHotspots1)
Hotspot2 <- unlist(ListeHotspots2)

x <- cbind(x, Hotspot1)
x <- cbind(x, Hotspot2)

################################################################################