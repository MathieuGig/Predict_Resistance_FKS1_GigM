## Importe les fichiers des sequences filtrees
SequencesFiltreesFKS1 <- readLines("FKS1_filtered.aln")
SequencesFiltreesFKS2 <- readLines("FKS2_filtered.aln")

## Formater les identifiants en utilisant des REGEX
# Pour FKS1
FKS1_IdentifiantFormate <- gsub("\\|.*", "", grep(">", SequencesFiltreesFKS1, value = TRUE))
FKS1_IdentifiantFormate <- gsub(" .*", "", FKS1_IdentifiantFormate)
FKS1_IdentifiantFormate <- gsub(">", "", FKS1_IdentifiantFormate)

# Pour FKS2
FKS2_IdentifiantFormate <- gsub("\\|.*", "", grep(">", SequencesFiltreesFKS2, value = TRUE))
FKS2_IdentifiantFormate <- gsub(" .*", "", FKS2_IdentifiantFormate)
FKS2_IdentifiantFormate <- gsub(">", "", FKS2_IdentifiantFormate)

## Formater les espèces en utilisant des REGEX
# Pour FKS1
FKS1_EspecesNonFormate <- gsub("^.*\\[", "", grep("\\[*\\]", SequencesFiltreesFKS1, value = TRUE))
FKS1_EspecesFormate <- gsub("\\]", "", FKS1_EspecesNonFormate)

# Pour FKS2
FKS2_EspecesNonFormate <- gsub("^.*\\[", "", grep("\\[*\\]", SequencesFiltreesFKS2, value = TRUE))
FKS2_EspecesFormate <- gsub("\\]", "", FKS2_EspecesNonFormate)

## Création du data frame
xFKS1 <- data.frame(Identifier = FKS1_IdentifiantFormate, Species = FKS1_EspecesFormate)

## Trouve hotspots1 et hotspots 2

# Création de liste
FKS1_ListeHotspots1 <- list()
FKS1_ListeHotspots2 <- list()

# Trouve les séquences correspondant au Hotspots en utilisant leurs positions fixes.
# Pour chaque ligne avec un chevron... donc pour chaque sequence de FKS1...
for (i in grep(">", SequencesFiltreesFKS1))
{
    # Ajoute à la liste des hotspots les 9 acides aminées de ces positions fixes.
    FKS1_ListeHotspots1 <- append(FKS1_ListeHotspots1, substr(SequencesFiltreesFKS1[i+17], start = 52, stop = 60))
    FKS1_ListeHotspots2 <- append(FKS1_ListeHotspots2, substr(SequencesFiltreesFKS1[i+38], start = 1, stop = 8))
}

# Nomme les colonnes
FKS1_Hotspot1 <- unlist(FKS1_ListeHotspots1)
FKS1_Hotspot2 <- unlist(FKS1_ListeHotspots2)

## Ajoute au data frame
xFKS1 <- cbind(xFKS1, FKS1_Hotspot1)
xFKS1 <- cbind(xFKS1, FKS1_Hotspot2)