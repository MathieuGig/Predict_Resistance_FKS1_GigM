
## Importe les fichiers de hotspots uniques dans des variables.
FKS1_Hotspots1 <- (readLines("FKS1_HS1_orthologs_unique.fa"))[-1]

FKS1_Hotspots2 <- (readLines("FKS1_HS2_orthologs_unique.fa"))[-1]

SequencesFiltrees <- readLines("FKS1_filtered.aln")

## La fonction qui accomplira le travail
Trouve_Identifiant_Hotspot <- function(sequences, hotspots)
{
    # Création de la liste des résultat de recherche des hotspots dans les séquences
    listeResultatGrep <- list()

    # Pour chaque séquence de hotspot, trouve les lignes du fichier où les sequences matchs
    for (i in hotspots)
    {
        listeResultatGrep <- append(listeResultatGrep, list(grep(i, sequences)))
    }
    names(listeResultatGrep) <- hotspots
    listeResultatGrep
}


# La commande
Trouve_Identifiant_Hotspot(SequencesFiltrees, FKS1_Hotspots1)