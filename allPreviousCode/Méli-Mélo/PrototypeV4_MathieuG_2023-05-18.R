###
### ExtractionHotspot(fichierAlignements)
###
##
## Trouve et extrait les hotspots 1 et 2 pour toutes les sequences du gène FKS 
## contenues dans un fichier *.aln. Affiche l'espèce et son identifiant.
##
##
## Arguments
##
## fichierAlignements: chaine de caractères
##
##
## Valeur
##
## Un data frame de 4 colonnes contenant les identifiants, les espèces,
## et les séquences de Hotspot 1 et 2 de ces espèces.
##
##
## Exemple
##
## ExtractionHotspot("FKS1_filtered.aln")
##

ExtractionHotspot <- function(fichierAlignements)
{
    # Validation du fichier input.
    stopifnot("Fichier invalide pour analyse. Besoin d'un fichier *.aln" = grepl(".aln", fichierAlignements))

    # Importe les fichiers des sequences filtrees.
    SequencesFiltrees <- readLines(fichierAlignements)

    # Formater les identifiants en utilisant des REGEX.
    IdentifiantFormate <- gsub("\\|.*", "", grep(">", SequencesFiltrees, value = TRUE))
    IdentifiantFormate <- gsub(" .*", "", IdentifiantFormate)
    IdentifiantFormate <- gsub(">", "", IdentifiantFormate)

    # Formater les espèces en utilisant des REGEX.
    EspecesNonFormate <- gsub("^.*\\[", "", grep("\\[*\\]", SequencesFiltrees, value = TRUE))
    EspecesFormate <- gsub("\\]", "", EspecesNonFormate)

    # Création du data frame.
    Tableau <- data.frame(Identifier = IdentifiantFormate, Species = EspecesFormate)

    # Création de liste des hotspots.
    ListeHotspots1 <- list()
    ListeHotspots2 <- list()

    # Trouve les séquences correspondant au Hotspots en utilisant leurs positions fixes.
    # Pour chaque ligne avec un chevron... donc pour chaque sequence...
    for (i in grep(">", SequencesFiltrees))
    {
        # Ajoute à la liste des hotspots les 9 acides aminées de ces positions fixes.
        ListeHotspots1 <- append(ListeHotspots1, substr(SequencesFiltrees[i+17], start = 52, stop = 60))
        ListeHotspots2 <- append(ListeHotspots2, substr(SequencesFiltrees[i+38], start = 1, stop = 8))
    }

    # Nomme les colonnes.
    Hotspot1 <- unlist(ListeHotspots1)
    Hotspot2 <- unlist(ListeHotspots2)

    # Ajoute les colonnes de Hotspots au data frame.
    Tableau <- cbind(Tableau, Hotspot1)
    Tableau <- cbind(Tableau, Hotspot2)

    # Retourne le tableau.
    Tableau
}

wow2 <- ExtractionHotspot("FKS2_filtered.aln")