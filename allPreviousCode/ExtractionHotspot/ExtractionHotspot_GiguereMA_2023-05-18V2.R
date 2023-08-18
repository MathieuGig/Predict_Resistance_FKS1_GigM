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

    # Importer les fichiers des sequences filtrées.
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
        # Si FKS1...
        if (grepl("FKS1", fichierAlignements))
        {
            # Ajoute à la liste des hotspots les 8 ou 9 acides aminées de ces positions fixes.
            ListeHotspots1 <- append(ListeHotspots1, substr(SequencesFiltrees[i+17], start = 52, stop = 60))
            ListeHotspots2 <- append(ListeHotspots2, substr(SequencesFiltrees[i+38], start = 1, stop = 8))   
        }
        # Si FKS2...
        else
        {
            # Ajoute à la liste des hotspots les 8 acides aminées de ces positions fixes.
            # Le Hotspot1 de FKS2 est une ligne plus loin et 1 aa plus petit que celui de FKS1
            # Le Hotspot2 de FKS2 est plus loin sur sa ligne que celui de FKS1
            ListeHotspots1 <- append(ListeHotspots1, substr(SequencesFiltrees[i+18], start = 52, stop = 59))
            ListeHotspots2 <- append(ListeHotspots2, substr(SequencesFiltrees[i+38], start = 33, stop = 40))
        }
    }

    # Nommer les colonnes.
    Hotspot1 <- unlist(ListeHotspots1)
    Hotspot2 <- unlist(ListeHotspots2)

    # Ajouter les colonnes de Hotspots au data frame.
    Tableau <- cbind(Tableau, Hotspot1)
    Tableau <- cbind(Tableau, Hotspot2)
    
    # Ajouter une colonne qui indique la présence de gaps dans les hotspots
    Gaps <- grepl("-", Tableau$Hotspot1) | grepl("-", Tableau$Hotspot2)
    Tableau <- cbind(Tableau, Gaps)

    # Retourner le tableau.
    Tableau
}

# Espace de Test
TableauFKS1 <- ExtractionHotspot("FKS1_filtered.aln")
TableauFKS2 <- ExtractionHotspot("FKS2_filtered.aln")
TableauFKS1$Species[grep("TRUE", TableauFKS1$Gaps)]
TableauFKS2$Species[grep("TRUE", TableauFKS2$Gaps)]
