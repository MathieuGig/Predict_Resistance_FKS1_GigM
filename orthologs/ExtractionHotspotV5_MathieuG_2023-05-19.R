###
### ExtractionHotspot(fichierAlignements)
###
##
## Find and extract hotspots 1 and 2 for each sequences of the FKS genes. 
## Returns a data frame with the hotspot sequences, the species and their identifier
## and other informations.
##
##
## Arguments
##
## fichierAlignements: character string, the name of a *.aln file.
##
##
## Value
##
## A data frame of 8 columns describing the species, their identifier,
## the hotspots sequences of these species, the presence of gaps in the hotspots,
## and the presence of the hotspot in the oPools.
##
##
## Example
##
## ExtractionHotspot("FKS1_filtered.aln")
##

ExtractionHotspot <- function(fichierAlignements)
{
    # Validation of the input file.
    stopifnot("Invalid file for analysis. Require a *.aln file" = grepl(".aln", fichierAlignements))
    
    # Read the file of filtered sequences.
    SequencesFiltrees <- readLines(fichierAlignements)
    
    # Extract and Format identifiers using des REGEX.
    IdentifiantFormate <- gsub("\\|.*", "", grep(">", SequencesFiltrees, value = TRUE))
    IdentifiantFormate <- gsub(" .*", "", IdentifiantFormate)
    IdentifiantFormate <- gsub(">", "", IdentifiantFormate)
    
    # Extract and Format the species using REGEX.
    EspecesNonFormate <- gsub("^.*\\[", "", grep("\\[*\\]", SequencesFiltrees, value = TRUE))
    EspecesFormate <- gsub("\\]", "", EspecesNonFormate)
    
    # Create data frame.
    Tableau <- data.frame(Identifier = IdentifiantFormate, Species = EspecesFormate)
    
    # Create hotspot lists.
    ListeHotspots1 <- list()
    ListeHotspots2 <- list()

    # Find the sequences corresponding to the Hotspots using their fixed positions.
    # For each header line (so for each sequence)
    for (i in grep(">", SequencesFiltrees))
    {
        # if FKS1
        if (grepl("FKS1", fichierAlignements))
        {
            # Add to the Hotspot list the 8 or 9 amino acids at these fixed positions.
            ListeHotspots1 <- append(ListeHotspots1, substr(SequencesFiltrees[i+17], start = 52, stop = 60))
            ListeHotspots2 <- append(ListeHotspots2, substr(SequencesFiltrees[i+38], start = 1, stop = 8))   
        }
        # if FKS2
        else
        {
            # Add to the Hotspot list the 8 or 9 amino acids at these fixed positions.
            # FKS2 Hotspot1 is 1 row farther and 1 aa smaller than the FKS1 Hotspot1.
            # FKS2 Hotspot2 is farther on its line than the FKS1 Hotspot2.
            ListeHotspots1 <- append(ListeHotspots1, substr(SequencesFiltrees[i+18], start = 52, stop = 59))
            ListeHotspots2 <- append(ListeHotspots2, substr(SequencesFiltrees[i+38], start = 33, stop = 40))
        }
    }

    # Add Hotspots columns to the data frame.
    Tableau <- cbind(Tableau, Hotspot1 = unlist(ListeHotspots1))
    Tableau <- cbind(Tableau, Hotspot2 = unlist(ListeHotspots2))

    ## Adding to the data frame columns indicating gaps in the hotspots
    GapsHotspot1 <- grepl("-", Tableau$Hotspot1)
    GapsHotspot2 <- grepl("-", Tableau$Hotspot2)
    Tableau <- cbind(Tableau, GapsHotspot1)
    Tableau <- cbind(Tableau, GapsHotspot2)

    ## Adding to the data frame columns indicating the presence of 
    ## hotspots sequences in the oPools.

    # Read the files carrying the sequences representing the oligos.
    FKS1_Hotspots1_uniques <- (readLines("FKS1_HS1_orthologs_unique.fa"))[-1]
    FKS1_Hotspots2_uniques <- (readLines("FKS1_HS2_orthologs_unique.fa"))[-1]
    FKS2_Hotspots1_uniques <- (readLines("FKS2_HS1_orthologs_unique.fa"))[-1]
    FKS2_Hotspots2_uniques <- (readLines("FKS2_HS2_orthologs_unique.fa"))[-1]
    
    # Create lists
    HS1oPools <- list()
    HS2oPools <- list()

    # if FKS1
    if (grepl("FKS1", fichierAlignements))
    {
        # For each hotspot sequence in the data frame, verify its existence in
        # the sequences used for the oligos, then add this logical value 
        # to a list.
        for (i in Tableau$Hotspot1)
        {
            HS1oPools <- append(HS1oPools, i %in% FKS1_Hotspots1_uniques)
        }
        for (i in Tableau$Hotspot2)
        {
            HS2oPools <- append(HS2oPools, i %in% FKS1_Hotspots2_uniques)
        }
    }

    # if FKS2
    else
    {
        for (i in Tableau$Hotspot1)
        {
            HS1oPools <- append(HS1oPools, i %in% FKS2_Hotspots1_uniques)
        }
        for (i in Tableau$Hotspot2)
        {
            HS2oPools <- append(HS2oPools, i %in% FKS2_Hotspots2_uniques)
        }
    }

    # Add columns to data frame.
    Tableau <- cbind(Tableau, HS1_in_oPools = unlist(HS1oPools))
    Tableau <- cbind(Tableau, HS2_in_oPools = unlist(HS2oPools))

    # Return data frame.
    Tableau
}

### Testing
TableauFKS1 <- ExtractionHotspot("FKS1_filtered.aln")
TableauFKS2 <- ExtractionHotspot("FKS2_filtered.aln")