# Create file for NCBI Common Taxonomy Tree tool.
write(c(TableauFKS1$Species, TableauFKS2$Species), file = "allspecies.txt")

# Read NCBI Common Taxonomy Tree output.
MyTree <- readLines("NewTree.txt")

# Format the file so it can be properly used.
MyTreeFormate <- gsub("[|+-]", "", MyTree)
MyTreeFormate <- gsub("\\", "", MyTreeFormate, fixed = TRUE)
MyTreeFormate <- gsub("^ +", "", MyTreeFormate)
MyTreeFormate <- gsub("\\[", "", MyTreeFormate)
MyTreeFormate <- gsub("\\]", "", MyTreeFormate)

# Create list and add the species to it.
TaxoTree <- list(Fungi = list(Blastocladiomycota = MyTreeFormate[4], Basidiomycota = MyTreeFormate[6:17], 
                              Ascomycota = MyTreeFormate[19:257], 
                              Mucoromycota = MyTreeFormate[259:260] ),
                 Viridiplantae = list(Streptophyta = MyTreeFormate[263],
                                      Chlorophyta = MyTreeFormate[265]))

# Specific cases.
TaxoTree$Fungi$Mucoromycota <- append(TaxoTree$Fungi$Mucoromycota, "Rhizopus oryzae")
TaxoTree$Fungi$Basidiomycota <- append(TaxoTree$Fungi$Basidiomycota, "Phanerochaete chrysosporium")
TaxoTree$Fungi$Basidiomycota <- append(TaxoTree$Fungi$Basidiomycota, "Melampsora larici-populina")
TaxoTree$Fungi$Basidiomycota <- append(TaxoTree$Fungi$Basidiomycota, "Melampsora laricipopulina")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Chaetomium thermophilum")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida glabrata")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida castellii")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida bracarensis")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida nivariensis")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Arthrobotrys oligospora")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Rutstroemia sp. NJR-2017a BBW")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Rutstroemia sp. NJR-2017a WRK4")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Marssonina brunnea")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Byssochlamys spectabilis")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Penicillium sp. 'occitanis'")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Talaromyces cellulolyticus")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Leptosphaeria maculans")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Pyrenophora tritici-repentis")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Acidomyces richmondensis")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Torrubiella hemipterigena")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Saccharomycetaceae sp. 'Ashbya aceri'")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Saccharomyces sp. 'boulardii'")

AddTaxonomy <- function(Vector_of_species, dataframe, TaxonomyTree = TaxoTree)
{
    # Create vectors
    Kingdom <- vector()
    Phylum <- vector()
    
    # Find the kingdom and phylum of every specie in the data.
    for (i in Vector_of_species)
    {
        if (i %in% TaxonomyTree$Fungi$Blastocladiomycota)
        {
            Kingdom <- append(Kingdom, "Fungi")
            Phylum <- append(Phylum, "Blastocladiomycota")
            next
        }
        if (i %in% TaxonomyTree$Fungi$Basidiomycota)
        {
            Kingdom <- append(Kingdom, "Fungi")
            Phylum <- append(Phylum, "Basidiomycota")
            next
        }
        if (i %in% TaxonomyTree$Fungi$Ascomycota)
        {
            Kingdom <- append(Kingdom, "Fungi")
            Phylum <- append(Phylum, "Ascomycota")
            next
        }
        if (i %in% TaxonomyTree$Fungi$Mucoromycota)
        {
            Kingdom <- append(Kingdom, "Fungi")
            Phylum <- append(Phylum, "Mucoromycota")
            next
        }
        if (i %in% TaxonomyTree$Viridiplantae$Streptophyta)
        {
            Kingdom <- append(Kingdom, "Viridiplantae")
            Phylum <- append(Phylum, "Streptophyta")
            next
        }
        if (i %in% TaxonomyTree$Viridiplantae$Chlorophyta)
        {
            Kingdom <- append(Kingdom, "Viridiplantae")
            Phylum <- append(Phylum, "Chlorophyta")
            next
        }
        else
        {
            Kingdom <- append(Kingdom, "-")
            Phylum <- append(Phylum, "-")
            next
        }
    }

    # Add to Kingdom and Phylum to FKS data frame.
    dataframe <- cbind(dataframe, Kingdom)
    dataframe <- cbind(dataframe, Phylum)

    # Set column order.
    col_order <- c("Identifier", "Species", "Kingdom", "Phylum",
                   "Hotspot1", "Hotspot2", "GapsHotspot1", "GapsHotspot2",
                   "HS1_in_oPools", "HS2_in_oPools")

    # Reorganize column order.
    dataframe <- dataframe[, col_order]

    dataframe
}


# Is human pathogen ? Basé sur le WHO fungal priority pathogens
# list to guide research, development and public health action.
IsPathogen <- function(Species)
{
    return (Species %in% c('Cryptococcus neoformans', 'Candida auris',
                           'Aspergillus fumigatus', 'Candida albicans',
                           'Nakaseomyces glabrata', 'Candida glabrata',
                           'Histoplasma spp.', 'Mucorales', 'Fusarium spp.',
                           'Candida tropicalis', 'Candida parapsilosis',
                           'Scedosporium spp.', 'Cryptococcus gattii',
                           'Lomentospora prolificans', 'Talaromyces marneffei',
                           'Coccidiodes spp.', 'Pneumocystis jirovecii',
                           'Pichia kudriavzeveii', 'Candida krusei',
                           'Paracoccidioides spp.'))
}


# Testing
TableauFKS1 <- AddTaxonomy(TableauFKS1$Species, TableauFKS1)
TableauFKS2 <- AddTaxonomy(TableauFKS2$Species, TableauFKS2)

CriticalHumanPathogen1 <- unlist(lapply(TableauFKS1$Species, IsPathogen))
CriticalHumanPathogen2 <- unlist(lapply(TableauFKS2$Species, IsPathogen))

TableauFKS1 <- cbind(TableauFKS1, Is_Human_Pathogen = CriticalHumanPathogen1)
TableauFKS2 <- cbind(TableauFKS2, Is_Human_Pathogen = CriticalHumanPathogen2)