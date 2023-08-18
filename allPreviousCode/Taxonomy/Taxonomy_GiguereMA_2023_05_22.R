write(VectorSpecies1, file = "Species.txt")

MyTree <- readLines("TaxoTree.txt")

MyTreeFormate <- gsub("[|+-]", "", MyTree)
MyTreeFormate
MyTreeFormate <- gsub("\\", "", MyTreeFormate, fixed = TRUE)
MyTreeFormate <- gsub("^ +", "", MyTreeFormate)
MyTreeFormate

TaxoTree <- list(Fungi = list(Basidiomycota = MyTreeFormate[3:12], 
                              Ascomycota = MyTreeFormate[14:197], 
                              Mucoromycota = MyTreeFormate[199:200] ),
                 Viridiplantae = list(Streptophyta = MyTreeFormate[204],
                                      Chlorophyta = MyTreeFormate[206]))

# CasSpécifiques
TaxoTree$Fungi$Mucoromycota <- append(TaxoTree$Fungi$Mucoromycota, "Rhizopus oryzae")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Chaetomium thermophilum")
TaxoTree$Fungi$Basidiomycota <- append(TaxoTree$Fungi$Basidiomycota, "Phanerochaete chrysosporium")
TaxoTree$Fungi$Basidiomycota <- append(TaxoTree$Fungi$Basidiomycota, "Melampsora larici-populina")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida glabrata")

# Création de liste
ListKingdom <- vector()
ListPhylum <- vector()

for (i in VectorSpecies1)
{
    if (i %in% TaxoTree$Fungi$Basidiomycota)
    {
        ListKingdom <- append(ListKingdom, "Fungi")
        ListPhylum <- append(ListPhylum, "Basidiomycota")
        next
    }
    if (i %in% TaxoTree$Fungi$Ascomycota)
    {
        ListKingdom <- append(ListKingdom, "Fungi")
        ListPhylum <- append(ListPhylum, "Ascomycota")
        next
    }
    if (i %in% TaxoTree$Fungi$Mucoromycota)
    {
        ListKingdom <- append(ListKingdom, "Fungi")
        ListPhylum <- append(ListPhylum, "Mucoromycota")
        next
    }
    if (i %in% TaxoTree$Viridiplantae$Streptophyta)
    {
        ListKingdom <- append(ListKingdom, "Viridiplantae")
        ListPhylum <- append(ListPhylum, "Streptophyta")
        next
    }
    if (i %in% TaxoTree$Viridiplantae$Chlorophyta)
    {
        ListKingdom <- append(ListKingdom, "Viridiplantae")
        ListPhylum <- append(ListPhylum, "Chlorophyta")
        next
    }
    else
    {
        ListKingdom <- append(ListKingdom, "-")
        ListPhylum <- append(ListPhylum, "-")
        next
    }
}
TaxoTableau <- data.frame(VectorSpecies1, ListKingdom, ListPhylum)