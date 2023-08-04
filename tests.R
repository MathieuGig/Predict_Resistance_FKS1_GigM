sample <- "BY4741"
gene <- "FKS1"
hotspot <- "HS1"

drug <- c("anidulafungin", "caspofungin", "micafungin")

nameOfFile <- paste(sample, "_", gene, "-", hotspot, "_single_ortho_", sep = "")

File <- paste(getwd(), "processed_data", paste(nameOfFile, drug, sep = ""),
              "selcoeff_all_libraries.csv", sep = "/")

Readfile <- list()

for (i in File)
{
    Readfile <- append(Readfile, list(read.csv(i)))
}

DataAni <- subset(Readfile[[1]], Readfile[[1]]$seq_type == "ortho")
DataCaspo <- subset(Readfile[[2]], Readfile[[2]]$seq_type == "ortho")
DataMica <- subset(Readfile[[3]], Readfile[[3]]$seq_type == "ortho")

## Compute median sel_coeffs.
MedianSelectionCoeffAni <- c()
MedianSelectionCoeffCaspo <- c()
MedianSelectionCoeffMica <- c()

MutationAA <- c()
AASequence <- c()

for (i in DataAni$X)
{
    MedianSelectionCoeffAni <- append(MedianSelectionCoeffAni, DataAni$median_s[DataAni$X == i])
    
    MedianSelectionCoeffCaspo <- append(MedianSelectionCoeffCaspo, DataCaspo$median_s[DataCaspo$X == i])
    
    MedianSelectionCoeffMica <- append(MedianSelectionCoeffMica, DataMica$median_s[DataMica$X == i])
    
    AASequence <- append(AASequence, DataAni$aa_seq[DataAni$X == i])
    MutationAA <- append(MutationAA, DataAni$Nham_aa[DataAni$X == i])
}

MutationAA <- as.factor(MutationAA)

## Make new data frame.
NewData <- data.frame(Anidulafungin = MedianSelectionCoeffAni,
                      Caspofungin = MedianSelectionCoeffCaspo,
                      Micafungin = MedianSelectionCoeffMica,
                      MutatedAminoAcids = MutationAA,
                      Sequence = AASequence)

BestData <- unique(NewData)