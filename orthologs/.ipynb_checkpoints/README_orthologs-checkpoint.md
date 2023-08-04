# Ortholog sequences (FKS HS DMS)

## Sequence retrieval procedure

FKS1 and FKS2 orthologous sequences were retrieved from [Metaphors](http://orthology.phylomedb.org/).

### FKS1

1. "FKS1" search in Metaphors
2. **M!42720474_ASPNG** (*Aspergillus niger*) was used as query --> 290 hits
3. Sequences were filtered based on length (need to be > 1,500 aa and < 1,970 aa). This corresponds to > 80% length of FKS1 from *S. cerevisiae* (1,877 aa) and < 105% length
4. Filtered sequences were aligned with MUSCLE and output in Pearson/FASTA format (> .aln)
5. Aligned sequences were trimmed using trimal (automated1 mode)
6. Hotspots were manually selected from the alignment
7. Unique hotspot sequences were curated manually

### FKS2

1. "FKS2" search in Metaphors
2. **M!115416204_CANGB** (*Candida glabrata*) was used as query
3. Sequences were filtered based on length (need to be > 1,517 aa and < 1,991 aa). This corresponds to > 80% length of FKS1 from *S. cerevisiae* (1,896 aa) and < 105% length
4. Three sequences from the same species were filtered out because they contained multiple stops throughout