
# Transcriptomics

Gene expression were estimated as FPKM for a set of organs in lamprey (Ljap), gar (Locu) and hagfish (Eata) in `*_organs_fpkms.tsv`. Gene expression specificity (tau) was calculated from these values using `tau_re.r`. Their distribution for ohnologues described in was estimated incorporating the gene family information reported in the files `*_ohnofams.txt` for each species. 

Gene expression pattern gain and losses were computed from FPKM using `subfunc.r`. 

The gene-cluster associated from the WGCNA analysis is also reported in `WGCNA`. 