

# Analysis of the brown hagfish *Eptatretus atami*

This repository contains the code employed to investigate the hagfish genome, as well as key data files. 

## Key supplementary data files

 * Supplementary File 1. Table gathering essential information for each gene model in E. atami including location, protein domains, gene family, gene expression cluster [Eptata_genes_filt.xlsx]

 * Supplementary File 3. Phylogenetic trees inferred for paralogons in each CLG assuming the C20+R
 model [paralogons_t6_c20.pdf]

 * Supplementary File 3. Phylogenetic tree of cadherin-related genes showing the tandem duplication
 patterns in gnathostomes [OG_4678_Cadh.tre.pdf]
 * Supplementary File 4. Functional enrichment for sets of paralogues showing distinct retention
 pattern after genome duplications in vertebrates (Figure 4b). [Vert2R_Go_enrich_wg_gds.txt].
 * Supplementary File 5. Functional enrichment in lamprey for gene families lost in hagfish (data for
 Figure 4e). [Loss_Lamprey_GO_gn_dsc.tsv]
* Supplementary File 6. Synteny-based paralogue classification for reconstructed gene families
 
## Synteny 

The folder contains multiple files describing the reciprocal coordinates for a number of genome comparisons based on mutual-best-hits performed using `mmseqs2`. Binned stacked barplot representation such as that in Figure 2a are based on the `gchkExp` function in `plot_clg_hagfish.r`, while direct plotting of individual gene linkage between chromosomes such as that of Figure 2b is performed using `Rideogram` such as described in `ideogram.r`. 

## 