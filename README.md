

# Analysis of the brown hagfish *Eptatretus atami*

This repository contains the code employed to investigate the hagfish genome, as well as key data files. 

## Key supplementary data files

 * Supp. File 1. Table gathering essential information for each gene model in E. atami including location, protein domains, gene family, gene expression cluster [Eptata_genes_filt.xlsx]()

 * Supp. File 3. Phylogenetic trees inferred for paralogons in each CLG assuming the C20+R
 model [paralogons_t6_c20.pdf](https://github.com/fmarletaz/hagfish/blob/main/Paralogons/paralogons_t6_c20.pdf)

 * Supp. File 3. Phylogenetic tree of cadherin-related genes showing the tandem duplication
 patterns in gnathostomes [OG_4678_Cadh.tre.pdf](https://github.com/fmarletaz/hagfish/blob/main/)
 * Supp. File 4. Functional enrichment for sets of paralogues showing distinct retention
 pattern after genome duplications in vertebrates (Figure 4b). [Vert2R_Go_enrich_wg_gds.txt](https://github.com/fmarletaz/hagfish/blob/main/Functional/Vert2R_Go_enrich_long_wg_gds.txt).
 * Supp. File 5. Functional enrichment in lamprey for gene families lost in hagfish (data for
 Figure 4e). [Loss_Lamprey_GO_gn_dsc.tsv](https://github.com/fmarletaz/hagfish/blob/main/Functional/Loss_Lamprey_GO_gn_dsc.tsv)
* Supp. File 6. Synteny-based paralogue classification for reconstructed gene families [Vert_Evt_OGrrA.txt](https://github.com/fmarletaz/hagfish/blob/main/Paralogons/Vert_Evt_OGrrA.txt)
 
## Synteny 

The folder contains multiple files describing the reciprocal coordinates for a number of genome comparisons based on mutual-best-hits performed using `mmseqs2`. Binned stacked barplot representation such as that in Figure 2a are based on the `gchkExp` function in `plot_clg_hagfish.r`, while direct plotting of individual gene linkage between chromosomes such as that of Figure 2b is performed using `Rideogram` such as described in `ideogram.r`. 

## Paralogons

The assembly of the concatenated paralogons was performed using the code in `Tetraploidy_cl.ipynb` and information gathered from several data files including reconciled tree files `brofams_rec.trees` and chromosome annotation in different species `geneInfo`. The file including paralogon taxa content (e.g. `paralogons_genes_rsGrRx2.txt`) is then converted into an alignment using the script `sel-pgon2.py`. Results of 'strict' (`test5`) and 'relaxed' (`test6`) analyses are also provided in corresponding folders. 

## Gene families

The code to perform gene gain, loss and duplication counts is provided in `gene_families.ipynb` along with necessaey files. Some selected reconciliated gene trees related to neural crest specification are included as well. 

## Functional 

Functional enrichment (Gene Ontology) was performed on lists of paralogues showing particular retention in regard to the vertebrate WGD in `ext_go.r`. The same file also contains GO enrichment analysis for gene families lost in hagfish. 

## Transcriptomics

Gene expression specificity (tau) was calculated using `tau_re.r` and gene expression pattern gain and losses using `subfunc.r`. Gene family and expression files required are also provided. 