

# Code and data for reconstruction of paralogons trees

The assembly of the concatenated paralogons was performed using the code in `Tetraploidy_cl.ipynb` and information from several data files including reconciled tree files `brofams_rec.trees` and chromosome annotation in different species `geneInfo`. Essentially for each gene family, it assigns its CLG and the duplication status of the different genes and species. A cutoff on the number of paralogons recovered in each gene family in gnathostomes is applied (at least 3 paralogon: `strict` or at least 2 paralogons: `relaxed`). 
The file including paralogon taxa content (e.g. `paralogons_genes_rsGrRx2.txt`) is then converted into an alignment using the script `sel-pgon2.py`. Alignments and resulting trees are included in the folder  `strict/alis_strict` or `relaxed/alis_relax` as well as resulting trees. 
For molecular dating, a custom set of calibration was generated for each CLG using speciation nodes as inferred using the 
`ete3` library using a set of calibration between species (`calibrations2r.txt`). Files and results are located in the `molecular dating` folder. Divergence dates for given nodes were extracted using the code in `plot_chrono.r` and summarised in `dat_res_clg.txt`. 

Dependencies: 
- (ete3)[http://etetoolkit.org/]