# Code and data for the vertebrate hox phylogeny (Figure S)

## Re-running the analysis

Alignments, concatenation and phylogeny inference can be re-run in a single command: `snakemake -s build_tree_concat_hox.smk --cores 30` (or `snakemake -s build_tree_concat_hox.smk --cores 30 -np` to print the list of commands without running them).

Dependencies:
- python (3.10)
- snakemake (7.17)
- mafft (v7.508)
- raxml-ng (v1.1)

## Data

We make input and output data available so that re-running the code is not necessary, see details below.

### Input

Sequences of hox and bystander genes in selected species are available in `concat_hox_clusters/seq/`
Gene names are preprended with species names: Braflo for amphioxus, Parata for hagfish, Petmar for lamprey, Lepocu for spotted gar, Galgal for chicken, Homsap for human and Musmus for mouse.

### Output

The concatenated alignment used for phylogenetic inference is available in `output/concatenated_hox_bystander_msa.fasta`
The reconstructed raxml tree is in: `output/concatenated_hox_bystander_msa.support`