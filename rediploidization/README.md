# Tests of ancestral and lineage-specific rediploidisation gene tree topologies

## Post-1R rediploidisation

### Data and results inspection

Multiple sequence alignments for each of the 1247 gene families considered in the analysis are provided in `results_redip_1r/ali/`.

Multifurcated constrained tree topologies under the ancestral rediploidisation and lineage-specific rediploidisation models are provided in `results_redip_1r/` (`aore_ctrees/` and `lore_ctrees/` subfolders, respectively).

Tree topologies obtained for the unconstrained ml trees, ancestral rediploidisation trees and lineage-specific rediploidisation trees are provided in `results_redip_1r/` (`ml/`, `aore/` and `lore/` subfolders, respectively.)

The results from likelihood AU tests are provided in the `results_redip_1r/lktests/` folder. These files correspond to raw outputs fron the CONSEL package (Shimodaira, H. & Hasegawa, 2001), where the ml tree is identified by the item id '1', lore tree by the item id '2' and aore trees by item ids starting from '3'.

### Reproducibility

The following commands allow to reproduce the presented analysis (barplot in Figure 3b):

`mamba env create -f redip.yaml`

`conda activtate redip`

`snakemake -s test_1r_ancestral_redip.smk --cores 50 -k`

`python scripts/parse_au.py -i output_1r_rediptest/lktests/`

`python scripts/plot_au.py -i out.pkl`

## Post-2R cyclostomes CLGB rediploidisation 

### Data and results inspection

Multiple sequence alignments for each of the 30 gene families considered in the analysis are provided in `results_redip_clgb_2rcyclo/ali/`.

Multifurcated constrained tree topologies under the ancestral rediploidisation and lineage-specific rediploidisation models are provided in `results_redip_clgb_2rcyclo/ctrees/` (`aore/` and `lore/` subfolders, respectively).

Tree topologies obtained for the unconstrained ml trees, ancestral rediploidisation trees and lineage-specific rediploidisation trees are provided in `results_redip_clgb_2rcyclo/` (`ml/`, `aore/` and `lore/` subfolders, respectively). 

The results from the likelihood AU tests are provided in the `results_redip_1r/lktests/` folder. These files correspond to raw outputs fron the CONSEL package (Shimodaira, H. & Hasegawa, 2001), where the unconstrained ml tree is identified by the item id '1', the lore tree by the item id '2' and the aore tree by the item id '3'.

### Reproducibility

The following commands allow to reproduce the gene tree topology tests for lineage-specific versus ancestral rediploidisation for genes descended from the ancestral chromosome CLG B 1R copy 1.

`mamba env create -f redip.yaml`

`conda activtate redip`

`snakemake -s test_2rcy_clgb1.smk --cores 50 -k`

`python scripts/parse_au.py -i output_cyclo-CLGB_rediptest/lktests/`

`python scripts/plot_au.py -i out.pkl`