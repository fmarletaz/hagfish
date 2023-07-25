# WGDs inference on the vertebrate phylogeny

## Reproducing the analysis

- Uncompress the `data_all` folder (which contains conditional clade probabilities for the 8931 gene families).

- Open an interactive julia session, using for instance 40 threads: `julia -t 40`

- Run the commands from the `whale_brspecDLWGD.jl` script.


Dependencies:

- julia 1.8.1
- WHALE v2.1.0 `https://github.com/arzwa/Whale.jl.git#master`

## WHALE Reference

Arthur Zwaenepoel & Yves Van de Peer. Inference of Ancient Whole-Genome Duplications and the Evolution of Gene Duplication and Loss Rates. Molecular Biology and Evolution. 2019.  https://doi.org/10.1093/molbev/msz088
