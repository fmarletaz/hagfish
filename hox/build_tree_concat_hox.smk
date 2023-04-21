# snakemake workflow to build phylogenies from concatenated hox clusters
# requires snakemake, python3, mafft and raxml-ng
# usage: snakemake -s build_tree_concat_hox.smk --cores 30 

GENES_ORDER = ['Copz', 'Nfe2', 'Hnrnpa', 'Cbx', 'Skap', 'Smug', 'Mtx2',
			   'Hox1', 'Hox2', 'Hox3', 'Hox4', 'Hox5', 'Hox6', 'Hox7', 'Hox8', 'Hox9', 'Hox10', 'Hox11', 'Hox12','Hox13','Hox14',
				'Evx', 'Lnp', 'Calcoco', 'Jazf', 'ATP5MC', 'Creb']

GENES_ORDER = [i.upper() for i in GENES_ORDER]

rule Target:
	input: 'concat_hox_clusters/hox_clusters_concat.raxml.bestTree'

rule mafft:
	input: 'concat_hox_clusters/seq/sequences_{gene}.fasta'
	output: 'concat_hox_clusters/ali/msa_{gene}.fasta'
	shell: 'mafft {input} > {output}'

rule concat_ali:
	input: expand('concat_hox_clusters/ali/msa_{gene}.fasta', gene=GENES_ORDER)
	output: 'concat_hox_clusters/concat_ali.fasta'
	shell: 'python concat_ali.py --input {input} --output {output}'

rule raxml:
	input: 'concat_hox_clusters/concat_ali.fasta'
	output: 'concat_hox_clusters/hox_clusters_concat.raxml.bestTree'
	params: pref = 'concat_hox_clusters/hox_clusters_concat', seed = '1671635227' #fix seed for reproducibility
	threads: 100
	shell: 'raxml-ng --msa {input} --model LG+G4+F --all --tree pars{{10}} --bs-trees 100 --prefix {params.pref} --seed {params.seed}'