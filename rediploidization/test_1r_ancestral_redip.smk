#`conda activate redip`
import os

from ete3 import Tree

PWD = os.getcwd()

RAXML_SEED = 8080

ali = 'data/ali'
trees = 'data/brofams_0221_rec.trees'
genes_to_clg = 'data/OG_CLG_table.txt'

out = 'output_1r_rediptest'

rule all:
    input: '.end'


wildcard_constraints:
    fam='OG_[0-9]+.CLG[A-Z1-9]*',
    i='[0-9]+'

checkpoint prepare_aore_lore_ctrees:
    input:  genes_to_clg
    output:
        a = directory(f'{out}/aore_ctrees'),
        l = directory(f'{out}/lore_ctrees')
    shell: "python scripts/get_informative_fams.py --input {input} -oa {output.a} -ol {output.l}"

rule add_outgr:
    input:
        lore = f'{out}/lore_ctrees/{{fam}}.nh',
        trees = trees
    output:
        lore = f'{out}/lore_ctrees/{{fam}}.ok.nh',
        outgr = f'{out}/outgr/{{fam}}.txt'
    shell: 'python scripts/add_outgr.py --fulltrees {input.trees} --ctree {input.lore} --out {output.outgr}'


rule add_outgr_aore:
    input: f'{out}/aore_ctrees/{{fam}}-{{i}}.nh', f'{out}/outgr/{{fam}}.txt'
    output: f'{out}/aore_ctrees/{{fam}}.{{i}}.ok.nh'
    run:
        with open(input[1], 'r') as infile:
            for line in infile:
                o = line.strip()
        t = Tree(input[0])
        for i in t.get_leaves():
            if i.name == 'outgr':
                i.name = o
                break

        t.write(outfile=output[0], format=9)

def get_ali(ali_folder, family):
    family = family.split('-')[0]
    return f'{ali_folder}/{family}.al.cl.fa'

rule extract_subali:
    input: t = f'{out}/lore_ctrees/{{fam}}.ok.nh', a = lambda w: get_ali(ali, w.fam)
    output: f'{out}/ali/{{fam}}.fa'
    shell: 'python scripts/get_subali.py --tree {input.t} --ali {input.a} --output {output}'


rule check_ali:
    input: f'{out}/ali/{{fam}}.fa'
    output: f'{out}/ali/{{fam}}.ok.fa'
    params: s = RAXML_SEED
    shell:
        "rm {out}/ali/RAxML_info.{wildcards.fam} || true && "
        "raxmlHPC -f c --print-identical-sequences -n {wildcards.fam} -m PROTGAMMAJTT "
        "-s {input} -w {PWD}/{out}/ali/ -p {params.s} >&2 && "
        "if [ ! -s {out}/ali/{wildcards.fam}.fa.reduced ]; "
        "then cp {input} {output}; else mv {out}/ali/{wildcards.fam}.fa.reduced {output};fi"


rule build_ml_constrained_lore:
    input: ctree = f'{out}/lore_ctrees/{{fam}}.ok.nh', ali =  f'{out}/ali/{{fam}}.ok.fa'
    output:  f'{out}/lore/{{fam}}.nh'
    params: s = RAXML_SEED
    shell:
        "rm {out}/ali/RAxML_info.{wildcards.fam}_lore || true && "
        "raxmlHPC -N 10 -g {input.ctree} -p {params.s} -n {wildcards.fam}_lore -m PROTGAMMAJTT "
        "-s {input.ali} -w {PWD}/{out}/ali/ && "
        "mv {out}/ali/RAxML_bestTree.{wildcards.fam}_lore {output}"       


rule build_ml_constrained_aore:
    input: ctree = f'{out}/aore_ctrees/{{fam}}.{{i}}.ok.nh', ali =  f'{out}/ali/{{fam}}.ok.fa'
    output:  f'{out}/aore/{{fam}}.{{i}}.nh'
    params: s = RAXML_SEED
    shell:
        "rm {out}/ali/RAxML_info.{wildcards.fam}_aore_{wildcards.i} || true && "
        "raxmlHPC -N 10 -g {input.ctree} -p {params.s} -n {wildcards.fam}_aore_{wildcards.i} -m PROTGAMMAJTT "
        "-s {input.ali} -w {PWD}/{out}/ali/ && "
        "mv {out}/ali/RAxML_bestTree.{wildcards.fam}_aore_{wildcards.i} {output}" 


rule build_ml_tree:
    input: ali = f'{out}/ali/{{fam}}.ok.fa'
    output: f'{out}/ml/{{fam}}.nh'
    params: s = RAXML_SEED
    shell:
        "rm {out}/ali/RAxML_info.{wildcards.fam}_ml || true && "
        "raxmlHPC -N 10 -p {params.s} -n {wildcards.fam}_ml -m PROTGAMMAJTT -s {input.ali} -w {PWD}/{out}/ali/ && "
        "mv {out}/ali/RAxML_bestTree.{wildcards.fam}_ml {output}"


def get_all_aore(wildcards, fam):
  co = checkpoints.prepare_aore_lore_ctrees.get(**wildcards).output[0]
  aore_ids, = glob_wildcards(f'{out}/aore_ctrees/{fam}-{{i}}.nh')
  return expand(f'{out}/aore/{fam}.{{i}}.nh', i=aore_ids)



rule lk_test:
    """
    Runs the likelihood AU test to compare likelihoods of the 3 constrained topologies and ML trees. 
    """
    input:
        ali = f'{out}/ali/{{fam}}.ok.fa', ml = f'{out}/ml/{{fam}}.nh', lore = f'{out}/lore/{{fam}}.nh',
        aore = lambda w: get_all_aore(w, w.fam)
    output: touch(f'{out}/lktest/{{fam}}.txt')
    shell:
        "bash scripts/prototype_au_test_many.sh {wildcards.fam} {input.ml} {input.ali} {input.lore} "
        "'{input.aore}' {output} PROTGAMMAJTT  && rm {out}/ali/RAxML_info.{wildcards.fam}_lktest"


def lktests(wildcards):
    co = checkpoints.prepare_aore_lore_ctrees.get(**wildcards).output[1]
    ids, = glob_wildcards(f'{out}/lore_ctrees/{{fam}}.nh')
    ids = [i for i in ids if '.' not in i.split('-')[0] and 'ok' not in i]
    return expand(f'{out}/lktest/{{fam}}.txt', fam=ids)


rule end:
    input: lktests
    output: touch('.end')