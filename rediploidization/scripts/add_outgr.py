import argparse
import sys
from ete3 import Tree

import numpy as np


def load_rectree(input_file, fam):
    fam = fam.split('-')[0]
    print(fam)
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('_reconciliated.nhx:')
            if line[0] == fam:
                tree = Tree(line[1], format=1)
                break
    return tree

def branch_length_closest(tree, gene, group_of_genes):

    """
    Finds the gene closest to `gene` in `tree` amongst a `group_of_genes`, i.e with the shortest
    branch length.
    Args:
        tree (ete3.Tree): input tree
        gene (ete3.TreeNode): gene or node for which to search a neighbour
        group_of_genes (list of str): list of candidate neighbour genes
    Returns:
        str: name of the closest neighbour, in terms of branch lengths
    """

    dist_min = np.inf
    min_gene_d = {}
    names = {}

    #iterate over genes in `group_of_genes` and compute branch-length distance
    for target in group_of_genes:

        dist = tree.get_distance(gene, target)

        if dist <= dist_min:
            dist_min = dist
            min_gene_d[target] = dist_min
            names[target] = target

    #arbitrary choice to ensure deterministic answer
    all_max_genes = []
    for hit_gene in min_gene_d:
        if min_gene_d[hit_gene] == dist_min:
            all_max_genes.append(hit_gene)
    best_gene = sorted(all_max_genes)[0]
    return names[best_gene]


def get_closest_outgr(ctree, tree, outgr, name):

    ctree = Tree(ctree)
    
    genes = {i.name for i in ctree.get_leaves() if i.name != 'outgr'}

    for leaf in tree.get_leaves():
        if leaf.name[-2:] == '.1' and leaf.S == 'Xentro':
            leaf.name = leaf.name[:-2]


    node = tree.get_common_ancestor(genes)

    for k in range(len(outgr)):
        outgr_genes = {i.name for i in tree.get_leaves() if outgr[k] in i.name}
        if outgr_genes:
            break

    if len(outgr_genes) > 1:
        outgr_gene = branch_length_closest(tree, node, outgr_genes)

    elif len(outgr_genes) == 1:
        outgr_gene, = outgr_genes

    # else:
    #     outgr_genes = {i.name for i in tree.get_leaves() if i.name.split('_')[0] in outgr}
    #     if outgr_genes:
    #         outgr_gene = branch_length_closest(tree, node, outgr_genes)
    #     else:
    else:
        sys.stderr.write(f'Error, no outgroup in tree {name}. Exiting.\n')
        return None

    return outgr_gene


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t', '--fulltrees', required=True)

    parser.add_argument('-c', '--ctree', required=True)

    parser.add_argument('--outgr', required=False, default='Bralan,Brabel,Braflo,Strpur,Ptyfla,Sackow,Acapla,Strpur')

    parser.add_argument('-o', '--out', required=False, default='out.txt')

    args = vars(parser.parse_args())

    ctree = args['ctree']

    name = ctree.split('/')[-1].replace('.nh', '')

    fulltree = load_rectree(args['fulltrees'], name)

    outgr_gene = get_closest_outgr(ctree, fulltree, args['outgr'].split(','), name)

    if not outgr_gene:
        sys.exit(1)

    out = ctree.replace('.nh', '.ok.nh')

    t = Tree(ctree)
    for i in t.get_leaves():
        if i.name == 'outgr':
            i.name = outgr_gene
            break

    t.write(outfile=out, format=9)

    with open(args['out'], 'w') as out:
        out.write(outgr_gene+'\n')
