import argparse
import os

from ete3 import Tree
import itertools


#homologs bootstrap > 70: for each amc chr, there are two groups of cyclostomata chromosomes that have strong support to be grouped together
#do not separate these when making the aore topologies
KEEP_TOGETHER = { 'CLGA1' : [{'Petmar_chrm17', 'Petmar_chrm33', 'Petmar_chrm53', 'Parata_chr4', 'Parata_chr10', 'Parata_chr1'}, {'Parata_chr14', 'Parata_chr13', 'Petmar_chrm31', 'Petmar_chrm29', 'Petmar_chrm6', 'Parata_chr9'}],
                  'CLGB' : [{'Petmar_chrm2', 'Parata_chr8', 'Petmar_chrm10', 'Parata_chr4', 'Parata_chr13', 'Petmar_chrm56'},
                             {'Petmar_chrm4', 'Parata_chr5', 'Petmar_chrm19', 'Parata_chr1', 'Petmar_chrm34'},
                              #{'Petmar_chrm4', 'Parata_chr5', 'Petmar_chrm19', 'Parata_chr1', 'Petmar_chrm34'}
                            ],
                  'CLGC' : [{'Petmar_chrm40', 'Parata_chr5', 'Parata_chr4', 'Petmar_chrm13', 'Parata_chr11', 'Parata_chr10', 'Petmar_chrm52'},
                             {'Petmar_chrm57', 'Parata_chr8', 'Petmar_chrm1', 'Parata_chr12'}],

                  'CLGD' : [{'Petmar_chrm25', 'Parata_chr10', 'Parata_chr4', 'Petmar_chrm63', 'Parata_chr9', 'Petmar_chrm45'},
                             {'Petmar_chrm54', 'Parata_chr8'}, {'Petmar_chrm20', 'Parata_chr11'}, {'Parata_chr14','Petmar_chrm69'}],

                  'CLGE' : [{'Petmar_chrm30', 'Parata_chr1', 'Parata_chr10', 'Petmar_chrm48', 'Parata_chr10', 'Petmar_chrm14'},
                             {'Petmar_chrm42', 'Parata_chr13', 'Petmar_chrm36', 'Petmar_chrm46', 'Parata_chr3'}],

                  'CLGF' : [{'Petmar_chrm61', 'Parata_chr16', 'Parata_chr17',  'Parata_chr13', 'Petmar_chrm18'},
                             {'Petmar_chrm16', 'Parata_chr3'}],

                  'CLGG' : [{'Petmar_chrm38', 'Parata_chr4', 'Parata_chr2'},
                             {'Petmar_chrm5', 'Parata_chr16'}],

                  'CLGH' : [{'Petmar_chrm43', 'Parata_chr17', 'Parata_chr8', 'Petmar_chrm26', 'Parata_chr2', 'Petmar_chrm3'},
                             {'Parata_chr7', 'Parata_chr1', 'Petmar_chrm37', 'Petmar_chrm58', 'Parata_chr14'}],

                  'CLGI' : [{'Petmar_chrm11', 'Parata_chr2'},
                             {'Petmar_chrm9', 'Parata_chr12'}],

                  'CLGJ' : [{'Petmar_chrm12', 'Parata_chr10', 'Parata_chr15', 'Petmar_chrm15'},
                             {'Parata_chr6', 'Parata_chr11', 'Petmar_chrm23', 'Petmar_chrm27', 'Parata_chr2', 'Petmar_chrm3'}],

                  'CLGK' : [{'Petmar_chrm51', 'Parata_chr2', 'Parata_chr14', 'Petmar_chrm60'},
                             {'Parata_chr16', 'Petmar_chrm49'}],

                  'CLGL' : [{'Petmar_chrm44', 'Parata_chr11', 'Parata_chr1', 'Petmar_chrm1'},
                             {'Parata_chr7', 'Petmar_chrm35'}],

                  'CLGM' : [{'Petmar_chrm59', 'Parata_chr15', 'Parata_chr2', 'Petmar_chrm50', 'Petmar_chrm74'},
                             {'Parata_chr6', 'Petmar_chrm22'}],

                  'CLGN' : [{'Petmar_chrm47', 'Parata_chr17', 'Parata_chr6', 'Petmar_chrm64', 'Parata_chr9'},
                             {'Petmar_chrm67', 'Parata_chr15', 'Parata_chr1', 'Petmar_chrm55'}, {'Petmar_chrm70', 'Parata_chr3'}],

                  'CLGO' : [{'Petmar_chrm24', 'Parata_chr6','Petmar_chrm21', 'Parata_chr15'},
                            {'Petmar_chrm65', 'Parata_chr5', 'Petmar_chrm7', 'Parata_chr7', 'Parata_chr1'}],

                  'CLGP' : [{'Petmar_chrm28', 'Parata_chr3'}, {'Petmar_chrm39', 'Parata_chr4'}],

                  'CLGQ' : [{'Petmar_chrm8', 'Parata_chr6'}, {'Petmar_chrm41', 'Parata_chr3'}],
}


def check_fam(fam): 

    dup = {i[1] for i in fam}

    if 'False' in dup:
        return False

    #needs at leat one 1R_1 and one 1R_2 in vertebrates
    dup = {i[1].replace('alpha', '').replace('beta', '') for i in fam if 'alpha' in i[1] or 'beta' in i[1]}
    if len(dup) != 2:
        return False

    #AT least one gene for hagfish and lamprey
    sp =  {i[2] for i in fam}
    if 'Parata' not in sp or 'Petmar' not in sp:
        return False

    # at least one vertebrates with two copies
    vert = {(i[2], i[1]) for i in fam if i[2] not in ['Petmar', 'Parata']}

    vert_1 = {i[0] for i in vert if i[1].replace('alpha', '').replace('beta', '')=='1'}#1R_1
    vert_2 = {i[0] for i in vert if i[1].replace('alpha', '').replace('beta', '')=='2'}#1R_2

    if not vert_1.intersection(vert_2):
        return False

    return True


def keep_ortho_together(groups, clg, to_keep_together):

    if len(to_keep_together[clg]) == 2:
        k1, k2 = to_keep_together[clg]
        gr1, gr2 = groups

        #if genes to keep together in k1 are found sep. into gr1 and gr2 -> do not retain the topology
        if k1.intersection(gr1) and k1.intersection(gr2):
            return False

        #if genes to keep together in k2 are found sep. into gr1 and gr2 -> do not retain the topology
        if k2.intersection(gr1) and k2.intersection(gr2):
            return False

    if len(to_keep_together[clg]) == 3:
        k1, k2, k3 = to_keep_together[clg]
        gr1, gr2 = groups

        #if genes to keep together in k1 are found sep. into gr1 and gr2 -> do not retain the topology
        if k1.intersection(gr1) and k1.intersection(gr2):
            return False

        #if genes to keep together in k2 are found sep. into gr1 and gr2 -> do not retain the topology
        if k2.intersection(gr1) and k2.intersection(gr2):
            return False

        #if genes to keep together in k3 are found sep. into gr1 and gr2 -> do not retain the topology
        if k3.intersection(gr1) and k3.intersection(gr2):
            return False


    if len(to_keep_together[clg]) == 4:
        k1, k2, k3, k4 = to_keep_together[clg]
        gr1, gr2 = groups

        #if genes to keep together in k1 are found sep. into gr1 and gr2 -> do not retain the topology
        if k1.intersection(gr1) and k1.intersection(gr2):
            return False

        #if genes to keep together in k2 are found sep. into gr1 and gr2 -> do not retain the topology
        if k2.intersection(gr1) and k2.intersection(gr2):
            return False

        #if genes to keep together in k3 are found sep. into gr1 and gr2 -> do not retain the topology
        if k3.intersection(gr1) and k3.intersection(gr2):
            return False

        #if genes to keep together in k4 are found sep. into gr1 and gr2 -> do not retain the topology
        if k4.intersection(gr1) and k4.intersection(gr2):
            return False

    return True

def load_fams(clg_file):

    og_cl_prev = ''
    tmp = []
    ok = {}
    oktmp = False

    res = {}

    with open(clg_file,'r') as infile:

        for line in infile:
            line = line.strip().split('\t')
            og = line[0]
            clg = line[1]

            dup, sp, gene = line[3:6]
            chrom = line[-2]

            og_cl = og + '-' + clg

            if '/' in dup: continue

            if og_cl_prev and og_cl_prev != og_cl:
                oktmp = check_fam(tmp)
                clgprev = og_cl_prev.split('-')[-1]
                if oktmp:
                    ok[clgprev] = ok.get(clgprev, 0) + 1
                    res[clgprev] = res.get(clgprev, {})
                    res[clgprev][og_cl_prev] = tmp

                oktmp = False
                tmp = []

            tmp.append([clg, dup, sp, gene, chrom])
            og_cl_prev = og_cl

    oktmp = check_fam(tmp)
    clgprev = og_cl_prev.split('+')
    if oktmp:
        ok[clgprev] = ok.get(clgprev, 0) + 1

    return res, ok

def make_lore_ctree(family):

    cyclo = {i[2]+'_'+i[3] for i in family if i[2] in ['Petmar', 'Parata']}
    vert = {(i[2]+'_'+i[3], i[1]) for i in family if i[2] not in ['Petmar', 'Parata']}

    vert_1 = [i[0] for i in vert if i[1].replace('alpha', '').replace('beta', '')=='1']#1R_1
    vert_2 = [i[0] for i in vert if i[1].replace('alpha', '').replace('beta', '')=='2']#1R_2

    outgr = 'outgr'

    tree = Tree()
    tree.add_child(name=outgr)

    next_node = tree.add_child(name="anc")
    gr1 = next_node.add_child(name="gr1")
    for i in cyclo:
        gr1.add_child(name=i)

    gr2 = next_node.add_child(name="gr2")
    vert_1_node = gr2.add_child(name="vert1")
    for i in vert_1:
        vert_1_node.add_child(name=i)
    
    vert_2_node = gr2.add_child(name="vert2")
    for i in vert_2:
        vert_2_node.add_child(name=i)


    tree.prune(tree.get_leaves())

    # print(tree)

    return tree


def make_tree(group1, group2):
    tree = Tree()
    tree.add_child(name='outgr')

    next_node = tree.add_child(name="anc")
    gr1 = next_node.add_child(name="gr1")
    gr1_vert = gr1.add_child(name="gr1_vert")

    for i in group1:
        if 'Parata' not in i and 'Petmar' not in i:
            gr1_vert.add_child(name=i)

    cyclo1 = {i for i in group1 if 'Parata' in i or 'Petmar' in i}
    if cyclo1:
        gr1_cyclo = gr1.add_child(name="gr1_cyclo")
        for i in cyclo1:
            gr1_cyclo.add_child(name=i)

    gr2 = next_node.add_child(name="gr2")
    gr2_vert = gr2.add_child(name="gr2_vert")
    for i in group2:
        if 'Parata' not in i and 'Petmar' not in i:
            gr2_vert.add_child(name=i)


    cyclo2 = {i for i in group2 if 'Parata' in i or 'Petmar' in i}
    if cyclo2:
        gr2_cyclo = gr2.add_child(name="gr2_cyclo")
        for i in cyclo2:
            gr2_cyclo.add_child(name=i)

    tree.prune(tree.get_leaves())
    return tree

def make_all_aore_ctree(family):

    trees = []

    petmar = {i[2]+'_'+i[-1] for i in family if i[2] in ['Petmar'] and 'unloc' not in i[-1] and 'scaf' not in i[-1]}
    parata = {i[2]+'_'+i[-1] for i in family if i[2] in ['Parata'] and 'scaf' not in i[-1]}

    clg = family[0][0]

    if len(petmar) > 6 or len(parata) > 6:
        return []

    vert = {(i[2]+'_'+i[3], i[1]) for i in family if i[2] not in ['Petmar', 'Parata']}

    vert_1 = [i[0] for i in vert if i[1].replace('alpha', '').replace('beta', '')=='1']#1R_1
    vert_2 = [i[0] for i in vert if i[1].replace('alpha', '').replace('beta', '')=='2']#1R_2

    #all possible ways to group cyclo chr with Vertebrates 1R_1 and 1R_2 (max 3 chr in each group)

    petmar_1 = [] #all_combination of 1, 2 and 3 chr
    for n in range(min(3, len(petmar))):
        petmar_1 += list(itertools.combinations(petmar, n + 1))

    for comb in petmar_1:
        petmar_2 = {i for i in petmar if i not in comb} #all not in comb

        parata_1 = [] #all_combination of 1, 2 and 3 chr
        for n in range(min(3, len(parata))):
            parata_1 += list(itertools.combinations(parata, n + 1))

        for comb2 in parata_1:
            parata_2 = {i for i in parata if i not in comb2}

            tmp_petmar_1 = [i[2]+'_'+i[3] for i in family if i[2]+'_'+i[-1] in comb]
            tmp_petmar_2 = [i[2]+'_'+i[3] for i in family if i[2]+'_'+i[-1] in petmar_2]
            tmp_parata_1 = [i[2]+'_'+i[3] for i in family if i[2]+'_'+i[-1] in comb2]
            tmp_parata_2 = [i[2]+'_'+i[3] for i in family if i[2]+'_'+i[-1] in parata_2]

            if keep_ortho_together((set(list(comb) + list(comb2)), set(list(petmar_2) + list(parata_2))), clg, KEEP_TOGETHER):
                tree = make_tree(vert_1 + tmp_petmar_1 + tmp_parata_1, vert_2 + tmp_parata_2 + tmp_petmar_2)
                trees.append(tree.copy())

                tree = make_tree(vert_2 + tmp_petmar_1 + tmp_parata_1, vert_1 + tmp_parata_2 + tmp_petmar_2)
                trees.append(tree.copy())

            if keep_ortho_together((set(list(comb) + list(parata_2)), set(list(petmar_2) + list(comb2))), clg, KEEP_TOGETHER):
                tree = make_tree(vert_1 + tmp_petmar_2 + tmp_parata_1, vert_2 + tmp_parata_2 + tmp_petmar_1)
                trees.append(tree.copy())

                tree = make_tree(vert_1 + tmp_petmar_1 + tmp_parata_2, vert_2 + tmp_parata_1 + tmp_petmar_2)
                trees.append(tree.copy())

    return trees

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', help='CLG file',
                        required=True)

    parser.add_argument('-oa', '--out_aore', required=False, default='1r_aore_ctrees/')

    parser.add_argument('-ol', '--out_lore', required=False, default='1r_lore_ctrees/')

    args = vars(parser.parse_args())

    fams, counts = load_fams(args['input'])

    os.makedirs(args['out_lore'], exist_ok=True)

    os.makedirs(args['out_aore'], exist_ok=True)

    for clg in fams:
        # if clg != 'CLGD':
        #     continue
        for fam in fams[clg]:
            lore = make_lore_ctree(fams[clg][fam])
            all_aore_ctree = make_all_aore_ctree(fams[clg][fam])

            # print(len(all_aore_ctree))

            if len(all_aore_ctree) == 0 or len(all_aore_ctree) > 10:

                # print('*****')

                counts[clg] = counts[clg] - 1

            else:
                lore_outfile = args['out_lore'] + '/' + fam + '.nh'
                lore.write(outfile=lore_outfile, format=9)
                for i, aore in enumerate(all_aore_ctree):
                    # print('AORe', aore)
                    aore_outfile = args['out_aore'] + '/' + fam + '-' + str(i) + '.nh'
                    aore.write(outfile=aore_outfile, format=9)

    # print(fams)

    # print(counts)
