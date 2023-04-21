import argparse
import sys
import glob
import numpy as np

from collections import defaultdict

import pickle

def concat_seq(seq_file, d, summary):
    seen = set()
    tmp = {}
    s = {}
    with open(seq_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                g, name = line.split('=')
                sp = g.split('_')[0][1:]

                if sp == 'Braflo':
                    cluster = ''
                elif '_' in name:
                    cluster = name.split('_')[-1]

                elif 'HOX' in name.upper():
                    cluster = name[3].lower()
                    print(name, cluster)

                else:
                    cluster = name[-2].lower()

                tmp[sp+'_'+cluster.strip()] = '' #to do select less gaps

                s[sp+'_'+cluster.strip()] = line.strip()
            else:
                tmp[sp+'_'+cluster.strip()] += line.strip()


    #add gaps for missing genes
    for key in tmp:
        lg = len(tmp[key])
        break

    for key in d:
        if key in tmp:
            d[key] += tmp[key]
            summary[key].append(s[key])
        else:
            d[key] += '-' * lg
            summary[key].append('None')
    return


def write_fasta(output_file, dseq):
    with open(output_file, 'w') as out:
        for key in dseq:
            out.write('>'+key+'\n')
            out.write(dseq[key]+'\n')
    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-o', '--output', required=True)

    parser.add_argument('-i', '--input', required=True, nargs='+')

    args = vars(parser.parse_args())

    dconcat = {}
    dsummary = defaultdict(list)
    clusters = ['a', 'b', 'c', 'd']
    vert = ['Galgal', 'Homsap', 'Musmus', 'Lepocu']
    for v in vert:
        for c in clusters:
            dconcat[v+'_'+c] = ''

    clusters = ['I', 'II', 'III', 'IV', 'V', 'VI']
    for c in clusters:
        dconcat['Parata_'+c] = ''

    clusters = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'dzeta']
    for c in clusters:
        dconcat['Petmar_'+c] = ''

    dconcat['Braflo_'] = ''

    for genefam in args['input']:
        concat_seq(genefam, dconcat, dsummary)


    write_fasta(args['output'], dconcat)

    with open('clusters_concat_summay.pkl', 'wb') as out:
        pickle.dump(dsummary, out)