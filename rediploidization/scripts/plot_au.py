import pickle
import argparse

from scipy.stats import fisher_exact
import pingouin as pg

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--input', type=str, required=True)

    # PARSER.add_argument('-l', '--lore', type=str, required=False)

    PARSER.add_argument('-o', '--output', required=False, default='out_rank1-ogra_tmp.svg')

    PARSER.add_argument('-o2', '--output2', required=False, default='out_sign-rank1_ogra_tmp.svg')

    ARGS = vars(PARSER.parse_args())

    with open(ARGS["input"], "rb") as infile:
        a = pickle.load(infile)

    data = {}
    records = []
    records2 = []
    k = 0
    counts = {}
    for key in a:
        rec = a[key]
        if len(rec) == 2:

            clg = key.split('-')[1]
            first = rec[0][0]
            records.append([key, clg[-1], first])

            if float(rec[0][3]) > 0.05 and float(rec[1][3]) < 0.05:
                records2.append([key, clg[-1], first])
                print(([key, clg[-1], first]))
                counts[clg[-1]] = counts.get(clg[-1], {})
                counts[clg[-1]][first] = counts[clg[-1]].get(first, 0) + 1


    ORDER = ['1', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q']
    rank_1 = pd.DataFrame.from_records(records, columns=['fam', 'CLG', 'hypothesis'])
    rank_1_sign = pd.DataFrame.from_records(records2, columns=['fam', 'CLG', 'hypothesis'])

    rank_1_sign.to_csv('family_sign_lore_aore.csv', index=False, sep='\t')
    # Make the barplot
    bar = sns.countplot(x="CLG", hue='hypothesis', data=rank_1, palette=['lightblue', 'pink'], order=ORDER);

    sns.despine()
    plt.tight_layout()
    plt.savefig(ARGS["output"])
    # plt.show()

    plt.close("all")


    # Make the barplot
    bar = sns.countplot(x="CLG", hue='hypothesis', data=rank_1_sign, palette=['steelblue', 'lightpink'], order=ORDER);

    sns.despine()
    plt.tight_layout()

    plt.savefig(ARGS["output2"])
    plt.close("all")

    print(counts)

    pvals = []
    for clg in ORDER:
        if clg not in counts:
            continue
        a = counts[clg].get('aore', 0)
        l = counts[clg].get('lore', 0)
        tot =  a + l
        oddsratio, pvalue = fisher_exact([[a, tot-a], [l, tot-l]])
        pvals.append(pvalue)

    reject, pvals_corr = pg.multicomp(pvals, method='fdr_bh')
    # for i, clg in enumerate(ORDER):
    #     print(clg, reject[i], pvals_corr[i])

