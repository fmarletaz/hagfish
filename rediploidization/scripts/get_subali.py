import argparse
import sys
from ete3 import Tree


def write_subali(alifile, genes, output_file):
    ok = False
    with open(alifile, 'r') as infile, open(output_file, 'w') as out:
        for line in infile:
            if line[0] == '>':
                if line[1:-1] in genes:
                    out.write(line)
                    ok = True
                elif line[1:-1].replace('.1', '') in genes: #Dirty fix for some xentro genes...
                    out.write(line.replace('.1', ''))
                    ok = True
                else:
                    ok = False

            else:
                if ok:
                    out.write(line)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t', '--tree', required=True)

    parser.add_argument('-a', '--ali', required=True)

    parser.add_argument('-o', '--output', required=True)

    args = vars(parser.parse_args())

    t = Tree(args['tree'])

    genes = {i.name for i in t.get_leaves()}

    write_subali(args['ali'], genes, args['output'])

