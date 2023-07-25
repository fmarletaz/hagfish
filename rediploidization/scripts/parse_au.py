#!/usr/bin/env python


"""
    Script to parse results of gene trees likelihood AU-test
    as written by CONSEL (Shimodaira, 2002).
"""

import os
import sys
import argparse
import pickle
import glob

def one_file_consel_4_trees(filename, alpha=0.05, item_dict={"1":"ml", "2":"lore"}, not_in_dict='aore'):

    """
    Parses a CONSEL result file for a comparison of 4 trees.

    Args:
        filename (str) : name of CONSEL file to parse
        item_dict (dict): correspondance between item in CONSEL and tree labels, one should be ml tree
    """

    tmp = []

    if  os.stat(filename).st_size != 0:

        
        i = 1

        seen = []

        with open(filename, 'r') as infile:

            for line in infile:

                if line[0] == '#' and line[0:3] != '# r':

                    res = line.strip().split()

                    #If no error occurred each result line in consel should contain more than 4 col
                    if len(res) > 4:

                        item = res[2]
                        if item in item_dict:
                            category = item_dict[item]
                        else:
                            category = not_in_dict

                        au_value = res[4]
                        lkdiff = res[3]

                        if float(au_value) < alpha:
                            au_test = 'reject'

                        else:
                            au_test = 'accept'

                        if category == 'ml' and float(au_value) < alpha:
                            tmp = ['convergence_pb']
                            break

                        if category != 'ml':
                            if i == 1:
                                lkref = float(lkdiff)
                                lkdiff = 0
                            else:
                                if lkref < 0:
                                    lkdiff = float(lkdiff)
                                else:
                                    lkdiff = float(lkdiff) - lkref
                            if category not in seen:
                                tmp.append((category, i, lkdiff, au_value))
                                i += 1
                                seen.append(category)

                    else:
                        tmp = ['error']
                        break
    return tmp



def summary_4trees(filenames, alpha=0.05, item_dict={"1":"ml", "2":"lore"}, parse_name=True, not_in_dict='aore'):

    """
    Parses all consel outputs in the input list and returns a summary, telling, for each tree, if
    the alt or hughes (or both) topologies can be rejected.
    Also prints statistics to stdout.

    Args:
        filenames (list of str): List of consel result files.
        alpha (float): alpha threshold for significance of the AU-test.
        item_dict (str, optional): tested tree consel label, the other is considered the reference.
        parse_name (bool, optional): parse filename (expects SCORPiOs naming pattern)

    Returns:
        dict: dictionary with the AU-tests results summary.

    """

    res_dict = {}
    all_res = {}

    for filename in filenames:

        #naming parsing is specific to SCORPiOs inputs and outputs
        #parsing below allows to retrieve the outgroup gene used as name for gene families.
        if parse_name:
            name = os.path.splitext(os.path.basename(filename))[0]
        else:
            name = filename

        au_result = one_file_consel_4_trees(filename, alpha, item_dict, not_in_dict=not_in_dict)
        res_dict[name] = au_result

    return res_dict

if __name__ == '__main__':


    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--input', type=str, help='Folder with all inputs', required=True)

    PARSER.add_argument('-o', '--output', help='Output file', default='out.pkl', required=False)


    ARGS = vars(PARSER.parse_args())


    ACCEPTED_TREES_1, ACCEPTED_TREES_2 = [], []


    INPUTS = glob.glob(ARGS['input']+'/*')
    SUMMARY_TREES = summary_4trees(INPUTS)

    print(SUMMARY_TREES)

    with open(ARGS['output'], "wb") as OUTFILE:
        pickle.dump(SUMMARY_TREES, OUTFILE)
        
