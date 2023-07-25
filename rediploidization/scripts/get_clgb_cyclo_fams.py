

import argparse
import os

from ete3 import Tree
import itertools

TO_TEST = {'Petmar_chrm10', 'Petmar_chrm2', 'Parata_chr4', 'Parata_chr8'}


def load_fams(clg_file, clg_ok='CLGB'):

	og_cl_prev = ''
	tmp = []
	ok = {}
	oktmp = False

	res = {}

	with open(clg_file,'r') as infile:

		for line in infile:
			line = line.strip().split('\t')
			og = line[0]
			clg = line[1][:4]
			if clg != clg_ok:
				continue

			dup, sp, gene = line[3:6]
			chrom = line[-2]

			og_cl = og + '-' + clg

			if '/' in dup: continue

			if og_cl_prev and og_cl_prev != og_cl:
				clgprev = og_cl_prev.split('-')[-1]
				ok[clgprev] = ok.get(clgprev, 0) + 1
				res[clgprev] = res.get(clgprev, {})
				res[clgprev][og_cl_prev] = tmp

				oktmp = False
				tmp = []

			tmp.append([clg, dup, sp, gene, chrom])
			og_cl_prev = og_cl

	clgprev = og_cl_prev.split('+')
	if oktmp:
		ok[clgprev] = ok.get(clgprev, 0) + 1

	return res, ok



def make_trees(family, to_test=TO_TEST):

	lamprey_genes = {i[2]+'_'+i[3] for i in family if i[2]+'_'+i[4] in to_test and i[2] == 'Petmar'}
	hagfish_genes = {i[2]+'_'+i[3]  for i in family if i[2]+'_'+i[4] in to_test and i[2] == 'Parata'}
	vert = {(i[2]+'_'+i[3], i[1]) for i in family if i[2] not in ['Petmar', 'Parata']}
	vert_1 = list({i[0] for i in vert if i[1].replace('alpha', '').replace('beta', '')=='1'})#1R_1

	if not vert_1:
		vert_1 = list({i[0] for i in vert })#1R_1

	if not vert_1:
		return False, False, False

	l_chr = {i[4] for i in family if i[2]+'_'+i[4] in to_test and i[2] == 'Petmar'}
	h_chr = {i[4] for i in family if i[2]+'_'+i[4] in to_test and i[2] == 'Parata'}

	if len(l_chr) != len(lamprey_genes) or len(hagfish_genes) != len(h_chr):
		if lamprey_genes and hagfish_genes:
			print('!***!')
			print(hagfish_genes, lamprey_genes)
			print('!***!')
		return 'NA', False, False

	if (len(l_chr) == 1 and len(h_chr) == 2) or (len(l_chr) == 1 and len(h_chr) == 2) or (len(l_chr) == 2 and len(h_chr) == 2):
		
		# LORE
		tree1 = Tree()
		og = tree1.add_child(name='outgr')
		for i in vert_1:
			og.add_child(name=i)
			break

		next_node = tree1.add_child(name="cyclo")
		hag = next_node.add_child(name="hagfish")
		for i in hagfish_genes:
			hag.add_child(name=i)

		lamp = next_node.add_child(name="lamprey")
		for i in lamprey_genes:
			lamp.add_child(name=i)

		tree1.prune(tree1.get_leaves())
		# print('lore', tree1)

		# AORE 1

		lamprey_genes = list(lamprey_genes)
		hagfish_genes = list(hagfish_genes)

		if len(lamprey_genes) == 2 and len(hagfish_genes) == 2:
			gr1 = [lamprey_genes[0], hagfish_genes[0]]
			gr2 = [lamprey_genes[1], hagfish_genes[1]]

		elif len(lamprey_genes) == 1:
			gr1 = [lamprey_genes[0], hagfish_genes[0]]
			gr2 = [hagfish_genes[1]]

		else:
			gr1 = [lamprey_genes[0], hagfish_genes[0]]
			gr2 = [lamprey_genes[1]]

		atree1 = Tree()
		out = atree1.add_child(name='outgr')
		for i in list(vert_1):
			out.add_child(name=i)
			break

		next_node = atree1.add_child(name="cyclo")
		g1 = next_node.add_child(name="gr1")
		for i in gr1:
			g1.add_child(name=i)

		g2 = next_node.add_child(name="gr2")
		for i in gr2:
			g2.add_child(name=i)

		atree1.prune(atree1.get_leaves())
		# print(atree1)


		# AORE 2

		if len(lamprey_genes) == 2 and len(hagfish_genes) == 2:
			gr1 = [lamprey_genes[1], hagfish_genes[0]]
			gr2 = [lamprey_genes[0], hagfish_genes[1]]

		elif len(lamprey_genes) == 1:
			gr1 = [hagfish_genes[0]]
			gr2 = [hagfish_genes[1], lamprey_genes[0]]

		else:
			gr1 = [lamprey_genes[0]]
			gr2 = [lamprey_genes[1], hagfish_genes[0]]

		atree2 = Tree()
		out = atree2.add_child(name='outgr')
		for i in vert_1[:1]:
			out.add_child(name=i)
			break

		next_node = atree2.add_child(name="cyclo")
		g1 = next_node.add_child(name="gr1")
		for i in gr1:
			g1.add_child(name=i)

		g2 = next_node.add_child(name="gr2")
		for i in gr2:
			g2.add_child(name=i)

		atree2.prune(atree2.get_leaves())
		# print(atree2)

		return tree1, atree1, atree2

	else:

		if lamprey_genes and hagfish_genes and vert:
			if len(lamprey_genes) > 2:
				print('***')
				# print(family, lamprey_genes, hagfish_genes, vert, vert_1)
				print([i[2]+'_'+i[4] for i in family if i[2]+'_'+i[4] in to_test and i[2] == 'Petmar'])
				print('***')

			if len(hagfish_genes) > 2:
				print('***')
				# print(family, lamprey_genes, hagfish_genes, vert, vert_1)
				print([i[2]+'_'+i[4] for i in family if i[2]+'_'+i[4] in to_test and i[2] == 'Parata'])
				print('***')
		return False, False, False


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__,
									 formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('-i', '--input', help='CLG file',
						required=True)

	parser.add_argument('-oa', '--out_aore', required=False, default='aore_clgb_ctrees/')

	parser.add_argument('-ol', '--out_lore', required=False, default='lore_clgb_ctrees/')

	args = vars(parser.parse_args())

	fams, counts = load_fams(args['input'])

	os.makedirs(args['out_lore'], exist_ok=True)

	os.makedirs(args['out_aore'], exist_ok=True)

	k = 0
	t = 0

	for clg in fams:

		for fam in fams[clg]:
			lore, aore1, aore2 = make_trees(fams[clg][fam])

			if not lore:
				continue

			if lore == 'NA':
				print(fam)
				t+= 1
				continue

			k += 1
			lore_outfile = args['out_lore'] + '/' + fam + '.nh'
			lore.write(outfile=lore_outfile, format=9)
			for i, aore in enumerate([aore1, aore2]):
				aore_outfile = args['out_aore'] + '/' + fam + '-' + str(i) + '.nh'
				aore.write(outfile=aore_outfile, format=9)

			# for aore in all_aore_ctree:
			#     print('**', aore)

			#
	# print(fams)

	print(k)

print(len(fams['CLGB']))
print(t)
#{'CLGB': 212, 'CLGA1': 157, 'CLGI': 150, 'CLGD': 81, 'CLGM': 93, 'CLGQ': 102, 'CLGL': 90, 'CLGC': 196, 'CLGJ': 101, 'CLGH': 100, 'CLGP': 56, 'CLGN': 91, 'CLGG': 89, 'CLGF': 159, 'CLGE': 132, 'CLGK': 116, 'CLGO': 109}


#{'CLGB': 102, 'CLGA1': 95, 'CLGI': 125, 'CLGD': 81, 'CLGM': 69, 'CLGQ': 89, 'CLGL': 72, 'CLGC': 133, 'CLGJ': 60, 'CLGH': 69, 'CLGP': 51, 'CLGN': 64, 'CLGG': 67, 'CLGF': 123, 'CLGE': 78, 'CLGK': 93, 'CLGO': 79}


#{'CLGB': 181, 'CLGA': 44, 'CLGI': 107, 'CLGD': 100, 'CLGM': 64, 'CLGQ': 75, 'CLGL': 71, 'CLGC': 170, 'CLGJ': 51, 'CLGH': 87, 'CLGP': 42, 'CLGN': 56, 'CLGG': 57, 'CLGF': 112, 'CLGE': 113, 'CLGK': 85, 'CLGO': 75}
#1,490 trees
#I had 5,559 for teleosts but let us see