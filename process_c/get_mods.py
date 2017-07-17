import pandas as pd
import sys
import numpy as np
from kn_tools.basic_tools import text_to_df
from kn_tools.rna_path_tools import get_rna_dfs,get_lookup2,get_pexons

rna_list = []
f = open('/ufrc/ewang/nguyenk/nmd_reg/candidates.txt')
for line in f:
	line = line.strip()
	if len(line) > 0:
		rna_list.append(line)
dir_ID = '1'
genome_df = text_to_df('/ufrc/ewang/nguyenk/nmd_reg/human.txt')
err_str = ''
for rna in rna_list:
	print rna
	#try:
	df = get_rna_dfs(rna,genome_df)
	flag_file = '/ufrc/ewang/nguyenk/nmd_reg/data_dir/' + rna + '-' + dir_ID + '_all_data/' + rna + '-' + dir_ID + '_analysis/flag_' + rna + '-' + dir_ID + '.csv'
	flag_df = text_to_df(flag_file,index='ID')
	### OUTPUT MODS DF ###
	lu2 = get_lookup2(df)
	nodes,se_out,names,inds = get_pexons(flag_df,lu_df=lu2)
	pexons = map(int,list(lu2.index))
	#reverse sort if it is minus
	is_minus = False
	if list(flag_df['strand'])[0] == '-':
		is_minus = True
	all_mod = []

	for i in range(len(nodes)):
		i_node = nodes[i]
		count = len(i_node)
		s = se_out[i][0]
		e = se_out[i][1]
		mod_list = map(int,list(-1*np.ones(len(lu2))))

		beg = -1
		end = -1

		print i_node
		# Go through each pexon from start to end to test if
		for m in range(1,len(pexons)+1):
			s_temp = lu2.loc[m,'start']
			e_temp = lu2.loc[m,'end']
			print s_temp,s,e_temp
			#finding where cdsStart is
			if s_temp <= s and s < e_temp:
				mod_list[m-1] = abs(e_temp - s) % 3
				if m in i_node:
					beg = i_node.index(m)
				break

		for n in range(len(pexons),0,-1):
			s_temp = lu2.loc[n,'start']
			e_temp = lu2.loc[n,'end']

			#finding where cdsEnd is
			if s_temp <= e and e < e_temp:
				mod_list[n-1] = abs(e - s_temp) % 3
				if n in i_node:
					end = i_node.index(n)
				break
		for ex_i in range(beg+1,end):
			mod_list[i_node[ex_i]-1] = lu2.loc[i_node[ex_i],'mod']

		nmd_ind = int(flag_df.loc[str(inds[i]),'basesFromJunct'])
		all_mod.append([inds[i]] + se_out[i] + mod_list + [nmd_ind])
		#all_mod.append([inds[i]] + mod_list + [nmd_ind])

	mod_df = pd.DataFrame(all_mod,columns=['ID','cdsStart','cdsEnd']+ pexons +['NMD_ind'])
	#mod_df = pd.DataFrame(all_mod,columns=['ID']+ pexons +['NMD_ind'])
	sys.stdout.write('Got mod_df for ' + rna + '\n')
	mod_df.insert(1,'name',names)
	mod_df.set_index('ID',inplace=True)
	mod_df.to_csv('/ufrc/ewang/nguyenk/nmd_reg/csv_dir/' + rna + '-' + dir_ID + '_nmd_ind.csv',sep='\t')
	print 'Output mod_df for ' + rna + '\n'
	#except:
		#err_str += rna + '\n'
f2 = open('/ufrc/ewang/nguyenk/nmd_reg/err_out.txt','w')
f2.write(err_str)
f2.close()