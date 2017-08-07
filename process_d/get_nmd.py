import pandas as pd
import numpy as np
import glob
import copy
import math
import itertools
import scipy.stats as stats

from kn_tools.basic_tools import text_to_df
from kn_tools.rna_path_tools import get_lookup2,get_rna_dfs
from kn_tools.analysis_tools import mutual_inf,get_splice_sites

data = glob.glob('/ufrc/ewang/nguyenk/nmd_reg/csv_dir/*')
to_df = []
splices = []
for csv in data:
	mod_csv = csv.split('/ufrc/ewang/nguyenk/nmd_reg/csv_dir/')[-1]
	name2 = '-'.join(mod_csv.split('_nmd')[0].split('-')[:-1])
	dir_id = mod_csv.split('_nmd')[0].split('-')[-1]
	
	# Read file and convert to DataFrame
	mod_df = text_to_df(csv,index='ID')
	mod_df.insert(1,'name2',name2)

	# Get lookup2 table to get psuedoexons
	genome_df = text_to_df('human.txt')
	name2_df = get_rna_dfs(name2,genome_df)
	strand = list(set(name2_df['strand']))[0]
	mod_df.insert(2,'strand',strand)
	lu2 = get_lookup2(name2_df)

	# Convert exon names in header to int, table values to int
	exon_cols = list(mod_df.columns)[5:-1]
	for col in exon_cols:
	    mod_df.insert(len(mod_df.columns)-1,int(col),'')
	    mod_df[int(col)] = map(int,mod_df[col])
	    del mod_df[col]
	mod_df['NMD_ind'] = map(int,mod_df['NMD_ind'])
	mod_df['cdsStart'],mod_df['cdsEnd'] = map(int,mod_df['cdsStart']),map(int,mod_df['cdsEnd'])
	exon_cols = map(int,exon_cols)
	splices.append(get_splice_sites(mod_df))

	# Cut out rows with NMD_ind > 50 and NMD_ind == -1
	temp_df = mod_df[mod_df['NMD_ind'] < 50]
	nmd_short_df = temp_df[temp_df['NMD_ind'] != -1]

	# Cut out columns that are 0 mod 3 or do not exist in isoforms
	redact_df = copy.deepcopy(nmd_short_df)
	for col in exon_cols:
		uniq_mods = sorted(list(set(nmd_short_df[col])))
		if uniq_mods == [0] or uniq_mods == [-1] or uniq_mods == [-1,0]:
			del redact_df[col]

	# Cut out columns that do not include alternative exons
	alt_df = copy.deepcopy(redact_df)
	alt_ones = []
	alt_twos = []
	rest_cols = list(alt_df.columns)[5:-1]
	for col in rest_cols:
		uniq_mods = sorted(list(set(nmd_short_df[col])))
		if not -1 in uniq_mods:
			del alt_df[col]
		elif len(uniq_mods) > 2:
			del alt_df[col]
		elif uniq_mods == [-1,1]:
			alt_ones.append(col)
			temp_mods = list(alt_df[col])
			temp_exists = [0 if mod == -1 else 1 for mod in temp_mods]
			alt_df[col] = temp_exists
		elif uniq_mods == [-1,2]:
			alt_twos.append(col)
			temp_mods = list(alt_df[col])
			temp_exists = [0 if mod == -1 else 1 for mod in temp_mods]
			alt_df[col] = temp_exists
	alt = alt_ones + alt_twos

	# Get all pairs of exons, further filter by nmd index and compute mutual information and fisher
	inds_set = list(set(list(alt_df['NMD_ind'])))
	mi_pairs = list(itertools.combinations(alt,2))
	for pair in mi_pairs: 
	    # Split up calculation by nmd index
	    for nmd_ind in inds_set:
	        ind_df = alt_df[alt_df['NMD_ind'] == nmd_ind]
	    	
	    	# Identify the mod of exon a
	        l_a = ind_df[pair[0]]
	        if pair[0] in alt_ones:
	            mod_a = 1
	        else:
	            mod_a = 2

	        # Identify the mod of exon b
	        l_b = ind_df[pair[1]]
	        if pair[1] in alt_ones:
	            mod_b = 1
	        else:
	            mod_b = 2
	        
	        # Existency pairs
	        pair_test = [(l_a[i],l_b[i]) for i in range(len(l_a))]

	        # Adding info that doesn't involve calculations
	        mi_toadd = [name2,pair[0],pair[1],lu2.loc[int(pair[0]),'start'],lu2.loc[int(pair[0]),'end']] + [lu2.loc[int(pair[1]),'start'],lu2.loc[int(pair[1]),'end']] + [mod_a,mod_b]

	        # Calculations
	        mi,det = mutual_inf(pair_test)
	        truth_tab = [[det[(1,1)],det[(0,1)]],[det[(1,0)],det[(0,0)]]]
	        oddsratio,fish_p = stats.fisher_exact(truth_tab)
	        
	        # Adding in calculations
	        mi_toadd += [det[(0,0)],det[(1,1)],det[(0,1)],det[(1,0)],mi,fish_p,nmd_ind]

	        to_df.append(mi_toadd)

mi_df = pd.DataFrame(to_df,columns=['name2','pexon1','pexon2','start1','end1','start2','end2','mod1','mod2',
                                    (0,0),(1,1),(0,1),(1,0),'MI','Fisher','NMD_Ind'])
		
writer = pd.ExcelWriter('mutual_inf.xlsx')
mi_df.to_excel(writer)
mi_df.to_csv('mutual_inf.csv')
writer.save()

all_splice = pd.concat(splices)
all_splice = all_splice.reset_index()
writer2 = pd.ExcelWriter('splicy.xlsx')
all_splice.to_excel(writer2)
all_splice.to_csv('splicy.csv')
writer.save()