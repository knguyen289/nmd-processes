import sys
import pandas as pd
import glob
from kn_tools.basic_tools import text_to_df

data = glob.glob('/ufrc/ewang/nguyenk/nmd_reg/data_dir/*')
candidates = ''
strange = ''


for dir_name in data:
	go = True
	sys.stdout.write(dir_name + '\n')
	try:
		rna = '-'.join(dir_name.split('/')[-1].split('_')[-3].split('-')[0:-1])
		dir_ID = dir_name.split('/')[-1].split('_')[-3].split('-')[-1]

		flag_file = dir_name + '/' + rna + '-' + dir_ID + '_analysis/flag_' + rna + '-' + dir_ID + '.csv'

		flag_df = text_to_df(flag_file,index='ID')
	
	except:
		sys.stdout.write('\tEXCEPT ' + dir_name + '\n\n')
		strange += rna + '\t' + dir_name + '\n'
		go = False

	if go:
		sys.stdout.write('\tWORKED ' + dir_name + '\n')
		has_3utr = False

		utr3_inds = list(set(flag_df['exonsFromEnd']))
		sys.stdout.write('\t\t' + ' '.join(utr3_inds) + '\n')
		for i in utr3_inds:
			sys.stdout.write('\t\tUTR Test Loop\n')
			if int(i) != 0:
				has_3utr = True
				sys.stdout.write('\t\tUTR Test Loop Before Break\n')
				break

		sys.stdout.write('\t\tHas 3\' UTR? ' + str(has_3utr) + '\n')
		nmd_diverse = False

		if has_3utr:
			sys.stdout.write('\t\tNMD Diverse Test Loop\n')
			nmd_inds = sorted(list(set(flag_df['basesFromJunct'])))
			under = False
			over = False
			for i in nmd_inds:
				try:
					if int(i) < 50 and int(i) != -1:
						under = True
					if int(i) >= 50:
						over = True
					if under and over:
						nmd_diverse = True
						break
				except:
					strange += rna + '\t' + dir_name + '\n'

		sys.stdout.write('\t\tNMD_Diverse? ' + str(nmd_diverse) + '\n')
		if nmd_diverse:
			candidates += rna + '\n'
		sys.stdout.write('\n')

f = open('/ufrc/ewang/nguyenk/nmd_reg/candidates.txt','w')
f.write(candidates)

f_strange = open('/ufrc/ewang/nguyenk/nmd_reg/strange_rna.txt','w')
f_strange.write(strange)
