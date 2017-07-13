import pandas as pd 
import sys
from optparse import OptionParser
import glob
import copy

from kn_tools.basic_tools import text_to_df
from kn_tools.analysis_tools import set_flags

def main():
	parser = OptionParser()
	
	parser.add_option("--rna", dest="rna", default=None,
		help="RNA Name for bed retrieval")
	parser.add_option("--id", dest="id",  default=None, help="Directory ID")
	parser.add_option("--dir", dest="dir", default=None,
		help="Directory of the ss csv")

	(options, args) = parser.parse_args()

	if options.rna == None:
		print "Error: need --rna"
		sys.exit(1)

	elif options.id == None:
		print "Error: need --id"

	elif options.dir == None:
		print "Error: need --dir"
	else:
		rna = str(options.rna)
		dir_num = str(options.id)
		loc = str(options.dir)
		fname = loc + '/ss_' + rna + '-' + dir_num + '.csv'
		oname = loc + '/flag_' + rna + '-' + dir_num + '.csv'

	ss_df = text_to_df(fname)

	ss_df.insert(len(ss_df.columns),'exonsFromEnd',0)
	ss_df.insert(len(ss_df.columns),'basesFromJunct','')
	ss_df.insert(len(ss_df.columns),'exists',0)
	ss_df.insert(len(ss_df.columns),'theorized_nmd',0)
	ss_df.insert(len(ss_df.columns),'erroneous_nmd',0)
	ss_df['ID'] = map(int,ss_df['ID'])
	ss_df = ss_df.set_index('ID')

	real_ids = [int(item.split('.')[1].split('_')[-1]) for item in glob.glob('./' + rna + 'bedinfo' + dir_num + '/*')]
	for num in real_ids:
		ss_df.loc[num,'exists'] = 1
	ss_df.to_excel(rna + '_pre_flag.xlsx')
	flag_df = set_flags(ss_df)

	flag_df.to_csv(oname,sep='\t')

if __name__ == '__main__':
	main()