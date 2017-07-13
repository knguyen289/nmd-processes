print '\tImporting Packages'
import sys
sys.stdout.write('\t. ')
import re
sys.stdout.write(' . ')
import glob
sys.stdout.write(' . ')
import pandas as pd
sys.stdout.write(' . ')
import numpy as np
sys.stdout.write(' . ')
from optparse import OptionParser
sys.stdout.write(' . ')
from kn_tools.analysis_tools import fetch_coords,seq_index,gene_index
sys.stdout.write(' .\n')

sys.stdout.write('\tImporting Packages Complete\n')

def main():
	parser = OptionParser()
	
	parser.add_option("--rna", dest="rna", default=None,
		help="RNA Name for bed retrieval")
	parser.add_option("--id", dest="id",  default=None, help="Directory ID")

	(options, args) = parser.parse_args()

	if options.rna == None:
		print "Error: need --rna"
		sys.exit(1)

	elif options.id == None:
		print "Error: need --id"

	else:
		rna = str(options.rna)
		dir_num = str(options.id)

	real_info = []

	csv_name = 'ss_' + rna + '-' + dir_num + '.csv'

	sys.stdout.write('\tCommand Line parameters used: ')
	print rna,dir_num
	sys.stdout.write('\tCSV Output to: ')
	print csv_name

	# Gets info from existant strands from bedinfo
	exist_dir = glob.glob(rna + 'bedinfo' + dir_num + '/*.txt')
	exist_inds = [int(item.split('_')[-1].split('.')[0]) for item in exist_dir]
	# Make the dataframe from info from the bed files
	to_convert = []
	for fname in glob.glob(rna + 'made_beds' + dir_num + '/*.txt'):
	    chrom = ''
	    ftemp = open(fname)
	    seq = ''
	    starts = []
	    ends = []
	    nodes = []
	    strand = ''
	    for line in ftemp:
	        line = line.rstrip()
	        if len(line) > 0:
	            if line[0] == '>':
	                info = line.split(':')
	                if len(chrom) == 0:
	                    chrom = info[0].split(';')[-1]
	                starts.append(int(info[1]))
	                ends.append(int(info[2]))
	                if len(strand) == 0:
	                    strand = info[3]
	            if line[0] != '>':
	                seq += line
	    strand_id = fname.split('.')[0].split('_')[-1]
	    nodes = sorted(starts + ends)
	    lens = [ends[i] - starts[i] for i in range(len(starts))]
	    print lens, sum(lens)
	    seqNodes = [seq_index(item,nodes,strand) for item in nodes]
	    to_append = [int(strand_id),'',chrom,strand,len(starts),','.join(map(str,starts))+',',','.join(map(str,ends))+',',rna,','.join(map(str,nodes)),','.join(map(str,seqNodes)),seq]

	    to_convert.append(to_append)

	info_df = pd.DataFrame(to_convert,columns=['ID','name','chrom','strand','exonCount','exonStarts','exonEnds','name2','Nodes','seqNodes','Seq'])
	info_df = info_df.sort_values('ID')
	info_df = info_df.set_index('ID')
	info_df.insert(6,'txStart','')
	info_df.insert(7,'txEnd','')
	info_df.insert(8,'cdsStart','')
	info_df.insert(9,'cdsEnd','')
	info_df.insert(len(info_df.columns),'regexStart','')
	info_df.insert(len(info_df.columns),'regexStop','')


	#Get the starts and ends for all of the isoforms
	to_drop = []
	for index,row in info_df.iterrows():
	    if int(index) in exist_inds:
	        finfo = open(rna + 'bedinfo' + dir_num + '/' + rna + '_info_' + str(index) + '.txt')
	        exists_name = finfo.read().split('\t')[0]
	        info_df.loc[index,'name'] = exists_name

	    if 'N' in row.get_value('Seq'):
	        to_drop.append(index)
	    else:
	        s,e = fetch_coords(row.get_value('Seq'))
	        nodes = map(int,row.get_value('Nodes').split(','))
	        strand = row.get_value('strand')
	        se_list = sorted([gene_index(item,nodes,strand) for item in [s,e]])
	        nodes = row.get_value('Nodes').split(',')
	        strand = row.get_value('strand')
	        info_df.loc[index,'regexStart'] = s
	        info_df.loc[index,'regexStop'] = e
	        info_df.loc[index,'txStart'] = nodes[0]
	        info_df.loc[index,'txEnd'] = nodes[-1]
	        info_df.loc[index,'cdsStart'] = se_list[0]
	        info_df.loc[index,'cdsEnd'] = se_list[1]
	for ind in to_drop:
	    info_df.drop(ind,inplace=True)


	info_df.to_csv(csv_name,sep='\t')

if __name__ == '__main__':
	main()
