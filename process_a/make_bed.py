import sys
from optparse import OptionParser

from kn_tools.basic_tools import text_to_df
from kn_tools.analysis_tools import go_to_bed

def main():
	parser = OptionParser()
	
	parser.add_option("--rna", dest="rna", nargs=1, default=None,
		help="RNA Name for bed retrieval")
	parser.add_option("--data", dest="data", nargs=1, default=None, help="The USCS browser .txt file")

	(options, args) = parser.parse_args()

	if options.rna == None:
		print "Error: need --rna"
		sys.exit(1)

	elif options.data == None:
		print "Error: need --data"

	else:
		rna = options.rna
		data = options.data

	ID = go_to_bed(rna,data)
	fname = 'id_' + rna + '.txt'

	f = open(fname,'w')
	f.write(str(ID))
	f.close()

if __name__ == '__main__':
	main()