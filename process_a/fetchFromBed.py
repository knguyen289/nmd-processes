import os, sys
from optparse import OptionParser
from string import maketrans


def fetchFromBed(bed_f, genomeDir, out_f):

    revcompl = maketrans("actguACTGU|","tgacaTGACA|")
    files = {}
    i = 0
    out = open(out_f,'w')
    for line in open(bed_f):
        chrom, start, end, name, qual, strand = line.strip().split("\t")[:6]
        if chrom not in files:
            x = open(os.path.join(genomeDir, chrom + ".fa"))
            x.readline()
            offset = x.tell()
            files[chrom] = [x, offset]
            file, offset = files[chrom]
        else:
            file, offset = files[chrom]
        try:
            start, end = map(int, [start, end])
            coords = [start, end]
            coords.sort()
            start, end = coords
            ostart = start + start / 50 + offset
            oend = end + end / 50 + offset
            file.seek(ostart)
            seq = file.read(oend - ostart).upper()
            seq = seq.replace("\n", "")
            if strand=='-':
                seq = seq.translate(revcompl)[::-1]
            out.write('>'+name+";"+":".join(map(str,\
                [chrom,start,end,strand]))+\
                "\n"+seq+"\n")
        except:
            print 'ERROR', chrom,start,end,strand
        i += 1
    out.close()


def main():
    parser = OptionParser()
    parser.add_option("--fetch", dest="fetch", nargs=3, default=None,
        help="Fetch sequences from regions defined by a bed file. Takes the arguments: (1) bed file, (2) directory of sequence files, and (3) output file")

    (options, args) = parser.parse_args()

    if options.fetch == None:
        print "Error: need --fetch"
        sys.exit(1)

    if options.fetch != None:
        bed_f = os.path.abspath(os.path.expanduser(options.fetch[0]))
        genomeDir = os.path.abspath(os.path.expanduser(options.fetch[1]))
        out_f = os.path.abspath(os.path.expanduser(options.fetch[2]))

        fetchFromBed(bed_f, genomeDir, out_f)

if __name__ == '__main__':
    main()
 
    

