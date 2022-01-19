import argparse
import gzip


parser = argparse.ArgumentParser(description='Remove samples and most INFO.')
parser.add_argument('-v', help='the VCF file', required=True)
parser.add_argument('-p', help='the prefix for the IDs', default='')
parser.add_argument('-o', help='the output VCF file', required=True)
args = parser.parse_args()

# open files
if args.v[-3:] == '.gz':
    inf = gzip.open(args.v, 'rt')
else:
    inf = open(args.v, 'r')

if args.o[-3:] == '.gz':
    outf = gzip.open(args.o, 'wt')
else:
    outf = open(args.o, 'w')

# read input
info_added = False
for line in inf:
    line = line.rstrip()
    # write headers
    if line[:2] == '##':
        if '##SAMPLE=' in line or '##FORMAT=' in line \
           or '##INFO=' in line or '##ALT=' in line \
           or '##FILTER=' in line or 'fileDate' in line:
            continue
        outf.write(line + '\n')
        continue
    # INFO field for the variant ID
    if not info_added:
        outf.write('##INFO=<ID=ID,Number=A,Type=String,' +
                   'Description="IDs of the original variants">\n')
        info_added = True
    # otherwise parse line and keep first columns
    line = line.split('\t')
    if line[0][0] != '#':
        # overwrite INFOs with ID info
        vid = line[2]
        if vid == '.':
            vid = '{}_{}_{}'.format(line[0], line[1], len(line[3]))
        line[7] = 'ID={}'.format(args.p + vid)
        # remove potiential QUAL and FILTER
        line[5] = '.'
        line[6] = '.'
    outf.write('\t'.join(line[:8]) + '\n')

inf.close()
outf.close()
