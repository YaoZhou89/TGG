## orignal code: https://github.com/vgteam/giraffe-sv-paper/blob/master/scripts/sv/merge-svs/
import argparse
from pyfaidx import Fasta


parser = argparse.ArgumentParser(description='Create lifted VCF.')
parser.add_argument('-v', help='the VCF file', required=True)
parser.add_argument('-l', help='the lifted over BED file', required=True)
parser.add_argument('-r', help='the reference fasta', required=True)
parser.add_argument('-o', help='the output VCF file', required=True)
args = parser.parse_args()

# Open reference file
fa = Fasta(args.r)

# Read lifted over BED file
lo = {}  # id -> coord
lof = open(args.l, 'r')
for line in lof:
    line = line.rstrip().split('\t')
    lo[line[3]] = line[:3]
lof.close()

# Parse VCF
outf = open(args.o, 'w')
filt_size = filt_unl = filt_base = end_updated = filt_seqn = passed = 0
contig_written = False
for line in open(args.v, 'r'):
    # write headers
    if line[0] == '#':
        line = line.replace('SL4.0ch0', '')
        line = line.replace('SL4.0ch', '')
        if 'contig=' in line:
            if contig_written:
                continue
            # write new contig lengths
            ctemp = '##contig=<ID={},length={}>\n'
            line = ''
            for chrn in fa.keys():
                line += ctemp.format(chrn, len(fa[chrn]))
            contig_written = True
        outf.write(line)
        continue
    # otherwise parse line
    line = line.rstrip().split('\t')
    vid = line[2]
    info = {}
    for infof in line[7].split(';'):
        infof = infof.split('=')
        if len(infof) == 1:
            infof.append('')
        info[infof[0]] = infof[1]
    symb = '<' in line[4]
    #print (symb)
    # if variant not in lifted over BED, skip
    if vid not in lo:
        filt_unl += 1
        continue
    los = int(lo[vid][1])
    loe = int(lo[vid][2])
    # if lifted to different chromosome, skip
    chrs = line[0].replace('SL4.0ch0', '')
    chrs = chrs.replace('SL4.0ch', '')
    if lo[vid][0] != chrs:
        filt_seqn += 1
        continue
    # if different size, skip
    sv_size = len(line[3])
    if symb:
        # if symbolic variant, use END
        if 'END' not in info:
            print('Error: Symbolic variant and no END information.')
            exit
        sv_size = int(info['END']) - int(line[1])
    if sv_size != (loe - los):
        filt_size += 1
        continue
    newseq = ''
    if not symb:
        # if first base different, skip
        newseq = fa[lo[vid][0]][los:loe]
        newseq = newseq.seq.upper()
        line[3] = line[3].upper()
        if line[3][0] != newseq[0]:
            filt_base += 1
            continue
    # update and write
    line[1] = str(los + 1)
    line[0] = chrs
    if not symb:
        line[3] = newseq
    # update END field if present
    if 'END' in info:
        end_updated += 1
        info['END'] = str(loe + 1)
        outinfo = []
        for field in info:
            if info[field] != '':
                info[field] = '=' + info[field]
            outinfo.append('{}{}'.format(field, info[field]))
        line[7] = ';'.join(outinfo)
    passed += 1
    outf.write('\t'.join(line) + '\n')

print('Filtered: unlifted:{}\tseqn:{}\tsize:{}\tbase:{}'.format(filt_unl,
                                                                filt_seqn,
                                                                filt_size,
                                                                filt_base))
print('Passed: {}'.format(passed))
print('End updated: {}'.format(end_updated))
outf.close()
