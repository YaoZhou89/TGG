import vcf
import argparse
from pyfaidx import Fasta
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.Seq import Seq


parser = argparse.ArgumentParser(description='Extract ref sequence and variants for a cluster')
parser.add_argument('-f', help='the reference genome fasta', required=True)
parser.add_argument('-v', help='the input VCF file.', required=True)
parser.add_argument('-t', help='the TSV with cluster information.', required=True)
parser.add_argument('-c', help='the cluster to extract.', required=True)
parser.add_argument('-ov', help='the output VCF file', required=True)
parser.add_argument('-of', help='the output FASTA file', required=True)
args = parser.parse_args()

## extract cluster information from tsv
tsv_in = open(args.t, 'r')
chr_name = ''
start_pos = 0
end_pos = 0
svids = []
for line in tsv_in:
    line = line.rstrip().split('\t')
    if(line[0] == args.c):
        chr_name = line[1]
        start_pos = int(line[2])
        end_pos = int(line[3])
        svids = line[4].split(',')

# retrieve reference sequence of the region
# Open reference fasta
ref = Fasta(args.f)
reg_seq = ref[chr_name][start_pos:end_pos]
reg_seq = reg_seq.seq

# read vcf
vcfi = open(args.v, 'r')
vcf_reader = vcf.Reader(vcfi)
vcf_out = []
for record in vcf_reader:
    # skip if not in variants of interest
    if(record.ID not in svids):
        continue
    # make a VCF record
    var_pos = record.POS - start_pos
    rec = [chr_name, str(var_pos), record.ID, str(record.REF), str(record.ALT[0]),
           '.', '.', '.']
    rec = '\t'.join(rec)
    vcf_out.append(rec)
vcfi.close()

# write VCF
# VCF header
vcf_h = '##fileformat=VCFv4.2\n'
vcf_h += '##contig=<ID={},length={}>\n'.format(chr_name, len(reg_seq))
vcf_h += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
with open(args.ov, 'w') as outf:
    outf.write(vcf_h + '\n'.join(vcf_out))

# write FASTA with the reference sequence
fa_rec = SeqRecord(MutableSeq(reg_seq.upper()), id=chr_name,
                   description='cl_' + args.c)
# write fasta
SeqIO.write(fa_rec, args.of, "fasta")
