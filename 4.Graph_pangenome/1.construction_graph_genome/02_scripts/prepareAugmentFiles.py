import vcf
import argparse
from pyfaidx import Fasta
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq

parser = argparse.ArgumentParser(description='Make fasta for each variant to align/augment.')
parser.add_argument('-v', help='the input VCF file.', required=True)
parser.add_argument('-r', help='the reference FASTA file.', required=True)
parser.add_argument('-s', help='the output FASTA file with SV sequence to align/augment', required=True)
parser.add_argument('-f', default=50000, type=int,
                    help='the flank size. Default 50000.')
args = parser.parse_args()

# get chromosome length
ref = Fasta(args.r)

# read vcf
vcfi = open(args.v, 'r')
vcf_reader = vcf.Reader(vcfi)
fa_outf = open(args.s, 'w')
tail_buff = 1000 # tail buffer: no sequence extracted from a buffer at the chunk tails to ensure they stay untouched
for record in vcf_reader:
    chr_len = len(ref[record.CHROM])
    # retrieve alt allele with flanks
    # left flank sequence
    fl1_e = record.POS - 1
    if fl1_e < tail_buff:
        l1_s = tail_buff / 2
    else:
        fl1_s = fl1_e - args.f
        fl1_s = max(0, fl1_s)  + tail_buff
    fl1_seq = ref[record.CHROM][fl1_s:fl1_e]
    fl1_seq = fl1_seq.seq
    # Get flank 2 sequence
    fl2_s = record.POS + len(record.REF) - 1
    if fl2_s > chr_len - tail_buff:
        fl2_e = (chr_len + fl2_s)/2
    else:
        fl2_e = fl2_s + args.f
        fl2_e = min(fl2_e, len(ref[record.CHROM])) - tail_buff
    fl2_seq = ref[record.CHROM][int(fl2_s):int(fl2_e)]
    fl2_seq = fl2_seq.seq
    # Fasta record
    oseq = fl1_seq + str(record.ALT[0]) + fl2_seq
    svid = '{}_{}_{}_{}'.format(record.CHROM, int(fl1_s), int(fl2_e), record.ID)
    orec = SeqRecord(MutableSeq(oseq.upper()), id=svid,
                     description='')
    SeqIO.write(orec, fa_outf, "fasta")
fa_outf.close()
vcfi.close()
