import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from pyfaidx import Fasta
import vcf
import subprocess
import os
import shutil


parser = argparse.ArgumentParser(description='Progressively align and augment'
                                 ' a chunked graph with sequences '
                                 'overlapping chunk.')
parser.add_argument('-r', help='the reference fasta', required=True)
parser.add_argument('-f', help='the input fasta file with all sequences',
                    required=True)
parser.add_argument('-c', help='the VCF with variants to vg construct',
                    required=True)
parser.add_argument('-o', help='the output augmented .vg graph', required=True)
parser.add_argument('-d', help='debug mode: keep temporary files',
                    action='store_true')
parser.add_argument('--noaugment', action='store_true',
                    help='construcion without augmentation (e.g. for'
                    ' chunk that take too long augmenting')
args = parser.parse_args()

vg = 'vg'

# temporary files
temp_fa = args.o + '_temp.fa'
temp_vcf = args.o + '_temp.vcf'
temp_vg = args.o + '_temp.vg'
temp_vg2 = args.o + '_temp2.vg'
temp_xg = args.o + '_temp.xg'
temp_gcsa = args.o + '_temp.gcsa'
temp_fastq = args.o + '_temp.fq'
temp_gam = args.o + '_temp.gam'
# log file for all the vg commands
log_err = open(args.o + '.log', 'w')

# chunk region from the output .vg filename
# (e.g. mpmap_augment_chunks/augmented/chunk_chr16_3839328_4394815.vg)
fn = args.o.rstrip('.vg').split('/')[-1]
fn = fn.split('_')
chunk_seqn = fn[1]
chunk_s = int(fn[2])
chunk_e = int(fn[3])
print('Chunk {} from {} to {}.'.format(chunk_seqn, chunk_s, chunk_e))

# extract reference sequence
ref = Fasta(args.r)
ref_seq = ref[chunk_seqn][chunk_s:chunk_e]
ref_farec = SeqRecord(MutableSeq(ref_seq.seq.upper()),
                      id=chunk_seqn,
                      description='')
SeqIO.write(ref_farec, temp_fa, "fasta")

# extract variants to add using vg construct
vcfi = open(args.c, 'r')
vcf_reader = vcf.Reader(vcfi)
# write the VCF for this chunk in a temp file
vcfo = open(temp_vcf, 'w')
vcf_writer = vcf.Writer(vcfo, vcf_reader)
nb_var_to_construct = 0
for record in vcf_reader:
    if record.CHROM == chunk_seqn and record.POS < chunk_e \
       and record.POS > chunk_s:
        record.POS = record.POS - chunk_s
        vcf_writer.write_record(record)
        nb_var_to_construct += 1
vcfi.close()
vcfo.close()

# vg construct
print('Make graph with vg construct adding ' +
      str(nb_var_to_construct) + ' variants...')
cons_cmd = [vg, 'construct', '-t', '16', '-r', temp_fa, '-a', '-f', '-S']
if nb_var_to_construct > 0:
    # add variants from the VCF
    cons_cmd += ['-v', temp_vcf]
with open(temp_vg2, 'w') as outf:
    subprocess.run(cons_cmd, stdout=outf, stderr=log_err, check=True)

# read fasta file with SV-containing sequence to align and augment
seqs = {}
sv_sizes = {}
for record in SeqIO.parse(args.f, "fasta"):
    # keep if region overlaps the chunk region
    # sequence named like chr16_0_1970434_sv0
    seqn = record.id.split('_')
    seq_s = int(seqn[1])
    seq_e = int(seqn[2])
    if seqn[0] == chunk_seqn and seq_s < chunk_e and seq_e > chunk_s:
        seqs[record.id] = record
        sv_sizes[record.id] = int(seqn[4])

# get number of nodes and edges in graph
stats_cmd = [vg, 'stats', '-z', temp_vg2]
stats_o = subprocess.check_output(stats_cmd, stderr=log_err)
stats_o = stats_o.decode('ascii')
stats_o = stats_o.split('\n')
nb_nodes = int(stats_o[0].split('\t')[1])
nb_edges = int(stats_o[1].split('\t')[1])
print('Graph: ' + str(nb_nodes) + ' nodes and ' +
      str(nb_edges) + ' edges.')

# if no overlapping variant, stop there
if len(seqs) == 0:
    print('No variants to augment in this chunk.')
else:
    print(str(len(seqs)) + ' overlapping variants to augment.')
    # for clusters with lots of variants, we need to prune
    # more agressively (for now)
    prune_max_degree = 32
    if(len(seqs) > 20):
        prune_max_degree = 4
    # order seq names by increasing sizes
    seqn_sorted = sorted(sv_sizes.keys(), key=lambda k: sv_sizes[k])
    # for each sequence: align and augment graph
    for seqn in seqn_sorted:
        print('Add ' + seqn + '...')
        # XG index
        idx_cmd = [vg, 'index', '-t', '1', '-x', temp_xg, temp_vg2]
        subprocess.check_output(idx_cmd, stderr=log_err)
        # prepare FASTQ input
        with open(temp_fastq, 'w') as outf:
            outf.write('@' + seqs[seqn].name + '\n' +
                       str(seqs[seqn].seq) + '\n+\n' +
                       '~' * len(seqs[seqn].seq) + '\n')
        # prune vg file
        prune_cmd = [vg, 'prune', '-M', str(prune_max_degree), temp_vg2]
        with open(temp_vg, 'w') as outf:
            subprocess.run(prune_cmd, stdout=outf, stderr=log_err, check=True)
        # GCSA index
        idx_cmd = [vg, 'index', '-t', '1', '-g', temp_gcsa,
                   '-b', 'temp', temp_vg]
        subprocess.check_output(idx_cmd, stderr=log_err)
        # map with mpmap in single-path mode
        mpmap_cmd = [vg, 'mpmap','-B','-F','GAM', '-S', '-t', '16', '-l','long','-f',
                     temp_fastq, '-x', temp_xg, '-g', temp_gcsa]
        #print(str(mpmap_cmd))
        with open(temp_gam, 'w') as outf:
            subprocess.run(mpmap_cmd, stdout=outf, stderr=log_err, check=True)
        # augment
        augment_cmd = [vg, 'augment', '-i', '-t', '16', temp_vg2, temp_gam]
        with open(temp_vg, 'w') as outf:
            subprocess.run(augment_cmd, stdout=outf, stderr=log_err,
                           check=True)
        # overwrite graph
        shutil.copy(temp_vg, temp_vg2)
        # get number of nodes and edges in graph
        stats_cmd = [vg, 'stats', '-z', temp_vg2]
        stats_o = subprocess.check_output(stats_cmd, stderr=log_err)
        stats_o = stats_o.decode('ascii')
        stats_o = stats_o.split('\n')
        new_nb_nodes = int(stats_o[0].split('\t')[1])
        new_nb_edges = int(stats_o[1].split('\t')[1])
        print('Graph: ' + str(new_nb_nodes) + ' nodes and ' +
              str(new_nb_edges) + ' edges.')
        # check that the graph got bigger
        if(new_nb_edges == nb_edges and new_nb_nodes == nb_nodes):
            print('Warning: same number of nodes and edges in augmented graph')
        # update current graph size
        nb_edges = new_nb_edges
        nb_nodes = new_nb_nodes

# close log
log_err.close()
print('Check the log at ' + args.o + '.log')

shutil.copy(temp_vg2, args.o)

# remove temporary files if not in "debug" mode
if not args.d:
    if(os.path.isfile(temp_fa)):
        os.remove(temp_fa)
    if(os.path.isfile(temp_fa + '.fai')):
        os.remove(temp_fa + '.fai')
    if(os.path.isfile(temp_vcf)):
        os.remove(temp_vcf)
    if(os.path.isfile(temp_vg)):
        os.remove(temp_vg)
    if(os.path.isfile(temp_vg2)):
        os.remove(temp_vg2)
    if(os.path.isfile(temp_xg)):
        os.remove(temp_xg)
    if(os.path.isfile(temp_gcsa)):
        os.remove(temp_gcsa)
    if(os.path.isfile(temp_gcsa + '.lcp')):
        os.remove(temp_gcsa + '.lcp')
    if(os.path.isfile(temp_fastq)):
        os.remove(temp_fastq)
    if(os.path.isfile(temp_gam)):
        os.remove(temp_gam)
