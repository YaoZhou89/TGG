#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os, re,datetime, glob
st = ""
ttt = 0
button = 1
list_key =[]
print('***********************************************************')
print('Start Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')
row1_list = []
if  len(sys.argv) == 4:
	ref_genome=sys.argv[1]
else:
	print('Use Error')
	print('Example : python3 merge_2_vg.py genome merge_vcf output')
	sys.exit()
merge_sv=open(sys.argv[3],'w')
false_sv=open("false_sv","w")

def DNA_complement(sequence):
	sequence=sequence[::-1]
	trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd','TGCAtgcaYRMKyrkmBVDHbvdh')
	string = sequence.translate(trantab)
	return string

false_ins_end_list = []

with open (sys.argv[2]) as f:
	for line in f:
##retain header file
		if "#" in line:
			st += line
			continue
		if button == 1:	
			merge_sv.write(st)
			button =0
		if '#' not in line:
			row_list = line.split()
			if "SVTYPE=DEL" in line and "INS" in row_list[2]:
				merge_sv.write(line)
				continue
			if "SVTYPE=TRA"  in line:
				merge_sv.write(line)
				continue
			if "["  in line:
				merge_sv.write(line)
				continue
			if "]" in line:
				merge_sv.write(line)
				continue
			if "<INS>" in line:
				pbsv_list = row_list[13].split(":")
				cutesv_list = row_list[11].split(":")
				sniffes_list = row_list[10].split(":")
				assemblic_list = row_list[9].split(":")
				

				
				if len(pbsv_list[9])>3:
					row_list[3] = pbsv_list[8]
					row_list[4] = pbsv_list[9]
					row_list='\t'.join(row_list) + '\n'
					merge_sv.write(row_list)
					continue
				elif len(cutesv_list[9]) >3:
					row_list[3] = cutesv_list[8]
					row_list[4] = cutesv_list[9]
					row_list='\t'.join(row_list) + '\n'
					merge_sv.write(row_list)
					continue
				elif len(sniffes_list[9]) >3:
					row_list[3] = sniffes_list[8]
					row_list[4] = sniffes_list[9]
					row_list='\t'.join(row_list) + '\n'
					merge_sv.write(row_list)
					continue
				elif len(assemblic_list[9]) >3:
					row_list[3] = assemblic_list[8]
					row_list[4] = assemblic_list[9]
					row_list='\t'.join(row_list) + '\n'
					merge_sv.write(row_list)
					continue
			if "SVTYPE=INS" in line:
				if len(row_list[3]) < len(row_list[4]):
					n=os.popen('samtools faidx '+ref_genome+' '+row_list[0]+ ':'+row_list[1]+'-'+row_list[1]).readlines()[1].strip()
					if row_list[3]=="N":
						row_list[3] = n[0]
					row_list='\t'.join(row_list) + '\n' 
					merge_sv.write(row_list)
				else:
					merge_sv.write(line)
			if "<INV>" in line:
				inv_pos=re.split(';|=',row_list[7])
				inv_idx=int(inv_pos.index('END'))+1
				inv_end=inv_pos[inv_idx]
				if (int(inv_end)-int(row_list[1]))>100000:
					merge_sv.write(line)
					continue
				n=os.popen('samtools faidx '+ref_genome+' '+row_list[0]+ ':'+row_list[1]+'-'+inv_end).readlines()[1:]
				seq=''
				for i in n:
					i=i.strip()
					seq+=i
				seq=seq.strip()
				n2=DNA_complement(seq)
				row_list[3]=seq
				row_list[4]=n2
				row_list='\t'.join(row_list) + '\n'
				merge_sv.write(row_list)
				continue

			if "SVTYPE=INV" in line: 
				merge_sv.write(line)	
			if "SVTYPE=DEL" in line:
				del_pos=re.split(';|=',row_list[7])
				del_idx=int(del_pos.index('END'))+1
				del_end=del_pos[del_idx]
				n=os.popen('samtools faidx '+ref_genome+' '+row_list[0]+ ':'+row_list[1]+'-'+del_end).readlines()[1:]
				seq=''
				for i in n:
					i=i.strip()
					seq+=i
				seq=seq.strip()
				if seq[3:-3] not in row_list[3]:
					false_list = line.split(":")
					for i in false_list:
						if row_list[1] in i:
							if "END" not in i:	
								false_list2 = i.split('_')
								if "/" in false_list2[-1]:
									false_ins_end = false_list2[-1][:-4]
									false_ins_end_list.append(false_ins_end)
									
								else:
									if "\n" in false_list2[-1]:
										false_ins_end = false_list2[-1].strip()
									else:
										false_ins_end = false_list2[-1]

									false_ins_end_list.append(false_ins_end)
									false_ins_end=max(false_ins_end_list)
					false_ins_end_list =[]
					seq1=''
					n1=os.popen('samtools faidx '+ref_genome+' '+row_list[0]+ ':'+row_list[1]+'-'+false_ins_end).readlines()[1:]
					seq1=''
					for i in n1:
						i=i.strip()
						seq1+=i
						seq1=seq1.strip()
					row_list[3]=seq1
					row_list[4]=seq1[0]
					row1_list='\t'.join(row_list) + '\n' 
					merge_sv.write(row1_list)
					seq1=''
					ttt = 1
				if ttt==0:
					
					row_list[3]=seq
					row_list[4]=seq[0]
					key = row_list[2]
					row_list='\t'.join(row_list) + '\n'
					if key not in  list_key:
						merge_sv.write(row_list)
						list_key.append(key)
				if ttt ==1 :
					ttt =0
					continue
			if "SVTYPE=DUP" in line: 
				dup_pos=re.split(';|=',row_list[7])
				dup_idx=(dup_pos.index('END'))+1
				dup_end=dup_pos[dup_idx]
				seq=''
				n=os.popen('samtools faidx '+ref_genome+' '+row_list[0]+ ':'+row_list[1]+'-'+row_list[1]).readlines()[1].strip()
				n1=os.popen('samtools faidx '+ref_genome+' '+row_list[0]+ ':'+row_list[1]+'-'+dup_end).readlines()[1:]
				seq=''
				for i in n1:
					i=i.strip()
					seq+=i
				seq=seq.strip()
				row_list[4]=seq
					
				
				row_list='\t'.join(row_list) + '\n' 
				merge_sv.write(row_list)

