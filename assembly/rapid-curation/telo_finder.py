#!/usr/bin/env python3


'''
MIT License
 
Copyright (c) 2022 Genome Research Ltd.
 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



'''

import os
import sys
import re
import argparse
import pyfastaq
import subprocess
from subprocess import check_output
from datetime import datetime
from time import gmtime, strftime
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import islice
import time


#other possibles to include - iyBuaLabo AAGAGG and/or GCTCCCTG.  Jo thinks it's the former


canonical=["TTAGG","TTGGG","TTGTGG","GAGCCTTGTTT","TCAGG","TTGCA", \
"TCTGGG","TTTGGATAGG","TTCGGG","ACTGGTGT","TTTAGGG","TGGGTC"]




def cprop(fasta):

	
	A_count=0
	C_count=0
	G_count=0
	T_count=0
	
	for base in fasta:
		b=base.upper()
		if b == 'A':
			A_count+=1
		elif b == 'C':
			C_count+=1
		elif b == 'G':
			G_count+=1
		elif b == 'T':
			T_count+=1

	ACGT=A_count+C_count+G_count+T_count
	C=C_count
	A=A_count
	cprop=(C/ACGT)*100


	return cprop






def top_and_tail(fastafile, size):

	pyfastaq.utils.syscall("rm -rf ends.fa")
	with open("ends.fa", 'w') as fout:
		for seq_record in SeqIO.parse(fastafile, "fasta"):
			scflen=len(seq_record.seq)
			s=int(size)+1
			sq1=seq_record.seq[0:s]
			sq2=seq_record.seq[scflen-s:scflen]
			fout.write('>'+seq_record.id+'L\n'+str(sq1)+"L\n"+'>'+seq_record.id+'R\n'+str(sq2)+"\n")
	#Below method rejected on grounds that line lengths may be 200mb and not necessarily 60bp!		
	#cmd="grep -A"+size+" -B"+size+" '>' "+fastafile+" > ends.fa"

	
	

def find_jellyfish():

	cmd="whereis jellyfish"
	loc=subprocess.getoutput(cmd)
	try:
		jellyfish=loc.split()[1]
		return jellyfish
	except:
		return None




def process_fasta(fasta, num):

	results={}
	final={}
	ends=[]
	with open(fasta,'r') as f:
		for line in f:
			x=line.strip().split()
			if x!=[]:
				if x[0][0]==">":
					header=x[0].replace(">","")[:-1]
					if header not in results:
						results[header]=[[],[]]
				elif results[header]==[[],[]]:
					results[header][0]=x
				elif results[header][0]!=[]:
					results[header][1]=x
					
					
				

					
	for item in islice(results, num):
		ends.append(item)
	print("\nLooking at head and tail of the following scaffs:\n")
	for k,v in results.items():
		if k in ends:
			print(k.replace(">",""))
			a="".join(v[0])
			b="".join(v[1])
			final[k+"L"]=a.upper()	#ensure bases are upper case
			final[k+"R"]=b.upper()


	return final
					


def write_split_fastas(splits):

	fastas={}
	for k,v in splits.items():
		with open(k.replace(">",""),'w') as fout:
			fout.write(">"+k+"\n")
			fout.write(v)
			fastas[k]=v
	return fastas		
			
			
def run_jelly(splits,it):

	result={}
	inter={}
	final={}
	for k in splits.keys():
		result[k]=[]
		sf=k.replace(">","")
		pyfastaq.utils.syscall("jellyfish count -L 3 -m "+it+" -s 100000 "+sf+" -o "+sf+".jf")
		a=subprocess.getoutput("jellyfish dump -c "+sf+".jf | sort -k2,2n").split("\n")
		for i in a:
			if len(i.split())==2:
				x=int(i.split()[1])
				result[k].append(i)
				


	for k,v in result.items():
		if len(v)>0:
			for i in v:
				kmer = i.split()[0]
				kmer_count = int(i.split()[1])
				if kmer[0]=="A" or kmer[0]=="C":
					kmer=rev_comp(kmer)
				if kmer not in inter:
					inter[kmer]=[kmer_count,1]
				else:
					inter[kmer][0]+=kmer_count
					inter[kmer][1]+=1	

	for k,v in inter.items():
		klen=len(k)
		khits=v[0]
		kends=v[1]
		if khits >5 and kends >2:
			final[k]=v



	#These are our telo candidate kmers before we check for tandem-ness
	return final
	
	
def rev_comp(dna):
	
	y=dna.upper()
	x=y.replace('A','t').replace('T','a').replace('C','g').replace('G','c')	 

	return x[::-1].upper()



def tmp_setup():

	pyfastaq.utils.syscall("mkdir -p tmp_jellyfish")
	os.chdir("tmp_jellyfish")


def check_tandemness(candidate_kmers, fastas):

	tandMers=[]
	results=[]
	temp={}
	
	#candidate_kmers - key = kmer, value = [total_counts,end_countss]
	
	
	for k,v in candidate_kmers.items():
		
		t=k*50
		tandem_bases=30
		tt=t[:tandem_bases]+"-"+k	#taking a slice of our 50kmer artificial tandem probe	
		tandMers.append(tt)
						
	for k,v in fastas.items():
		for tm in tandMers:
			t=tm.split("-")
			if cprop(t[0])>30:	#ensure telomeres are flipped towards TG composition rather than CA
				t=[rev_comp(t[0]),rev_comp(t[1])]
			#record number of matches of tandem-mer in temp
			#t[0]=tandem-mer
			#v.count tandem-mer is how many times it's found in that end
			hits=v.count(t[0])+rev_comp(v).count(t[0])
			if t[0] not in temp:
				temp[t[0]]=hits
			else:
				temp[t[0]]+=hits
				
	
					
	#if matches of tandem-mer exceeds 10, it's unlikely to be simply seq-error				
	for k,v in temp.items():
		if v > 10:
			for tm in tandMers:
				if k==tm.split("-")[0]:
					results.append(tm.split("-"))	
	
	check=[]
	check2=[]
	check3=[]
	removed=[]
	tandems={}
	tandemstring=""
	div=[2,3,5,7]
	for r in results:
		check.append(r)
	#removing kmers nested within longer kmers
	for r in results:
		for c in check:
			for d in div:
				if len(r[1])==len(c[1])/d and c[1].find(r[1]) != -1:
					#print("( removed candidate: ",c,")\n")
					removed.append(c)	
				elif c not in removed:
					if c not in check2:
						check2.append(c)
						
	#checking we don't duplicate tandem-mers
	for c in check2:	
		if tandemstring=="":
			tandemstring=c[0]
			tandems[c[0]]=[]
		if c[0][8:22] not in tandemstring and rev_comp(c[0][8:22]) not in tandemstring:
			tandems[c[0]]=[]
			tandemstring+=c[0]

				
	for c in check2:
		for k,v in tandems.items():
			if c[1] in k:
				if c[1] not in v:
					v.append(c[1])

	canonical_output=[]
	can=[]
	noncanonical_output=[]
	
	print("="*54)
	print("\n      ","R","E","S","U","L","T","S","\t", sep=" -- ")
	print()
	for k,v in tandems.items():
		m=len(min(v, key=len))
		for ca in canonical:
			if ca in k or rev_comp(ca) in k:
				canonical_output.append("="*54+"\n")
				line="LIKELY_TELO-REPEAT:\t"+k
				if line not in canonical_output:
					canonical_output.append(line)
				can.append(k)
				for i in v:
					khits=candidate_kmers[i][0]
					kends=candidate_kmers[i][1]
					if len(i) == m:
						line1="telo_kmer:\t\t"+i
						line2=i+"_count_at_ends:\t"+str(khits)
						line3=i+"_scaff_ends_count:\t"+str(kends)
						if line in canonical_output:
							if line1 not in canonical_output:
								canonical_output.append(line1)
							if line2 not in canonical_output:	
								canonical_output.append(line2)
							if line3 not in canonical_output:
								canonical_output.append(line3)
					
						
	noncanonical_output.append("\nOTHER_RESULTS:\n")				
	for k,v in tandems.items():
		m=len(min(v, key=len))
		if k not in can:
			noncanonical_output.append("-"*54+"\n")
			line="candidate_telo-repeat:\t"+k
			if line not in noncanonical_output:
				noncanonical_output.append(line)
			for i in v:
				khits=candidate_kmers[i][0]
				kends=candidate_kmers[i][1]
				if len(i) == m:
					line1="candidate_telo_kmer:\t"+i
					line2=i+"_count_at_ends:\t"+str(khits)
					line3=i+"_scaff_ends_count:\t"+str(kends)
					if line1 not in noncanonical_output:
						noncanonical_output.append(line1)
					if line2 not in noncanonical_output:	
						noncanonical_output.append(line2)
					if line3 not in noncanonical_output:	
						noncanonical_output.append(line3)
						
	

	########### short block of code removes tandem-mers based on overly long kmers that
	#happen to occur in the sequence fasta due to errors, typically slurring of mono-runs
	#In testing, first tandem-mer is the high confidence one, 2nd one is the one with seq errors which we throw away
	s=[i for i, n in enumerate(canonical_output) if n == "="*54+"\n"]
	for i in s:
		if canonical_output[i-1][0:6]=="LIKELY":
			del(canonical_output[i-3:i])
	###########						
						
	with open("../canonical.txt","w") as fout:
		if len(canonical_output)>0:
			for c in canonical_output:
				print(c)
				if not "====" in c:
					fout.write(c+"\n")
		else:
			print("="*54+"\n")
			print("\nNo canonical telomere motif found\n")
	with open("../non_canonical.txt","w") as fout:
		if len(noncanonical_output)>1:
			for n in noncanonical_output:
				print(n)
				if not "----" in n and not "====" in n:
					fout.write(n+"\n")
		else:
			print("-"*54+"\n")
			print("\nNo non-canonical telomere motif found\n")
					
						

		
def location(f):

	cmd="readlink -f "+f
	b=subprocess.getoutput(cmd)
	
	return b
	

				
def main():


	parser = argparse.ArgumentParser(description='Finds most likely telomere motif in Hi-fi or equivalent quality assembly where telomeres are expected to occur at the ends of multiple scaffolds',formatter_class=argparse.ArgumentDefaultsHelpFormatter) 

	#positional args
	parser.add_argument('fasta', metavar='fasta', type=str, help='assembly fasta')

	
	#optional args
	parser.add_argument('--size', metavar='size', type=int, help='top and tail scf bp', default=200)
	parser.add_argument('--klo', metavar='klo', type=int, help='min kmer', default=4)
	parser.add_argument('--khi', metavar='khi', type=int, help='max kmer', default=15)
	parser.add_argument('--ends', metavar='ends', type=int, help='ends to scan', default=1000)

	#display help when misusage
	if len(sys.argv) <1: 
		parser.print_help()
		

	args = parser.parse_args()  #gets the arguments
	start_time = datetime.now()
	
	print("""\n\n\nFinds most likely telomere motif in Hi-fi or equivalent quality assembly where telomeres \n
	are expected to occur at the ends of multiple scaffolds\n\n""")

	
	fasta=location(args.fasta)
	
	print("\nTop and tailing input scaffold fasta...\n")
	
	#Writes fasta file of just ends (ends.fa) from input fasta file (args.fasta)
	top_and_tail(fasta, args.size)
	fasta=location("ends.fa")
	


	if find_jellyfish() != None:
		jellyfish=find_jellyfish
	else:
		print("\njellyfish is not in your path!\n")

	candidate_kmers={}
	tmp_setup()
	if args.ends:
		ends=int(int(args.ends)/2)

	split_fastas=process_fasta(fasta, ends)

	num_ends=int(len(split_fastas)/2)
	

	fastas = write_split_fastas(split_fastas)
	for i in range(int(args.klo),int(args.khi)+1):
		start = datetime.now()
		print("\nRunning jellyfish for head and tail of "+str(num_ends)+" scaffolds...kmer length "+str(i)+"\n")
		kmers=run_jelly(split_fastas, str(i))
		for k,v in kmers.items():
			if cprop(k)>30:
				k=rev_comp(k)
				if k not in candidate_kmers:
					candidate_kmers[k]=v
			else:
				if k not in candidate_kmers:
					candidate_kmers[k]=v
			#found case where jellyfish was calling a 10mer as an 11mer and so missing the array when I built it
			#so here knock longer kmers down by 1
			if len(k)>9:
				candidate_kmers[k[:-1]]=v
		end = datetime.now()
		print(end-start)
		
	#for k,v in candidate_kmers.items():
	#	print(k,v)

	
			
	print("\nchecking tandemness of candidates...\n")	
	check_tandemness(candidate_kmers, fastas)
	print("\ncleaning up\n")
	pyfastaq.utils.syscall("cd ../ ; rm -rf tmp_jellyfish")
	
	end_time=datetime.now()
	print('\n\nFINISHED:\t{}'.format(end_time - start_time)+"\n\n")
	

if __name__ == '__main__':
	main()





