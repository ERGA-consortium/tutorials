#!/usr/bin/env python
#Alan Tracey 2021

'''
MIT License
 
Copyright (c) 2021 Genome Research Ltd.
 
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

import sys
import pyfastaq
import argparse


parser = argparse.ArgumentParser(description='Produces a pseudo-AGP file from a TPF, giving the user cumulative coordinates for their TPF') 

parser.add_argument('tpf', metavar='tpf', type=str, help='curated tpf')

#display help when misusage
if len(sys.argv) <1: 
	parser.print_help()

args = parser.parse_args()  #gets the arguments

prevscaff=""
count=0
with open(args.tpf, 'r') as f:

	cmd="perl /software/grit/projects/vgp_curation_scripts/test_tpf_sanity.pl -scafflevel "+args.tpf
	pyfastaq.utils.syscall(cmd)
	for line in f:
		x=line.strip().split()
		if x[0] !="GAP":
			scaff=x[2]
			if  prevscaff!=scaff:
				count=1
				total=1		#If new scaffold, set coords to one and count up again
				lo=int(x[1].split(":")[1].split("-")[0])
				hi=int(x[1].split(":")[1].split("-")[1])
				tscaff=x[1].split(":")[0]
				if x[3]=="PLUS":
					orient="+"
				else:
					orient="-"
				amount=hi-lo
				x.append(str(total))	#start coord
				total+=amount
				prevscaff=scaff
				x.append(str(total))	#end coord
				new1=[scaff,x[-2],x[-1],str(count),"W",tscaff,str(lo),str(hi),orient]
				new="\t".join(new1)
				print(new)
			else:			#otherwise, add one and keep counting up from previous value
				count+=1
				total+=1
				lo=int(x[1].split(":")[1].split("-")[0])
				hi=int(x[1].split(":")[1].split("-")[1])
				tscaff=x[1].split(":")[0]
				if x[3]=="PLUS":
					orient="+"
				else:
					orient="-"
				amount=hi-lo
				x.append(str(total))
				total+=amount
				prevscaff=scaff
				x.append(str(total))
				new1=[scaff,x[-2],x[-1],str(count),"W",tscaff,str(lo),str(hi),orient]
				new="\t".join(new1)
				print(new)
		else:
			count+=1
			total+=1
			#x.append(str(total))
			try:
				total+=int(x[2])
				new1=[prevscaff,str(total+1),str(total+int(x[2])),str(count),"U",x[2],"scaffold\tyes\tcuration_agp"]
			except:
				total+=100
				
			#x.append(str(total))
			#new="\t".join(x)
				new1=[prevscaff,str(total+1),str(total+100),str(count),"U","100","scaffold\tyes\tcuration_agp"]
			new="\t".join(new1)
			print(new)
