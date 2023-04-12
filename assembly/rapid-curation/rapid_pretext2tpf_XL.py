#!/usr/bin/env python3
#Alan Tracey 2021

import sys
import argparse
#import re
import pyfastaq
import subprocess
from subprocess import check_output
from datetime import datetime
from time import gmtime, strftime
import time
import os


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



#This 'XL' version of rapid_pretext2tpf is written to deal with larger genomes (>2gb) and fragmented genomes (>400 breaks)
#where pretext resolution may confound attempts to correctly deduce break coordinates purely from inspection of the
#TPF.  This version expects the curator to mark the GAP lines where breaks are intended by inserting a '>' at the beginning
#of the GAP line.

'''


#Global parameters
#=================

netsize=3	#throw a wide net to find tpf coords matching agp breaks - but not too wide (in testing >4 misses breaks in highly fragmented genomes - this just means that the tpfchunks stay together rather than splitting fully to match the agp)
lowcutoff=1.5	#how many texels smaller than the smallest contig in this scaffold we are prepared to go down to look for small contigs (if too small we may be misled by agp fragment artifacts).  Had to raise to 1.2 for idMelMell2_1
prefix="RL_"
sex=["X","Y","Z","W"]
borderlen=80
errors={}
warnings={}
info={}
#NB there is no wholesale texel length cutoff, lowcutoff performs this role more intelligently on a per scaffold basis


def append_dict(k,v,d):

	if k not in d:
		d[k]=[v]
	else:
		d[k].append(v)
	return d



#Takes dictionary as input
def genome_size(tpfdict):


	total=0
	totals={}
	for k,v in tpfdict.items():
		for line in v:
			x=line.split()
			if not "GAP" in line:
				scaff=x[1].split(":")[0]
				hi=int(x[1].split("-")[1])
				lo=int(x[1].split(":")[1].split("-")[0])
				amount=(hi-lo)+1
				if scaff not in totals:
					totals[scaff]=amount
				else:
					totals[scaff]+=amount
				
	for k,v in totals.items():
		total+=v
		
	texel=int(round(total/32768,0))	#Not necessarily the true texel length.  This is used as a standard unit of the genome
	#This is not the length reported in the stats or used to calculate resolution

	return total,texel
	
	

#Takes line list as input	
def genome_size2(output):

	total=0
	totals={}
	for k,v in output.items():
		for line in v:
			x=line.split()
			if not "GAP" in line:
				scaff=x[1].split(":")[0]
				hi=int(x[1].split("-")[1])
				lo=int(x[1].split(":")[1].split("-")[0])
				amount=(hi-lo)+1
				if scaff not in totals:
					totals[scaff]=amount
				else:
					totals[scaff]+=amount

				
	for k,v in totals.items():
		total+=v

	return total
	
	
	
	
def complement_scaffold(scaff_line_list):
	x=scaff_line_list[::-1]
	newlines=[]
	for line in x:
		if "PLUS" in line:
			newlines.append(line.replace("PLUS","MINUS"))	
		elif "MINUS" in line:
			newlines.append(line.replace("MINUS","PLUS"))
		else:
			newlines.append(line)	
			

	return newlines
	

	

def scaffs_from_agp(agp,fragsize,scafflens,ctg_lengths,texel):


	sscaffdict={}
	check={}
	discards={}
	agplines={}
	tagdict={}
	sexchrm={}
	painted_scaffs={}
	allsex=[]
	basesex=[]
	nohap_sex=""
	scaff_with_haplo=[]
	pv_vers=""
	pv_res=0	#If older version of PretextView, no value will be set so make a default here

	
	with open(agp, 'r') as f:
		for line in f:
			if line[0] != "#" and line != "\n":
				x=line.strip().split()
				superscaff=x[0]
				cumulativelow=x[1]
				cumulativehigh=x[2]
				frag=int(cumulativehigh)-int(cumulativelow)
				scaff=x[5]
				orientation=x[8]
				entry=x[4]
				if entry != "U":	#is not a gap
					low=int(x[6])
					high=int(x[7])
					vals=[scaff,low,high,orientation,cumulativelow,cumulativehigh]
					append_dict(superscaff,line,agplines)
					if frag > 10*texel:	#always take agp frags above this size - never get an artefact bigger than 10 texels
						append_dict(superscaff,vals,sscaffdict)
					else:	#If fragment is small...
						if frag > min(ctg_lengths[scaff])-lowcutoff*texel:	#take if small but bigger than smallest contig by a margin
							append_dict(superscaff,vals,sscaffdict)
						else:			
							append_dict(scaff,low,discards)		#discard everything else
							
					#Get the tags
					s=[]
					if len(x)>9:
						unqkey="/".join([scaff,str(low),str(high)])
						for i in range(len(x[9::])):
							if unqkey not in tagdict:
								tagdict[unqkey]=[[x[9::][i].upper()]]
							else:
								tagdict[unqkey][0].append(x[9::][i].upper())
							s.append(x[9:][i].upper())
							
						if "HAPLOTIG" in s and "UNLOC" in s:
							msg=scaff+" is both \'Haplotig\' and \'unloc\' - needs to be one or the other\n"
							append_dict("hap_unloc",msg,errors)
						if "HAPLOTIG" in s:
							scaff_with_haplo.append(superscaff)
							
						for i in s:
							if "W" in i[0]:
								nohap_sex="W"
							if "Y" in i[0]:
								nohap_sex="Y"	
							for sx in sex:
								if sx in i[0]:
									if i not in allsex:
										allsex.append(i)
										
									if i not in sexchrm:
										sexchrm[i]=[superscaff]	#Starting a check to see if eg Z chr is referenced in >1 chrm
									else:
										if superscaff not in sexchrm[i]:
											sexchrm[i].append(superscaff)
						
						if "PAINTED" in s:	
							painted_scaffs[superscaff]="PAINTED"	
						
								
						
						
						if superscaff not in check:
							check[superscaff]=[s,[scaff]]
						else:
							for i in x[9:]:
								if i.upper() not in check[superscaff][0]:
									check[superscaff][0].append(i.upper())
							check[superscaff][1].append(scaff)
							
			else:	#AGP header
				if "DESCRIPTION: Generated by PretextView Version" in line:
					pv_vers=line.split()[-1]
				if "HiC MAP RESOLUTION:" in line:
					pv_res=round(float(line.split()[-2]),2)	#bp/texel
					

								
						
	for k,v in check.items():
		if "UNLOC" in v[0] and len(v[1])==1 and len(painted_scaffs)>=1:
			print("\nUnlocalised scaff is only scaff in this chromosome\n")
			print("AGP scaffold "+k+" , TPF scaffold "+v[1][0]+"\n")
			sys.exit()			
						

	
	
					
						
	for a in allsex:
		if a[0] not in basesex:
			basesex.append(a[0])					
					
	if "Z" in basesex and "X" in basesex:
		msg = "Bad sex combination found "+" ".join(allsex)
		append_dict("badsex",msg,errors)
	elif "W" in basesex and "X" in basesex:
		msg = "Bad sex combination found "+" ".join(allsex)
		append_dict("badsex",msg,errors)
	elif "Z" in basesex and "Y" in basesex:
		msg = "Bad sex combination found "+" ".join(allsex)
		append_dict("badsex",msg,errors)
	elif len(sexchrm)==2:
		for k,v in sexchrm.items():
			if v[0] in scaff_with_haplo:
				append_dict("hetero_sex_haplo",k,errors)						


			
								
	#print(tagdict)						
	multisex2={}
	for k,v in sexchrm.items():
		if len(v)>1:
			length=str(len(v))
			msg=k+" chromosome is referenced in "+length+" AGP chromosomes:\n"
			append_dict("multiple_sex",msg,errors)
			for i in v:
				append_dict("multiple_sex",i,errors)
		for i in v:
			append_dict(i,k,multisex2)
	
	for k,v in multisex2.items():
		if len(v)>1:
			msg=k+" chromosome is referenced as both "+v[0]+" and "+v[1]+" in AGP:\n"
			append_dict("multiple_sex2",msg,errors)
			





	return sscaffdict, discards, agplines, tagdict, sexchrm, painted_scaffs, pv_vers, pv_res


				
	


def tpf_sanity(tpfout):
	
	try:
		cmd="perl "+os.path.dirname(sys.argv[0])+"/test_tpf_sanity.pl -scafflevel "+tpfout
		pyfastaq.utils.syscall(cmd)
	except:
		print("\nUnable to find test_tpf_sanity.pl - please copy this to script location\n")



def location(f):

	cmd="readlink -f "+f
	b=subprocess.getoutput(cmd)
	
	return b
	
	
#Gets tpf lines and checks for typos in coordinates (some coordinates may be added manually so we check them)
def parse_tpf(tpf):

	prevscaff=""
	result={}
	coordtest={}	#Checking for obvious coord typos in the input tpf
	errors2=[]
	with open(tpf,'r') as f:

		for line in f:
			if not "gap" in line.lower():
				x=line.split()
				scaff=x[2]
				lo=x[1].split(":")[1].split("-")[0]
				hi=x[1].split(":")[1].split("-")[1]
				if scaff in coordtest:
					if int(lo)<=max(coordtest[scaff]):	#If coord is less than previous line
						if line not in errors2:
							errors2.append(line)
							
					if int(hi)<=max(coordtest[scaff]):
						if line not in errors2:
							errors2.append(line)
							
					if int(lo)>=int(hi):
						if line not in errors2:
							errors2.append(line)
							
				append_dict(scaff,int(hi),coordtest)
				if scaff not in result:
					result[scaff]=[line]
					prevscaff=scaff
				else:
					result[scaff].append(line)
					prevscaff=scaff
			else:
				result[prevscaff].append(line)
				

				

	return result, errors2



def report_errors(errors2):

	if len(errors2)>0:
		print("coordinate errors detected amongst the following tpf lines - presumed typos:\n")
		for e in errors2:
			print(e)
		sys.exit()


	
def components_from_dict(tpf):

	check=[]
	for k,v in tpf.items():
		for l in v:
			if not "GAP" in l:
				comp=l.split()[1]
				check.append(comp)

	return check
	
	
#The coordinates where our agp fragments terminate
def agp_dividers(agpdict):

	tdivs={}
	result={}
	for k,v in agpdict.items():
		for i in v:
			tscaff=i[0]
			tmin=int(i[1])
			tmax=int(i[2])
			#print(k,tscaff,tmax)
			vals=[tmin,tmax]
			append_dict(tscaff,vals,tdivs)
				
	
	for k, v in tdivs.items():
		result[k]=sorted(v)
		
		
	return result





def lens(tpfdict):

	lens={}
	for k,v in tpfdict.items():
		for i in v:
			#print(k,i)
			if not "GAP" in i:
				hi=int(i.split()[1].split("-")[1])
				if k not in lens:
					lens[k]=hi
				elif hi > lens[k]:
					lens[k]=hi
					
	
	return lens
	


def write_output_tpf(outlinesfull,outfile):
 
	with open(outfile,'w') as fout:
		for k,v in outlinesfull.items():
			for line in v:
				if "AMBIG" not in k:
					fout.write(line+"\n")
	with open(outfile,'a') as fout:
		for k,v in outlinesfull.items():
			for line in v:
				if "AMBIG" in k:
					fout.write(line+"\n")




def compare_scaff(tpfdict,agp):


	scaffs=[]
	probs=[]
	with open(agp, 'r') as f:
		for line in f:
			if line[0] != "#" and line != "\n":
				x=line.strip().split()
				scaff=x[5]
				entry=x[4]
				if entry != "U":	#is not a gap
					scaffs.append(scaff)

	for s in scaffs:
		if s not in tpfdict:
			probs.append(s)
	if len(probs) > 0:
		print("\nagp and tpf not in sync, ",probs[0],"not in tpf for example\n")
		sys.exit()
		

def commas(number):

	n=list(str(number))
	for i in range(len(n))[::-3]:
    		n.insert(i+1,",")
	if n[-1]==",":
		del n[-1]
	if n[0]==",":
		del n[0]
	
	return "".join(n)
	
	
def contig_lens(tpfdict):

	ctg_lengths={}
	for k,v in tpfdict.items():
		for i in v:
			if not "GAP" in i:
				x=i.strip().split()
				hi=int(x[1].split(":")[1].split("-")[1])
				lo=int(x[1].split(":")[1].split("-")[0])
				length=hi-lo+1
				append_dict(k,length,ctg_lengths)
			
	return ctg_lengths
	
def report_agp_discards(discards,agplines):


	for k,v in agplines.items():
		for line in v:
			i=line.strip().split()
			scaff=i[5]
			coord=i[6]
			if scaff in discards:
				#print(scaff,coord)
				#print(type(discards[scaff][0]),type(coord))
				if discards[scaff][0]==int(coord):
					print("AGP_discard",line.strip())



def check_input_tpf(tpf):

	#We make a temporary tpf with ">" characters removed so we can check tpf sanity (">" are needed to mark breaks, but would make tpf appear corrupt)
	tempout="tempfile"
	with open(tpf,'r') as f:
		with open(tempout,'w') as fout:
			for line in f:
				if line[0]==">":
					newline=line[1:]
					fout.write(newline)
				else:
					fout.write(line)
	
	print("\nChecking "+tpf+" sanity...\n\n")
	tpf_sanity(tempout)
	cmd="rm -rf "+tempout
	pyfastaq.utils.syscall(cmd)	#Get rid of the temporary file we've created

	
	
#Gets breaks the curator has flagged up with ">" in the GAP line	
def get_breaks(tpf):

	breaks={}
	lastline=[]
	with open(tpf,'r') as f:
		scaff=""
		for line in f:
			x=line.strip().split()
			if not "GAP" in line:
				scaff=x[2]
				mx=int(x[1].split(":")[1].split("-")[1])
				if len(lastline)==0:
					lastline.append(x)
				else:
					lastline[0]=x
			else:
				if line[0]==">":
					append_dict(scaff,mx,breaks)
					
	if len(breaks)==0:
		msg='>>>  Expected '+tpf+' to contain breaks marked with ">" but none were found\t<<<\n>>>\t\t\tremember this is the XL script\t\t\t\t\t<<<'
		append_dict("no_breaks",msg,warnings)

	return breaks
				
			
#output will be a dict of numbered chunks to a list of scaff lines, eg scaffold_44%2 ['?\tscaffold_44:20173615-20272426\tscaffold_44\tPLUS\n']			
def breaktpf(tpfdict,breaksnew):

	#for k,v in breaksnew.items():
		#print(k,v)	
		
	tpfchunks={}
	for k,v in tpfdict.items():
		iteration=0
		if k not in breaksnew:
			tpfchunks[k+"%1"]=v	#We add all our unbroken tpfchunks
		else:	
			iteration+=1
			key=k+"%"+str(iteration)
			for line in v:
				append_dict(key,line.strip(),tpfchunks)
				x=line.strip().split()
				if not "GAP" in line:
					hi=int(x[1].split(":")[1].split("-")[1])
					if hi in breaksnew[k]:
						iteration+=1
						key=k+"%"+str(iteration)
						
	for k,v in tpfchunks.items():
		if "GAP" in v[0]:
			del v[0]
		if "GAP" in v[-1]:
			del v[-1]

					
					
	return tpfchunks
	
	
def ranges(tpfchunks):

	check=[]
	keep=[]
	results={}
	final={}
	for k,v in tpfchunks.items():
		if k not in results:
			results[k]=[[],[]]
		for line in v:
			if not "GAP" in line:
				x=line.strip().split()
				lo=int(x[1].split(":")[1].split("-")[0])
				hi=int(x[1].split(":")[1].split("-")[1])
				#print(k,lo,hi)
				if results[k][0]==[]:
					results[k][0]=lo
				elif lo < results[k][0]:
					results[k][0]=lo
				if results[k][1]==[]:
					results[k][1]=hi
				elif hi > results[k][1]:
					results[k][1]=hi
				#print(results[k][1],hi)
					

		c1=k.split("%")[0]
		if c1 not in check:
			check.append(c1)
		else:
			keep.append(c1)
	
	for k,v in results.items():
		if k.split("%")[0] in keep:
			final[k]=v
	#for k,v in final.items():
	#	print(k,v)		

	#final only contains tpfchunks resulting from breaks		
	return final
				
	
	
#Each one of these must have only 1 chunk so unique occurences in agp can always be paired
def get_unbroken_scaff(tpfchunks):

	counts=[]
	results=[]
	for k,v in tpfchunks.items():
		counts.append(k.split("%")[0])
	for count in counts:
		if counts.count(count)==1:
			results.append(count+"%1")
	return results	#List of all unbroken scaffolds	
		


#Attempt to uniquely assign our tpf chunks to agp and report any issues
def maptpf2agp(cleanedranges,tpfchunks,agpcomps_map,prefix,unbroken,scaff_dict,texel,scaffold_lengths,shrapnel,tpf,tpfout):



	#for k,v in cleanedranges.items():
	#	print(k,v)



	notused=[]	#captures any chunks that fall through the net
	tpfchunklist=[]
	ambiguous=[]
	unambigagp={}
	#for k,v in tpfchunks.items():
		#print(k,v)
		#tpfchunklist.append(k)
	for chunk in unbroken:
		if chunk not in tpfchunklist:
			tpfchunklist.append(chunk)
			
	#for t in tpfchunklist:
		#print(t)

	old_lines=[]
	unambig={}
	
	
	ambig=list(cleanedranges.keys())	#assume everything is ambiguous to begin with
	for u in unbroken:
		if u not in ambig:
			ambig.append(u)
	#print(ambig)
	iteration=1
	prevagp=""
	factor=0.499	#0.499 == the middle
	ln=0		#ln is the agp line number
	prevchr=""
	agplines={}
	newambig={}
	for k, v in agpcomps_map.items():
		for i in v:
			chrm=k
			chrnum=k.split("_")[1]
			pre=prefix+chrnum
			agpscaff=i[0]
			agplow=min(i[1],i[2])
			agphigh=max(i[1],i[2])
			orientation=i[3]
			if prevchr==k:
				ln+=1
				link=k+":"+orientation+":"+str(ln)	#ln is the agp line number
				prevchr=k
			else:
				ln=1
				link=k+":"+orientation+":"+str(ln)
				prevchr=k
			#print(link)
			append_dict(link,i,agplines)
			for a, b in cleanedranges.items():	#cleanedranges are the tpf chunks ranges where tpf is broken
				tpfscaff=a.split("%")[0].replace("_ctg1","")
				if tpfscaff==agpscaff:		#scaffold_21
					#print(k,a,b)		
					tpfchunk=a		#scaffold_95_ctg1%1
					tpfchunklow=b[0]
					tpfchunkhigh=b[1]
					size=tpfchunkhigh-tpfchunklow
					adjust=int(round(size*factor,0))	#wiggle room as we know pretext break coord won't be accurate (adjusting below 0.05 misses some for sure, higher is more permissive)
					tpflow=tpfchunklow+adjust	
					tpfhigh=tpfchunkhigh-adjust
					deltalo=abs(agplow-tpfchunklow)
					deltahi=abs(agphigh-tpfchunkhigh)
					#print(a,tpflow,tpfhigh,tpfchunklow,tpfchunkhigh,agplow,agphigh,texel,deltalo, deltahi)
					if tpflow > agplow and tpfhigh < agphigh:	#This is near the middle
						#print(a,tpflow,tpfhigh,tpfchunklow,tpfchunkhigh,agplow,agphigh,texel,deltalo, deltahi)
						if deltalo < 100*texel and deltahi < 100*texel:	#if pretext breaks are less than 10 texels from our actual tpf breaks...
							#Putting this value up to 100 (effectively removing the check - no obvious problems with this found and results are improved with iyDolMedi)
							#Previous value was 5 in line above which worked for 5 test datasets.  This failed to find some chunks in larger genome (idCorMarg)
							#Going to 8 on idCorMarg dataset caused 2nd match warning A and the script breaks.  Fixed this by changing 'size' to tpfchunk size (was agp chunk size)
							#Such code could (instead of delta filter above, check on the expected size of the agp chunk and ensure that (for just cases that fail above check)
							#that the tpf chunk goes with the closest agp chunk size		
							#print(a, tpflow, tpfhigh, agplow, agphigh)
							if a not in unambig:
								vals=[tpflow, tpfhigh, agplow, agphigh]
								#print("T",a,vals)
								unambig[a]=[vals]
								unambigagp[link+"#"+a]=[vals]
								ambig.remove(a)
							else:
								#unambig[a].append(vals)
								print("2nd match - warning A", a)	#No cases of this found - it may never happen - wait on more test data
								newambig[a]=unambig[a]		#Recording the details of any theoretical 2nd match cases we might come across for further investigation
								newambig[a].append(vals)
								if a not in ambig:
									ambig.append(a)
								del unambig[a]
								if link+"#"+a in unambigagp:
									del unambigagp[link+"#"+a]


			test=agpscaff+"%1"
			if test in unbroken:	#These entries are unique so the coord check could actually be turned off below
				#print("here",test,i,agplow)	#works
				scaffkey=test.split("%")[0]
				if test not in unambig:
					#print(test)
					size=agphigh-agplow
					#print("A",size)
					adjust=int(round(size*factor,0))
					tpflow=1+adjust
					tpfhigh=scaffold_lengths[scaffkey]-adjust-size
					#print(scaffold_lengths[scaffkey])
					#print(link,i,tpflow,agplow,tpfhigh,agphigh)
					if tpflow > agplow and tpfhigh < agphigh:
						#print(test)
						vals=[1,scaffold_lengths[scaffkey],agplow,agphigh]
						#print(vals)
						unambig[test]=[vals]
						unambigagp[link+"#"+test]=[vals]
						ambig.remove(test)
						#print(test,1,scaffold_lengths[scaffkey])
				else:
					size=agphigh-agplow
					adjust=int(round(size*factor,0))
					tpflow=1+adjust
					tpfhigh=scaffold_lengths[scaffkey]-adjust-size
					#print("RR",k,test,tpflow,agplow,tpfhigh,agphigh)
					if tpflow > agplow and tpfhigh < agphigh:
						print("2nd match - warning B", test)
						vals=[1,scaffold_lengths[scaffkey],agplow,agphigh]
						if test not in ambig:
							ambig.append(test)

						del unambig[test]
						try:
							del unambigagp[link+"#"+test]
						except:
							pass
	
	#for k, v in unambig.items():
		#print("unambig", k,v)
				
				
	for k,v in newambig.items():	#This is a check that has never been triggered
		print("newambig",k,v)
		
	unbrokentest=[]
	for u in unbroken:
		if u not in unbrokentest:
			unbrokentest.append(u)
	
	for a in ambig:
		#print("ambig",a)
		for k,v in agpcomps_map.items():
			for i in v:
				#print(a,i)
				test=a.split("%")[0].replace("_ctg1","")
				if test== i[0]:
					if test in unbrokentest:
						unbrokentest.remove(test)
	for u in unbrokentest:
		if u in ambig:
			ambig.remove(u)	#Can't be ambiguous if it's not been broken - just put it in where scaffold is mentioned in agp
			
	if len(ambig)==0:
		msg="ambiguous tpf to agp (AMBIGs):\t0"
		info["ambig_tpf_to_agp"]=[msg]
	
	else:
		msg="ambiguous tpf to agp (AMBIGs):\t"+str(len(ambig))+"    <--- (check "+tpf+" & "+tpfout+" tail)"
		info["ambig_tpf_to_agp"]=[msg]
		msg="\nTPF chunks that don't match to AGP:\n"
		warnings["ambig"]=[[msg],[],[]]
		for a in ambig:
			scaff=a.split("%")[0].replace("_ctg1","")
			tpfvals=[scaff,cleanedranges[a]]
			lo=int(cleanedranges[a][0])
			hi=int(cleanedranges[a][1])
			diff=hi-lo
			mid=(diff/2)+lo
			warnings["ambig"][1].append(tpfvals)
			for k,v in agpcomps_map.items():
				for i in v:
					if scaff==i[0] and mid > i[1] and mid < i[2]:
						warnings["ambig"][2].append(i[0:3])

	

	msg="AMBIG"
	warnings[msg]=ambig

	#print("both ambig and unambig:\t",len(ambig)+len(unambig.keys()))
	#print("tpfchunklist:\t",len(tpfchunklist))
	
	msg="tpfchunks master list:\t\t"+str(len(tpfchunks.keys()))
	info["masterlist"]=[msg]
	
	msg="in tpf but not agp (shrapnel):\t"+str(len(shrapnel))
	info["shrapnel"]=[msg]

	msg="ambig+unambig tpf to agp:\t"+str(len(set(ambig+list(unambig.keys()))))
	info["ambig_plus_unambig"]=[msg]
	
	known=set(ambig+list(unambig.keys())+shrapnel)
	
	unknown=[]
	for t in tpfchunklist:
		if t not in known:
			unknown.append(t)
	#print(unknown)


	accounted=len(set(ambig+list(unambig.keys())+shrapnel))
	msg="ambig+unambig+shrap (1st pass):\t"+str(accounted)

	info["ambig_plus_unambig_plus_shrapnel1"]=[msg]
	

	#for k,v in unambigagp.items():
		#print("unambiguous agp",k,v)
	
		
	#for a in ambig:
		#print(a)
		
	#for s in shrapnel:
		#print(s)
	

	
	#print("unknown:\t", len(unknown))
	
	ulist=[]
	for u in unknown:
		b=u.split("_")[0:2]
		base="_".join(b)
		ulist.append(base)
		#print(ulist)
	
	#IMPORTANT - ADDING IN THE SMALL PIECES WE COULDN'T GET UNAMBIGUOUSLY BY COORDINATES
	ulist4agp=[]
	oldaccounted=accounted
	if len(unknown)>0:
		with open("2nd_pass_pickup.info","w") as fout:
			for u in unknown:
				b=u.split("_")[0:2]
				base="_".join(b)
				if ulist.count(base)==1:	#This has only 1 chunk so we can take the biggest hit in the agp and place it	
					for k,v in agpcomps_map.items():
						for i in v:
							if base==i[0]:
								#print(k,u,i) 
								#Need to place in the largest match in the agp
								ulist4agp.append(i[0])	#this i[0] is what we need
								if base in ulist:
									ulist.remove(base)
									accounted+=1
									for line in tpfchunks[u]:
										fout.write(line+"\n")
									#Here we update the stats to reflect new tiny scaffold inclusions
							
	msg="plus "+str(accounted-oldaccounted)+" (2nd pass):\t\t"+str(accounted)
	info["ambig_plus_unambig_plus_shrapnel2"]=[msg]
							
	
	#for k,v in cleanedranges.items():
		#print(k,v)
		
	#for u in unbroken:	#list of chunk %1 no range
		#print(u)
	
	#for k,v in scaffold_lengths.items():
		#print(k,v)
		
	line=0
	prevchr=""
	for k,v in agpcomps_map.items():
		#print(k,v)
		for i in v:
			agpscaff=i[0]
			orientation=i[3]
			if prevchr==k:
				line+=1
				#print(k,orientation,str(line))
				link=k+":"+orientation+":"+str(line)
				prevchr=k
			else:
				line=1
				#print(k,orientation,str(line))
				#print(k)
				link=k+":"+orientation+":"+str(line)
				prevchr=k
				#print(link)
			if agpscaff in ulist4agp:
				newkey=link+"#"+agpscaff+"_ctg1%1"
				unambigagp[newkey]=["unq match, not found by coords"]
				newkey2=agpscaff+"_ctg1%1"
				unambig[newkey2]=["found"]
				

				
	msg="not in output tpf:\t\t"+str(len(ulist))
	info["not_in_out_tpf"]=[msg]
	#print(ulist)


			
		
	return ambig, unambig, unambigagp, agplines, ulist	#These are tpfchunks assigned unambiguously or not to agp line, agplines gives order in agp
	


def tag_tpfchunks(tagdict,unambigagp):	
		
	#for k,v in tagdict.items():
	#	print(k,v)
	tpfchunktags={}
	for k,v in unambigagp.items():
		tpfchunkkey=k.split("#")[1]
		for a,tags in tagdict.items():
			scaff=a.split("/")[0]
			lo=a.split("/")[1]
			hi=a.split("/")[2]
			if scaff == tpfchunkkey.split("%")[0].replace("_ctg1",""):
				alo=str(v[0][2])
				ahi=str(v[0][3])
				if lo == alo and hi == ahi:
					#print(k,v,tpfchunkkey,tags[0])
					append_dict(tpfchunkkey,tags[0],tpfchunktags)
				elif v[0]=="unq match, not found by coords":
					append_dict(tpfchunkkey,tags[0],tpfchunktags)
				#elif tpfchunkkey not in tpfchunktags and tags[0][0] == "PAINTED":
				#	append_dict(tpfchunkkey,tags[0],tpfchunktags)
				#	#print(lo,alo,hi,ahi,k,v,tpfchunkkey,tags[0])
					
	#Scaffold_2:-:25#scaffold_2%4 [[708303706, 708374746, 689703245, 726572175]] scaffold_2%4 ['PAINTED']	

	#{'scaffold_38_ctg1%1': [['W']], 'scaffold_35_ctg1%1': [['W']], 'scaffold_60_ctg1%1': [['W']], 'scaffold_33_ctg1%1': [['W']], 'scaffold_49_ctg1%1': [['UNLOC', 'Z']], 'scaffold_34_ctg1%1': [['HAPLOTIG']]}
	
	#for k,v in tpfchunktags.items():
	#	print(k,v)
			
	return tpfchunktags
	
	
#breaks dict has component as key with break count per component as value
#unambigagp is agplines with snap mode artefacts removed
def break_count_vs_agp_line_count(unambigagp,tpfchunks,minscaffagpcount):

	agpcount={}	#agplines per scaffold
	tpfcount={}	#tpfchunks per scaffold
	missingbreaks=[]
	probs={}
	probsfinal={}
			

	for k,v in minscaffagpcount.items():
		agpcount[k]=len(v)
		

			
			
	for k,v in tpfchunks.items():
		scaff=k.split("%")[0].replace("_ctg1","")
		if scaff not in tpfcount:
			tpfcount[scaff]=1
		else:
			tpfcount[scaff]+=1		
	
	for k,v in agpcount.items():
		if v != tpfcount[k]:
			if v>tpfcount[k]:
				missingbreaks.append(k)
				
	
	for scaff in missingbreaks:
		temp=[]
		for i in range(len(minscaffagpcount[scaff])):
			val=[minscaffagpcount[scaff][i][1]]
			temp.append(val)
		for k,v in unambigagp.items():
			uscaff=k.split("#")[1].replace("_ctg1","").split("%")[0]
			if scaff==uscaff:
				hi=str(v[0][3])
				for t in temp:
					if hi in t:
						t.remove(hi)	#removing all the non-problems
		for t in temp:
			if t!=[]:
				if scaff not in probs:
					probs[scaff]=[int(t[0])]
				else:
					if int(t[0])> probs[scaff][-1]:
						probs[scaff].append(int(t[0]))
					if int(t[0])<probs[scaff][0]:
						probs[scaff].insert(0,int(t[0]))
	
	#using the probs ID'd on the AGP max coord, we are now going to make a string of lo+"-"+hi to use in message to curator	
	for k,v in probs.items():
		for prob in v:
			for i in minscaffagpcount[k]:
				#print(k,i,prob)
				if str(prob)==str(i[1]):
					vals=i[0]+"-"+i[1]
					append_dict(k,vals,probsfinal)
		
		
	if len(probsfinal)>0:
		append_dict("missing_breaks",probsfinal,warnings)
		

		

	
	return probsfinal	
	


#For problems we've already found, we try to report the cause and location of these
def locate_problems(tpf_vs_agp_discreps, unambigagp,tpfchunks,cleanedranges,texel):

	scaffs=[]	#The scaffs that contain some issues
	probs={}	#difference between this and notprobs is the actual probs and their coords
	notprobs={}

	for i in tpf_vs_agp_discreps:
		scaffs.append(i[0])
	for s in scaffs:
		for k,v in tpfchunks.items():
			tscaff=k.split("%")[0].replace("_ctg1","")
			if s == tscaff:
				if s not in probs:
					probs[s]=[cleanedranges[k]]
				else:
					probs[s].append(cleanedranges[k])




	for s in scaffs:
		for k,v in unambigagp.items():
			ascaff=k.split("#")[1].split("%")[0].replace("_ctg1","")
			if s == ascaff:
				alo=v[0][2]
				ahi=v[0][3]
				if s in probs:
					for p in probs[s]:
						#print(s,p)
						tlo=p[0]
						thi=p[1]
						difflo=abs(alo-tlo)
						diffhi=abs(ahi-thi)
						#print(difflo,diffhi,7*texel)
						if difflo < 7*texel and diffhi < 7*texel:	#The good matches, ie not problems
							if s not in notprobs:
								notprobs[s]=[[tlo,thi]]
							else:
								notprobs[s].append([tlo,thi])
								




	for k,v in probs.items():
		for i in v:
			if k in notprobs:
				if i not in notprobs[k]:	#The real problems
					msg="\nTPF chunks that don't match well with AGP lines:\n"
					if "TPF_problem_locations" not in warnings:
						warnings["TPF_problem_locations"]=[msg]
						warnings["TPF_problem_locations"].append([k,commas(i[0])+"-"+commas(i[1])])
					else:
						warnings["TPF_problem_locations"].append([k,commas(i[0])+"-"+commas(i[1])])
						
			else:	#believed to be shrapnel only so we can just let these pass as we haven't broken these in the AGP
				pass
			
				#msg="\nTPF chunks that don't match well with AGP lines:\n"
				#if "TPF_problem_locations" not in warnings:
				#	warnings["TPF_problem_locations"]=[msg]
				#	warnings["TPF_problem_locations"].append([k,commas(i[0])+"-"+commas(i[1])])
				#else:
				#	warnings["TPF_problem_locations"].append([k,commas(i[0])+"-"+commas(i[1])])



	
#We check that our pieces are bigger than 10 texels, if they are, we check the tpf length is within 10% of the agp length	
def additional_check1(unambigagp,cleanedranges,texel,agplines,tpfchunks,minagpscaffcount, pv_res):


	#unambigagp - unq key, tpf low and high after *0.499 mutliplication (these values are useless), agp low and high - useful
	#Scaffold_2:+:3#scaffold_2_ctg1%3 [[17598223, 17666271, 640459, 34655914]]
	missing_breaks=[]
	excess_breaks=[]
	for k,v in unambigagp.items():
		tkey=k.split("#")[1]
		akey=k.split("#")[0]
		alo=v[0][2]
		ahi=v[0][3]
		tpfchunkstartcomp=tpfchunks[tkey][0].split()[1]
		tpfchunkendcomp=tpfchunks[tkey][-1].split()[1]
		scaff=tpfchunks[tkey][0].split()[1].split(":")[0]
		if tkey in cleanedranges:	#These are our broken scaffolds	
			tlo=cleanedranges[tkey][0]
			thi=cleanedranges[tkey][1]
			tsize=thi-tlo
			asize=ahi-alo
			diff=abs(tsize-asize)	#very useful - used
			maxsize=max(tsize,asize)
			pcntdiff=int(diff/maxsize*100)	#not used
			tex=round(diff/texel,1)		#how many texels we're out by
			truetex=round(diff/pv_res,1)
			if diff>8*texel:	#Non-problematic differences seem generally below 5 texels; Ncomps check - we're not interested in reporting N-only components
				#print(k,tex,tsize,asize,diff,agplines[akey],tpfchunkstartcomp,tpfchunkendcomp,str(tlo),str(thi),str(alo),str(ahi))
				vals=[scaff+"\t"+str(commas(tlo))+"--"+str(commas(thi))+"\t"+str(commas(alo))+"--"+str(commas(ahi))+"\t"+commas(diff)+"\t"+str(truetex)]
				if tsize>asize:
					missing_breaks.append(vals)
				elif asize>tsize:
					#print(vals)
					excess_breaks.append(vals)
				



		else:
			if "unq" not in v[0]:	#These are our unbroken scaffolds
				tlo=v[0][0]
				thi=v[0][1]
				tsize=thi-tlo
				asize=ahi-alo
				diff=abs(tsize-asize)
				maxsize=max(tsize,asize)
				pcntdiff=int(diff/maxsize*100)	#not useful
				tex=round(diff/texel,1)		#how many texels we're out by
				truetex=round(diff/pv_res,1)
				#print(cleanedranges[tkey])
				if diff>8*texel:
					#print(tsize,diff)
					vals=[scaff+"\t"+str(commas(tlo))+"--"+str(commas(thi))+"\t"+str(commas(alo))+"--"+str(commas(ahi))+"\t"+commas(diff)+"\t"+str(truetex)]
					#print(k,tex,tsize, asize, diff, agplines[akey],tpfchunkstartcomp,tpfchunkendcomp,str(tlo),str(thi),str(alo),str(ahi))	#These are unbroken scaffolds (but if break is missing the size will be wrong)								
					if tsize>asize:
						missing_breaks.append(vals)
						#Since these are unbroken scaffolds, we don't need to check for excess breaks
						
	
	msg="\nMismatch between TPF and AGP - missing or inaccurate breaks:\n\nScaff\t\tTPF-range\tAGP-range\tsize diff\ttexels diff\n"
	if len(missing_breaks)>0:
		warnings["missing_breaks2"]=[msg]
		for m in missing_breaks:
			warnings["missing_breaks2"].append(m)

			
	#Full set of excess breaks will only be reported once no missing breaks remain
	msg="\nMismatch between TPF and AGP - excess or inaccurate breaks:\n\nScaff\t\tTPF-range\tAGP-range\tsize diff\ttexels diff\n"
	if len(excess_breaks)>0:
		warnings["excess_breaks"]=[msg]
		for e in excess_breaks:
			warnings["excess_breaks"].append(e)	
	




def prepare_output_tpf(agplines,unambigagp,tpfchunks,tpfout,shrapchunks,ambigtcs):

			

	lastl=""
	prevchr=""
	outlines={}
	unfilled_agp=[]
	unused=[]
	joins=0
	newold={}
	for k,v in tpfchunks.items():
		unused.append(k)

		
	for l in agplines.keys():
		#print(l)
		chrm=l.split(":")[0]
		iteration=chrm.split("_")[1]
		for k,v in unambigagp.items():
			#print(l,k,v)
			a=k.split("#")[0]
			#print(l,a,k,v)
			if l==a:
				tpfchunkkey=k.split("#")[1]
				orientation=k.split(":")[1]
				if prevchr==chrm and "GAP" not in outlines[chrm][-1]:	#ensures we don't add a double gap where a N chunk has been deleted
					outlines[chrm].append("GAP\tTYPE-2\t200")
					joins+=1
				if orientation == "-":
					unused.remove(tpfchunkkey)
					for line in complement_scaffold(tpfchunks[tpfchunkkey]):
						if chrm not in outlines:
							x=line.strip().split()
							x[2]=prefix+iteration
							newline="\t".join(x)
							outlines[chrm]=[newline]
							
							
							#print(newline)
						else:
							x=line.strip().split()
							if not "GAP" in line:
								x[2]=prefix+iteration
								newline="\t".join(x)
								outlines[chrm].append(newline)
							else:
								outlines[chrm].append(line.strip())
							
							#print(newline)
				else:
					unused.remove(tpfchunkkey)
					for line in tpfchunks[tpfchunkkey]:
						if chrm not in outlines:
							x=line.strip().split()
							x[2]=prefix+iteration
							newline="\t".join(x)
							outlines[chrm]=[newline]
							
							#print(newline)
						else:
							x=line.strip().split()
							if not "GAP" in line:
								x[2]=prefix+iteration
								newline="\t".join(x)
								outlines[chrm].append(newline)
							else:
								outlines[chrm].append(line.strip())
							#print(newline)


				prevchr=chrm
				
				#print(k,tpfchunks[tpfchunkkey])
				lastl=l	#<---
		

		if l !=lastl:			#<---
			#print("unfilled AGP lines:",l)		#<---Gives us unfilled lines
			unfilled_agp.append(l)
			

	
	
	count=0
	for a in ambigtcs:
		if a in unused:
			count+=1
			if a in tpfchunks:	#Ensures N-only components aren't used which would cause a crash
				for line in tpfchunks[a]:
					if not "GAP" in line:
						x=line.strip().split()
						x[2]="AMBIG_"+str(count)	
						newline="\t".join(x)
						append_dict(a,newline,outlines)
					else:
						append_dict(a,line.strip(),outlines)
			unused.remove(a)
						
	for s in shrapchunks:
		if s in unused:
			lines=tpfchunks[s]
			for line in lines:
				append_dict(s,line.strip(),outlines)					
			unused.remove(s)
	
	#for k,v in outlines.items():
	#	for line in v:
	#		print(k,line)
	
	remove_excess_gaps(outlines)
	
	for k,v in outlines.items():
		for line in v:
			#print(k,line)
			if not "GAP" in line:
				x=line.strip().split()
				new=x[2]
				if new not in newold:
					newold[new]=k	#new=eg RXL_1, k = eg 'Scaffoldl_1'
	#newold					
	#{'RXL_1': 'Scaffold_1', 'RXL_2': 'Scaffold_2', 'RXL_3': 'Scaffold_3'}

				
	return unfilled_agp, outlines, newold
		

#get haplines and removes them from output	
def get_haps(unambigagp, tpfchunktags, tpfchunks, output, cleanedranges, scaffoldlengths, texel):

	scaff_dict={}
	hapcomponents=[]
	haptpfchunklens={}
	final={}
	
	for k,v in tpfchunktags.items():
		for tag in v[0]:
			if "HAPLOTIG" in tag:
				#print(k,tag)
				#if this tpfchunk results from a break...
				if k in cleanedranges:
					length=cleanedranges[k][1]-cleanedranges[k][0]+1
					haptpfchunklens[k]=length
					for line in tpfchunks[k]:
						if not "GAP" in line:
							hapcomponents.append(line.split()[1])	#add component
				#if the tpf scaffold is unbroken...
				else:
					scaff=k.split("%")[0]
					length=scaffoldlengths[scaff]
					haptpfchunklens[k]=length
					for line in tpfchunks[k]:
						if not "GAP" in line:
							hapcomponents.append(line.split()[1])	#add component
	
	lengthcheck=[]			
	for k,v in haptpfchunklens.items():
		if v>500*texel:
			if k in cleanedranges:
				start=str(cleanedranges[k][0])
				end=str(cleanedranges[k][1])
				vals=[k.split("%")[0]+"   span: "+start+"-"+end+"   length:   "+str(commas(v))+" bp"]
				lengthcheck.append(vals)
			else:
				print("\nLOOKS LIKE BREAK MISSING - IS NEEDED FOR HAPLOTIG DESIGNATION:\t",k)
	if len(lengthcheck)>0:
		warnings["hap_unusually_large"]=lengthcheck

	prev=""
	chrm=""
	for k,v in output.items():
		for line in v:
			if line.split()[1] not in hapcomponents:
				if not "GAP" in line:
					chrm=line.split()[2]
					#print(chrm,line)
					append_dict(chrm,line,scaff_dict)
					prev=chrm
				else:
					if chrm !="":
						append_dict(prev,line,scaff_dict)
						
	prev=""
	for k,v in scaff_dict.items():
		for line in v:
			if not "GAP" in line:
				chrm=line.strip().split()[2]
				append_dict(k,line,final)
				prev=chrm
			elif "GAP" not in final[k][-1]:
				append_dict(k,line,final) 
				
	#for k,v in final.items():
	#	for line in v:
	#		print(k,line)
			
	#for k,v in haptpfchunklens.items():
	#	print(k,v)				


			
	#scaff_dict is output minus haps
	#haptpfchunklens contains the hap tpfchunk key and its length so we can name haps by size
	return final, haptpfchunklens
					
	

	
def get_unlocs(unambigagp,tpfchunktags,tpfchunks,output,cleanedranges,scaffoldlengths, sex_chrms):

	scaff_dict={}
	unloctpfchunklens={}
	chrms_containing_unlocs=[]
	check={}
	
	for k,v in tpfchunktags.items():
		#print(k,v)
		for tag in v[0]:
			if "UNLOC" in tag:
				#print(k,tag)
				#if this tpfchunk results from a break...
				if k in cleanedranges:
					length=cleanedranges[k][1]-cleanedranges[k][0]+1
					unloctpfchunklens[k]=length
				#if the tpf scaffold is unbroken...
				else:
					scaff=k.split("%")[0]
					length=scaffoldlengths[scaff]
					total=sum(length)	#01/03/2022
					unloctpfchunklens[k]=total #01/03/2022
					#unloctpfchunklens[k]=length[0] #removed - not correct total length
					
	

							
	named_unlocs = name_unlocs(unloctpfchunklens,tpfchunktags,unambigagp,sex_chrms)

	#print(named_unlocs)


	
	return named_unlocs


#Here shrapnel is defined as anything in tpf but not agp (removing agp lines helps us define chrs as what we leave behind)	
def get_shrapnel(sscaffdict, tpfchunks):

	ascflist=[]
	tscflist=[]
	shrapnel=[]
	shrapchunks=[]
	for k,v in sscaffdict.items():
		for i in v:
			ascf=i[0]
			if ascf not in ascflist:
				ascflist.append(ascf)
	for k,v in tpfchunks.items():
		t=k.split("_")[0:2]
		tscf="_".join(t).split("%")[0]
		#print(k,tscf)
		if tscf not in tscflist:
			tscflist.append(tscf)
	for t in tscflist:
		if t not in ascflist:
			shrapnel.append(t)
			
			
	for k,v in tpfchunks.items():
		t=k.split("_")[0:2]
		tscf="_".join(t).split("%")[0]
		if tscf in shrapnel:
			shrapchunks.append(k)

			
	return shrapchunks
	


def remove_excess_gaps(outlines):



	#Fixing up excess GAP lines
	
	#Terminal gap lines
	for k,v in outlines.items():	
		#print(k,v)
		#removing gaps from start/end of chromosomes
		#if "GAP" in v:
		while "GAP" in v[0]:
			del v[0]
		while "GAP" in v[-1]:
			del v[-1]
			
	#internal gap duplicates
	lastline=""
	for k,v in outlines.items():
		for e,line in enumerate(v):
			if line==lastline and "GAP" in line:
				#print("here", line)
				lastline=line
				del v[e]
			else:
				lastline=line
							

					
	return outlines	
	


	

def name_unlocs(unloctpfchunklens,tpfchunktags,unambigagp,sex_chrms):


	sex2={}
	sizes_={}	#chrm key to size list in descending order of unlocs
	sizes={}
	named_unlocs={}
	iteration=0
	for k,v in sex_chrms.items():
		for i in v:
			if i not in sex2:
				sex2[i]=k
			
	for k,size in unloctpfchunklens.items():
		for a,b in unambigagp.items():
			if k == a.split("#")[1]:	 #Capture the precise tpfchunk
				#print(k,a)
				root=a.split(":")[0]
				if root not in sizes_:
					sizes_[root]=[size]
				else:
					sizes_[root].append(size)
					
	
	for k, v in sizes_.items():
		#print(k,v)
		sv=sorted(v)[::-1]
		sizes[k]=sv




	for k,v in sizes.items():
		#print(k,v)
		for a, b in unloctpfchunklens.items():
			#print(k,a,v)
			for i in v:
				if i==b:	#if it's the right size
					#print(k,sex2)
					if k not in sex2:
						#print("a",k,a,i,str(v.index(i)+1), k.replace("Scaffold",prefix)+"_unloc_"+str(v.index(i)+1))
						named_unlocs[a]=k.replace("Scaffold_",prefix)+"_unloc_"+str(v.index(i)+1)
					else:
						#print("b",k,a,i,str(v.index(i)+1), sex2[k]+"unloc_"+str(v.index(i)+1))
						named_unlocs[a]=sex2[k]+"_unloc_"+str(v.index(i)+1)
				#else:
					#print(k,a,i,b)
				
						
	

	#{'scaffold_49_ctg1%1': 'Z_unloc_1', 'scaffold_54_ctg1%1': 'Z_unloc_2', 'scaffold_57_ctg1%1': 'Z_unloc_3'}		
	return named_unlocs
	
	
#all the components that belong to unloc scaffolds 	
def get_unloc_comps(named_unlocs,tpfchunks):

	results={}
	for k,v in named_unlocs.items():
		#print(tpfchunks[k],v)
		for line in tpfchunks[k]:
			comp=line.split()[1]
			append_dict(comp,v,results)
			
	return results
	




def update_output_unlocs(unloc_comps,output,sex_chrms,texel):



	results={}
	check={}
	prob_comps={}
	final={}
	missed={}
	missed_filtered=[]
	
	
	for k,v in output.items():
		for line in v:
			#print(k,line)
			if not "GAP" in line:
				x=line.strip().split()
				comp=x[1]
		
				if comp in unloc_comps:
					#print(line,unloc_comps[comp])
					newline=line.replace(x[2],unloc_comps[comp][0])
					append_dict(k,newline,results)
				else:
					append_dict(k,line,results)
			else:
				append_dict(k,line,results)
				

				
	#Making a check on results to ensure no internal unlocs
	for k,v in results.items():
		for line in v:
			#print(k,line)
			if not "GAP" in line:
				append_dict(k,line,check)
				

				

	#if the only scaffold in a chromosome with unlocs has major TPF/AGP discrep(s) it will break the program - all we can do is report the offending chr in question and exit
	unloc_count=0
	for k,v in check.items():
		matching = [s for s in v if "unloc" in s]
		if len(matching)==len(v):
			print("\nSignificant TPF/AGP discrepancy somewhere in AGP Scaffold_"+k.replace("RL_","")+" (ie, chromosome "+k.replace("RL_","")+")")
			print("\nThe problem is in the main chromosome, not the unlocs")
			print("\nPlease check the manually assigned breaks carefully for this scaffold\n")
			sys.exit()
			
	#Remove legitimate terminal unlocs...
	for k,v in check.items():
		while "unloc" in v[0]:
			del v[0]
		while "unloc" in v[-1]:
			del v[-1]


				
	#...then look to see if any remain
	for k,v in check.items():	
		#if "unloc" in v:
		for line in v:
			if "unloc" in line:
				#print(k,line)
				comp=line.split()[1]
				append_dict(k,comp,prob_comps)
	
	iteration=0
	for k,v in check.items():	
		for line in v:
			comp=line.split()[1]
			if "unloc" not in line and k in prob_comps:
				name=k+"."+str(iteration)
				l=int(comp.split(":")[1].split("-")[0])
				h=int(comp.split(":")[1].split("-")[1])
				lo=min([l,h])
				hi=max([l,h])
				length=hi-lo
				#print(line)
				if name not in missed:
					missed[name]=[[length],[comp]]
				else:
					missed[name][0][0]+=length
					missed[name][1].append(comp)
			else:
				iteration+=1
				
				
	for k,v in missed.items():
		#Filter so we only report the true problems not the main chromosome
		if v[0][0] < 70*texel:
			missed_filtered.append(v[1])
				
	#if len(prob_comps)>0:
	for k,v in prob_comps.items():
		while len(v)>2:
			del v[1]
			


	#print(prob_comps)
	for k,v in output.items():
		if k in prob_comps:
			if len(prob_comps[k])>0:
				msg=[k,prob_comps[k],missed_filtered]
				append_dict("internal_unloc",msg,errors)
				
				
	#Ensure all unlocs get put in separate scaffolds (up to this point, they're joined into non-unloc scaffolds)
	for k,v in results.items():
		for line in v:
			#print(k,line)
			#scaff=""
			x=line.split()
			if not "GAP" in line:
				scaff=x[2]
				append_dict(scaff,line,final)
				#print(k,scaff,line)
				prev=scaff
			else:
				append_dict(prev,line,final)

				
		

	return final




	
	
	
def write_messages_to_screen(cleanedranges,sexchrm,painted_scaffs):

	

	if len(errors.keys())>0:
		print("\n"+borderlen*"=")
		print("\n\t\t\t\tERRORS\n")
		for k,v in errors.items():
			if k=="2nd_match":
				if len(set(v[1:]))==1:
					print(v[0])
				else:
					print("These occur only once in TPF - Either break missing or pretext editing error\n")
				for i in set(v[1:]):
					print(i+"  occurs"+str(v.count(i)+1)+"  times in AGP")
			if k=="hap_unloc":
				for i in v:
					print(i)
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")
				sys.exit()
			if k=="multiple_sex":
				for i in v:
					print(i)
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")
				sys.exit()
			if k=="multiple_sex2":
				for i in v:
					print(i)
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")
				sys.exit()
			if k=="internal_unloc":
				#Keep this - although it seems redundant it will catch ALL problems
				print("\nUnloc scaffs internal to painted chromosomes.  Have you missed tagging an unloc near...?\n")
				for i in v:
					if len(i[1])==1:
						print(i[0],i[1][0])
					else:
						print(i[0],i[1][0]+" to "+i[1][1])
				#....whereas this is filtered to try to catch just the useful information and could miss something
				if  len(v[0][2]) > 0:	#These are small missed unlocs that have passed our filter of n*texels
					print("\nLocations where 'unloc' seems to be missing from AGP:\n")	
					for i in v[0][2]:
						print(i)
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")
				sys.exit()	
			if k=="hetero_sex_haplo":
				print("Doesn't make sense - haplotig painted into heterogametic sex chromosome:\n")
				for i in v:
					print(i[0])
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")	
				sys.exit()	
			if k=="badsex":
				for i in v:
					print(i)
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")	
				sys.exit()	
			if k=="broken_shrapnel":
				print("Scaffold has been broken by curator but child scaffolds are not present in AGP\n")
				for i in v:
					print(i)
				print("\n\t\t>>> PLEASE FIX AGP/BREAKS AND RERUN <<<\n")
				sys.exit()
	
					
							
			else:
				if len(v)>1:
					for i in v:
						print(i)
					print("\n\t\t>>> PLEASE FIX BREAKS FILE AND RERUN <<<\n")
					sys.exit()	###IMPORTANT - turn this back on once program is written	
				else:
					print(v[0])
					print("\n\t\t>>> PLEASE FIX BREAKS FILE AND RERUN <<<\n")
					sys.exit()	###IMPORTANT - turn this back on once program is written


	outlist=[]

	
	
	
	for k,v in warnings.items():
		if k=="missing_breaks":
			outlist.append("\nBreaks missing from TPF (AGP coords with no corresponding TPF match):\n")
			#for a,b in v[0].items():
			outlist.append(v)
		if k=="ambig":
			outlist.append(v[0][0])
			outlist.append("\nProblem TPF coords:\n")
			for i in v[1]:
				outlist.append(i)
			outlist.append("\nAGP region:\n")
			for i in v[2]:
				if i not in outlist:
					outlist.append(i)


		if k=="TPF_problem_locations":
			outlist.append(v[0])
			for i in v[1:]:
				outlist.append("\t".join(i))
		if k=="accidentally_unpainted":
			outlist.append("\nYou may have accidentally not painted the following?\n")
			for i in v:
				line_num=str(i[0]+1)
				scaff=i[1].split("%")[0]
				inf="AGP line",line_num,"->",scaff
				outlist.append(" ".join(inf))	
					

		msgd="\n(if diff is small no action may be needed for the below)"			
		
		if k=="missing_breaks2":
			if msgd not in outlist:
				outlist.append(msgd)
			outlist.append(v[0])
			for i in v[1:]:
				outlist.append("\t".join(i))	

		if k=="excess_breaks":
			if msgd not in outlist:
				outlist.append(msgd)
			outlist.append(v[0])
			for i in v[1:]:
				outlist.append("\t".join(i))
				
		if k=="hap_unusually_large":
			outlist.append("\nUnusually large haplotig - is this correct?\n")
			for i in v:
				outlist.append(i[0])

	if len(outlist)!=0:
		print("\n\n"+borderlen*"=")
		print("\n\t\t\t\tWARNINGS\n")
		for o in outlist:
			print(o)

		print("\n")
		if len(warnings["AMBIG"])>0:
			print("AMBIGs represent TPF chunks that can't be placed to the AGP line, could be due to\nmissing or excess breaks somewhere in the scaffold:\n")
			lo=str(commas(cleanedranges[warnings["AMBIG"][0]][0]))
			hi=str(commas(cleanedranges[warnings["AMBIG"][0]][1]))
			print("AMBIG "+warnings["AMBIG"][0].split("%")[0]+" is TPF chunk number "+warnings["AMBIG"][0].split("%")[1]+" - "+lo+"-"+hi+"\n")
		if "no_breaks" in warnings:
			print(warnings["no_breaks"][0])


	
	
	print("\n"+borderlen*"=")

	print("\n\t\t\t\tINFO\n")


	if len(sexchrm.keys())==1:
		print("Sex chromosome:\t\t\t"+"".join(sexchrm.keys()))
	elif len(sexchrm.keys())==0:
		print("No sex chromosomes defined")		
	else:
		print("Sex chromosomes:\t\t"+",".join(sorted(sexchrm.keys())))
	if len(painted_scaffs.keys())>0:	
		print("Autosomes:\t\t\t"+str(len(painted_scaffs.keys())-len(sexchrm.keys()))+"\n")
	else:
		print("\n>>> !! NO PAINTED SCAFFOLDS !! <<<\n")	
	print("\t".join(info["masterlist"]))
	print("\t".join(info["ambig_tpf_to_agp"]))
	print("\t".join(info["shrapnel"]))
	print("\t".join(info["ambig_plus_unambig"]))
	if "ambig_plus_unambig_plus_shrapnel1" in info:
		print("\t".join(info["ambig_plus_unambig_plus_shrapnel1"]))
	if "ambig_plus_unambig_plus_shrapnel2" in info:
		print("\t".join(info["ambig_plus_unambig_plus_shrapnel2"]))	
	if "all_components_successfully_used" in info:
		print(info["all_components_successfully_used"][0])



	#This function will collate final information and write it to screen


	
def name_haps(haptpfchunklens):


	sizes=[]
	named_haps={}
	iteration=0
	for k,size in haptpfchunklens.items():
		if size not in sizes:
			sizes.append(size)
	sortedsizes=sorted(sizes, reverse=True)
	for size in sortedsizes:
		for k,v in haptpfchunklens.items():
			if size==v:
				iteration+=1
				name="H_"+str(iteration)
				named_haps[k]=name			
			
			
	return named_haps




def prepare_haps_tpf(named_haps, tpfchunks):

	outlines={}


	for k,v in named_haps.items():
		for line in tpfchunks[k]:
			if not "GAP" in line:
				x=line.strip().split()
				x[2]=v
				newline="\t".join(x)
				append_dict(v,newline,outlines)
			else:
				append_dict(v,line.strip(),outlines)
			
	return outlines
	
	
			
def write_hap_tpf(hapoutlines,tpfout):

	haptpf="haps_"+tpfout
	with open(haptpf,'w') as fout:
		for k,v in hapoutlines.items():
			for line in v:
				fout.write(line+"\n")

	
	
	
#Ensuring that all members of a chr are labelled with the sex even if all members not labelled in AGP				
def sex_components(tpfchunks,tpfchunktags,unambigagp,sex_chrms):

	comp_sex={}

	for k,v in unambigagp.items():
		chrm=k.split(":")[0]
		tpc=k.split("#")[1]
		for a,b in sex_chrms.items():
			if chrm==b[0]:
				for line in tpfchunks[tpc]:
					comp=line.split()[1]
					comp_sex[comp]=a
								
	return comp_sex
	
	
def apply_sex(hapoutlines,comp_sex,sex_chrms):
	
	sex_prtxt={}
	for k,v in sex_chrms.items():
		s=prefix+v[0].split("_")[1]
		sex_prtxt[s]=k
	#print(sex_prtxt)
	
	outlines={}
	for k,v in hapoutlines.items():
		for line in v:
			x=line.strip().split()
			scaff=x[2]
			if scaff in sex_prtxt:
				x[2]=sex_prtxt[scaff]
				newline="\t".join(x)
				append_dict(k,newline,outlines)
			else:
				append_dict(k,line,outlines)
				

				

	return outlines	
	
	
#remove original painted scaff key and replace keys with new keys derived from the new scaffold name (ie unloc keys)	
def update_chr_keys(outlines):

	updated={}
	chrm=""
	for k,v in outlines.items():
		for line in v:
			if not "GAP" in line:
				newkey=line.split()[2]
				chrm=newkey
				append_dict(newkey,line.strip(),updated)
			else:
				append_dict(chrm,line.strip(),updated)
				
	remove_excess_gaps(updated)

			
	return updated




		
#checking all components we think should be in output are in output - belt and braces
#This check takes  longer than any other part of the program unfortunately, but is very important
#and does prevent things going wrong for definite
def check_tpfout_sanity(output,tpfchunks,tpfout,haptpfchunklens):

	old_check=0
	check=[]
	check2=[]
	dups=[]
	probslog={}
	
	for k,v in tpfchunks.items():
		for line in v:
			if not "GAP" in line:
				cmpnt=line.strip().split()[1]
				if cmpnt not in check:
					check.append(cmpnt)
	old_check=len(check)

	

	
		

	for a,lines in output.items():
		for line in lines:
			if not "GAP" in line:
				cmpnt=line.strip().split()[1]
				if cmpnt not in check2:
					#print(a,cmpnt)
					check2.append(cmpnt)
					check.remove(cmpnt)
				else:
					dups.append(cmpnt)
			else:
				pass
				
	#Don't count hap components we know we're writing to haps as missing components		
	for k in haptpfchunklens.keys():
		for line in tpfchunks[k]:
			if not "GAP" in line:
				cmpnt=line.strip().split()[1]
				if cmpnt in check:
					check.remove(cmpnt)
	
	
				
	#This code was taking ages until I reduced the number of loops performed on tpfchunks to only those involving duplicate lines
	for k,v in tpfchunks.items():
		for line in v:
			for dup in dups:
				if dup in line:
					if k not in probslog:
						probslog[k]=v
					

	if len(probslog)>0:
		print("\nTPFchunk(s) used twice in agp - missing break(s) from breaks file?")
		for k,v in probslog.items():
			print("\n"+k)
			for i in v:
				print(i)	#We collect and report all the cases found at once
		print("\nPlease fix the problem and rerun\n")
		sys.exit()		#...and importantly then exit as we don't want to write a corrupted tpf			
	
	#print(check)
	#print(len(check))
	if len(check)==0:
		msg="\nAll "+str(old_check)+" tpf components used in outfile:   "+tpfout
		append_dict("all_components_successfully_used",msg,info)
	else:
		msg="\n\nERROR "+str(len(check))+" tpf components not used in output file!!\n"
		print(msg)
		print(check)
		print(msg)
						
								
	
#As there are no "GAP" lines in original.tpf, any "GAP" lines in output tpf must be joins - we simply count them
def count_joins(unambigagp,tpfchunktags):

	results={}
	joins=0
	unjoinedtchunks=[]
	
	for k,v in tpfchunktags.items():
		for t in v:
			if "UNLOC" in t or "HAPLOTIG" in t:
				unjoinedtchunks.append(k)
	
	for k in unambigagp.keys():
		sscaff=k.split(":")[0]
		tchunk=k.split("#")[1]
		if tchunk not in unjoinedtchunks:
			append_dict(sscaff,tchunk,results)
		
	for v in results.values():
		joins+=len(v)-1
		
		
		
	return joins	

						

#If it's in the breaks file, it should be in the AGP - error if this isn't true
def broken_but_not_in_agp(breaksnew,agpdict):


	used=[]
	for k,v in agpdict.items():
		for i in v:
			tscaff=i[0]
			if tscaff not in used:
				used.append(tscaff)
	for k in breaksnew.keys():
		if k not in used:
			append_dict("broken_shrapnel",k,errors)
			



def get_unpainted_comps(finalout2,tpfchunktags,tpfchunks):
	
	unpainted=[]
	unpainted_comps={}
	for k,v in tpfchunktags.items():
		#print(k,v[0])
		if "PAINTED" not in v[0]:
			unpainted.append(k)
			#print(k,v[0])
			
	#Some unpainted	are not necessarily in tagdict - we now add those too:	
	for k in tpfchunks.keys():
		if k not in tpfchunktags:
			unpainted.append(k)
	
	for u in unpainted:
		for line in tpfchunks[u]:
			if not "GAP" in line:
				x=line.split()
				comp=x[1]
				append_dict(u,comp,unpainted_comps)
			
	#unpainted_comps is tpfchunkkey to component list, no gap lines
	#for k,v in unpainted_comps.items():
		#print(k,v)
			


	return unpainted_comps
		

	
	
	
#Only an issue if flood fill hasn't been used, but still, let's check...	
def check_for_accidentally_unpainted(tpfchunktags,unambigagp,agplines):

	check=[]
	check2={}
	check3=[]
	for a,b in agplines.items():
		for k,v in unambigagp.items():
			tkey=k.split("#")[1]
			akey=k.split("#")[0]
			if a==akey:
				#print(a,tkey)
				if tkey in tpfchunktags and tpfchunktags[tkey][0].count("PAINTED")==1:	#must be painted
					check.append([tkey,"PAINTED"])
				else:
					#print(tkey)
					check.append([tkey,"NOT"])
		

	#Assume anything left unpainted at r/h end of map is intentional 
	while check[-1][1]=="NOT":
		del[check[-1]]
	
	#for k in check:
	#	print(k)
	

	for i,c in enumerate(check):
		if c[1]=="NOT":
			check3.append([i,c[0]])
			
	if len(check3)>0:
		warnings["accidentally_unpainted"]=check3



#Here we are only dealing with chromosomal entities
def write_chr_csv(finalout3):

	temp={}
	temp2={}
	temp3={}
	results={}
	#Make dict of chr to unloc
	for k,v in finalout3.items():
		for line in v:
			#print(k,line)
			if not "GAP" in line:
				#print(k,line)
				x=line.split()
				chrm=x[2]
				base="_".join(chrm.split("_")[0:2]).replace("_unloc","")	#eg: Y,RL_6,scaffold_72
				#filtering for non-chromosomal entities
				if "affol" not in chrm and "AMBIG" not in chrm and "tg" not in chrm and "tig" not in chrm:
					if base not in temp:
						temp[base]=[chrm]
					elif chrm not in temp[base]:
						temp[base].append(chrm)
						
					
	
						
						
	for k,v in temp.items():
		for i in v:
			#print(k,i)
			x=i.split("_")
			#print(x)
			#get the chromosomes
			if len(x)<3:	#chrs excluding unlocs, eg 'RL', '12'
				results[k]=[i]
			else:
				append_dict(k,i,temp2)
				
	#for k,v in results.items():
	#	print(k,v)


	for k,v in temp2.items():
		iterations=[]
		for i in v:
			#print(k,i)
			x=int(i.split("_")[-1])
			iterations.append(x)
			sorts=sorted(iterations)
			temp3[k]=sorts

	#chr keys with a list of unloc_iteration
	for k,v in temp3.items():
		for s in v:
			#print(k,s)
			#print(k,v)
			new=k+"_unloc_"+str(s)
			if new not in results[k]:
				results[k].append(new)
				
	with open("chrs.csv","w") as fout:
		for k,v in results.items():
			line=",".join(v)
			fout.write(line+"\n")








#Takes RX_ prefixed scaffolds and reverts names to scaffold_ names for non-painted entities
def reset_names_unpainted(finalout2,unpainted_comps,painted_scaffs,newold):

	results={}
	painted={}
	unpainted={}
	ambigs={}
	check={}
	
	
	for k,v in finalout2.items():
		if k[0] in sex:
			for line in v:
				append_dict(k,line,painted)	#These are defined sex chromosomes and are therefore painted - we can keep line  unchanged
		else:
			if len(k.split("_"))>1:
				agpscaff="Scaffold_"+k.split("_")[1]
				if agpscaff in painted_scaffs and "AMBIG" not in k:
					for line in v:
						append_dict(k,line,painted)	#These are painted chrs - we keep line unchanged
				elif "AMBIG" in k:
					for line in v:
						append_dict(k,line,ambigs)	#These are ambigs - we keep line unchanged		
				else:
					for line in v:
						append_dict(k,line, unpainted)	#These are unpainted - we will reset scaff names below
			else:	#unlocs have to be in painted chrms
				for line in v:
					append_dict(k,line,painted)		

						
	#Now we need to reset our unpainted names and ensure they're unique
	#Unpainted  entities will only occur on a single AGP line so we can take a base name from the component and just add an iterator if we have multiples
	for k,v in unpainted.items():
		#print(k,v[0])
		#we only need consider first line in our line list since all scaffolds with be the same in unpainted line list
		base=v[0].split()[1].split(":")[0]	#k = RXL_698, base = scaffold_281
		append_dict(base,k,check)		#Now we can get a unique iterator for our base names from their index in check
		
		
	
		
	#Now we can create our final table with reset and guaranteed unique names
	#checking here for scaffs that exist in chrs but that also exist in shrapnel
	for k,v in finalout2.items():
		if k in painted:
			for line in v:
				if not "GAP" in line:
					x=line.strip().split()
					base=line.split()[1].split(":")[0]
					if base in check:
						check[base].append(k)
		elif k in ambigs:
			for line in v:
				if not "GAP" in line:
					x=line.strip().split()
					base=line.split()[1].split(":")[0]
					if base in check:
						check[base].append(k)
						
	
	#checking here for scaffs that exist duplicate in shrapnel only
	for k,v in finalout2.items():
		if k in painted:
			for line in v:
				append_dict(k,line,results)
		elif k in ambigs:
			for line in v:
				append_dict(k,line,results)
		else:
			for line in v:
				if not "GAP" in line:
					x=line.strip().split()
					base=line.split()[1].split(":")[0]
					if len(check[base])==1:				#unique scaffolds - just add them
						x[2]=base
						newline="\t".join(x)
						append_dict(k,newline,results)
					else:
						#print(check[base])
						new=base+"_"+str(check[base].index(k)+1)	#add an iterator if base is not unique to that agp scaff
						x[2]=new
						newline="\t".join(x)
						#print(k,newline)
						append_dict(k,newline,results)
				else:
					append_dict(k,line,results)
					

				
	return results





def check_tagdict(tagdict):



	unpainted_sex=[]

	for k,v in tagdict.items():
		for tlist in v:
			for tag in tlist:
				if tag[0] in sex and "PAINTED" not in tlist:
					unq=tag+"@"+k.replace("/","\t")
					unpainted_sex.append(unq)
				

					

	if len(unpainted_sex)>0:
		print("Sex chromosome not painted - please paint\n")
		for i in unpainted_sex:
			print(i)
		print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")	
		sys.exit()	
		


#Getting any AGP scaff > 10 texels on the basis that this is too big to be an artefact
def count_minimum_expected_agp_scaffs(agp,texel):


	results={}
	with open(agp,'r') as f:
		for line in f:
			if line[0] != "#":
				x=line.strip().split()
				if x != [] and x[4] != "U":	#is not a gap
					ascaff=x[5]
					alo=min(int(x[6]),int(x[7]))
					ahi=max(int(x[6]),int(x[7]))
					lngth=ahi-alo
					meta=[]
					for i in range(len(x[9::])):
						if x[9::][i] not in meta:
							meta.append(x[9::][i].upper())
					if lngth > 10*texel:	
						if meta.count("PAINTED")>0:
							p="PAINTED"
						else:
							p="NOT"
						vals=[str(alo),str(ahi),p]
						append_dict(ascaff,vals,results)

	#scaffold_11 [['3497070', '6680030', 'PAINTED'], ['2533805', '3497069', 'PAINTED'], ['1', '2533804', 'NOT']]
	return results




def count_output_chromosomes(finalout3, painted_scaffs, sexchrm):
	
	outchrs=[]
	for k,v in finalout3.items():
		for line in v:
			if not "GAP" in line:
				scf=line.split()[2]
				if not "ffol" in scf:
					if scf not in outchrs:
						outchrs.append(scf)
	

	sx=len(sexchrm.keys())
	total=len(painted_scaffs.keys())
	autos=total-sx
	actual=len(outchrs)
	if total > actual:
		print("\n\n>>>>>  WARNING - fewer output chromosomes than expected!! <<<<<\n\n")
		print("\tPainted  sex chrs:\t"+str(sx))
		print("\tPainted  autosomes:\t"+str(autos))
		print("\tTotal chrs in output:\t"+str(actual)+"\n")
		print("\n>>>>>  Accidental  hap tag on whole chromosome? <<<<<\n\n")
				






		


def main():



	parser = argparse.ArgumentParser(description='Designed to take pretext generated AGP and fit your assembly TPF to it.') 

	#positional args
	parser.add_argument('tpf', metavar='tpf', type=str, help='assembly TPF with gaps that are needed to allow rearrangement to match the edited PretextView map annotated with a ">" as the first character of the GAP line.')
	parser.add_argument('agp', metavar='agp', type=str, help='Pretext agp')
	#parser.add_argument('breaks', metavar='breaks', type=str, help='breaks file')
	#parser.add_argument('fasta', metavar='fasta', type=str, help='original assembly fasta')
	


	#display help when misusage
	if len(sys.argv) <2: 
		parser.print_help()
		

	args = parser.parse_args()  #gets the arguments
	start_time = datetime.now()


	print("\n")


	check_input_tpf(args.tpf)

	tpfdict, errors2=parse_tpf(args.tpf)

	
	compare_scaff(tpfdict,args.agp)

	ctg_lengths=contig_lens(tpfdict)


	gsize,texel=genome_size(tpfdict)
	scafflens=lens(tpfdict)
	
	
	fragcutoff=1*texel
	agpdict,discards,agplines,tagdict,sexchrm,painted_scaffs, pv_vers, pv_res = scaffs_from_agp(args.agp,fragcutoff,scafflens,ctg_lengths,texel)	#All the agp order and orientation information

	#where pv_res hasn't been set due to old version of pretextView without info in header
	if pv_res==0:
		pv_res=gsize/32767



	check_tagdict(tagdict)


	report_errors(errors2)
	checkin=components_from_dict(tpfdict)

	dividers=agp_dividers(agpdict)
	breaksnew=get_breaks(args.tpf)	#Unique to the XL version - these are curator directed tpf break coords
	
	broken_but_not_in_agp(breaksnew,agpdict)

	tpfchunks=breaktpf(tpfdict,breaksnew)

	unbroken=get_unbroken_scaff(tpfchunks)

	tpfchunk_ranges=ranges(tpfchunks)

	
	shrapnel=get_shrapnel(agpdict, tpfchunks)
	#print(shrapnel)
	

	tpfout="rapid_prtxt_XL.tpf"
	
	ambigtcs, unambigtcs, unambigagp, agplines, ulist = maptpf2agp(tpfchunk_ranges,tpfchunks,agpdict,prefix,unbroken,tpfdict,texel,scafflens,shrapnel,args.tpf,tpfout)



	tpfchunktags = tag_tpfchunks(tagdict,unambigagp)
	
	minagpscaffcount = count_minimum_expected_agp_scaffs(args.agp,texel)
	
	painted=False
	for k,v in tpfchunktags.items():
		if "PAINTED" in v[0]:
			painted=True
	
	if painted==True:
		check_for_accidentally_unpainted(tpfchunktags,unambigagp,agplines)	#only an issue if flood fill hasn't been used in pretext
	
	tpf_vs_agp_discreps = break_count_vs_agp_line_count(unambigagp,tpfchunks,minagpscaffcount)
	locate_problems(tpf_vs_agp_discreps, unambigagp,tpfchunks,tpfchunk_ranges,texel)
	
	additional_check1(unambigagp,tpfchunk_ranges,texel,agplines,tpfchunks,minagpscaffcount, pv_res)	#Reports missing and excess breaks 
	
	
	unfilled,output,newold_superscaffs = prepare_output_tpf(agplines,unambigagp,tpfchunks,tpfout,shrapnel,ambigtcs)

	haplessoutput, haptpfchunklens = get_haps(unambigagp, tpfchunktags, tpfchunks,output, tpfchunk_ranges, scafflens, texel)
	
	
	
	
	named_unlocs = get_unlocs(unambigagp,tpfchunktags,tpfchunks,haplessoutput,tpfchunk_ranges,ctg_lengths, sexchrm)
	

	
	unloc_comps = get_unloc_comps(named_unlocs,tpfchunks)
	#print(unloc_comps)

	unlochaplessoutput = update_output_unlocs(unloc_comps, haplessoutput, sexchrm,texel)
	


	write_messages_to_screen(tpfchunk_ranges,sexchrm,painted_scaffs)	#Moved this line from below check_tpfout_sanity to make errors fail before output is written


	named_haps = name_haps(haptpfchunklens)

	hapoutlines = prepare_haps_tpf(named_haps, tpfchunks)
	
	comp_sex = sex_components(tpfchunks,tpfchunktags,unambigagp, sexchrm)
	sexedunlochaplessoutput = apply_sex(unlochaplessoutput,comp_sex, sexchrm)
	
	finalout1 = update_chr_keys(sexedunlochaplessoutput)
	
	finalout2 = remove_excess_gaps(finalout1)

	unpainted_comps = get_unpainted_comps(finalout2,tpfchunktags,tpfchunks)		#NEW
	
	if painted==True:
		finalout3 = reset_names_unpainted(finalout2,unpainted_comps, painted_scaffs,newold_superscaffs)		#NEW
	else:
		finalout3=finalout2
	
	
	count_output_chromosomes(finalout3, painted_scaffs, sexchrm)
	
	
	check_tpfout_sanity(finalout3,tpfchunks,tpfout,haptpfchunklens)
	write_output_tpf(finalout3,tpfout)
	
	write_chr_csv(finalout3)		#NEW
	
	#if len(hapoutlines)>0:
	write_hap_tpf(hapoutlines,tpfout)
	

		
	gsize2=genome_size2(finalout2)	#Haps removed
	gsize3=genome_size2(hapoutlines)	#just haps

	joins = count_joins(unambigagp,tpfchunktags)
	
	breakscount=0
	for i in breaksnew.values():
		for val in i:
			breakscount+=1
	
	print("\nAGP from PretextView Version\t"+pv_vers)
	print("Texels in map:\t\t\t"+str(commas(int(gsize/pv_res))))
	print("Texel length:\t\t\t"+str(commas(int(pv_res)))+"bp")
	print("Input genome size:\t\t"+str(commas(gsize))+"bp")
	print("Output genome size:\t\t"+str(commas(gsize2+gsize3))+"bp\t\t(primary + haps)")
	print("Break count:\t\t\t"+str(commas(breakscount)))
	print("Join count:\t\t\t"+str(commas(joins)))
	print("Haps count:\t\t\t"+str(commas(len(named_haps))))
	print("Unlocs count:\t\t\t"+str(commas(len(named_unlocs))))
	breaks_and_joins=joins+breakscount
	changes=round(breaks_and_joins/(gsize2/1000000000),2)
	if changes==0:
		print("No manual interventions!\n")
	elif changes >=1 and changes <= 49:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t\t(very low)\n")
	elif changes >=30 and changes <= 49:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t\t(below average)\n")
	elif changes >=50 and changes <= 199:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t\t(around average)\n")
	elif changes >=200 and changes <= 349:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t\t(above average)\n")
	elif changes >=350 and changes <= 399:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t\t(moderately high)\n")
	elif changes >=400 and changes <= 549:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t\t(high)\n")
	elif changes >=550 and changes <= 949:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t\t(very high)\n")
	elif changes >=950:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t\t(unusually high)\n")			
	print(borderlen*"=")
	print("\nChecking "+tpfout+" sanity...")
	tpf_sanity(tpfout)
	if len(hapoutlines)>0:
		print("\nChecking "+"haps_"+tpfout+" sanity...")
	print("\nWritten:\n")
	print(location(tpfout))
	print(location("haps_"+tpfout))
	print(location("chrs.csv"))		#NEW
	end_time=datetime.now()
	print('\n\nFINISHED:\t{}'.format(end_time - start_time)+"\n\n")
	
	



if __name__ == '__main__':
	main()
