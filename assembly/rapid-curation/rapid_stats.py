#!/usr/bin/env python
#Alan Tracey Last updated 16/07/2021

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
import argparse
from datetime import datetime
import pyfastaq
import re

# get scaffold ids
def get(fl):
    res={}
    with open(fl,'r') as f:
        for line in f:
            x=line.split()
            if not "GAP" in line:
                base,coords=x[1].split(":")
                if not base in res:
                    res[base]=[]
                res[base].append(coords)
    return res  
    

# merge 2 hashes of lists
def combine(cur,haps):

    for k,v in haps.items():
        if k not in cur:
            cur[k]=v
        else:
            for i in v:
                if i not in cur[k]:
                    cur[k].append(i)
    return cur      


# compare two hashes of lists, to get an approx number of breaks
def compare(curH,orig):
    count=0
    for k,v in curH.items():
        for i in v:
            if i not in orig[k]:
                count+=1    #Number of changed 'components' in curated tpf
    return round(count/2)
                
# count regular gaps in a TPF
def count_gaps(orig):
    count=0
    with open(orig,'r') as f:
        for line in f:
            if "GAP" in line:
                count+=1
    return count

# count all non-200bp gaps
def count_gaps2(cur):
    countother=0
    with open(cur,'r') as f:
        for line in f:
            x=line.split()
            if x[0]=="GAP" and x[2] !="200":
                countother+=1
    return countother   

# get pairs from the same parents
def get_pairs(file):
    pairs=[]
    singletons=[]
    last1=''
    with open(file,'r') as f:
        for line in f:
            # if the last has the same parent push a+b into the pairs
            if "GAP" in line:
                continue
            a=line.split()
            b=last1.split()
            if b and a[2] == b[2]:
                pairs.append(a[1]+b[1])
            last1=line
            singletons.append(a[1])
    return pairs,singletons

# grab gap lines from the original TPF and return them as concatenated triples
def get_gaps(orig):
    gaps=[]
    last2=''
    last1=''
    with open(orig,'r') as f:
        for line in f:
            if "GAP" in last1:
                gaps.append(transform_id_before(last2)+last1+transform_id_after(line))
                gaps.append(transform_id_after(line)+last1+transform_id_before(last2))
            last2=last1
            last1=line
    return gaps

# get the id and the coord nearest to the gap
def transform_id_before(line):
    x=line.split()
    name,start,end = re.split(':|-',x[1])
    if x[3] == "PLUS":
        return name+":"+end
    else:
        return name+":"+start

def transform_id_after(line):
    x=line.split()
    name,start,end = re.split(':|-',x[1])
    if x[3] == "PLUS":
        return name+":"+start
    else:
        return name+":"+end

# compare a TPF with a list of gap triples
def compare_gaps(tpf,gaps):
    last2=''
    last1=''
    both =0
    n    =0
    with open(tpf,'r') as f:
        for line in f:
            if "GAP" in last1:
                if  transform_id_before(last2)+last1+transform_id_after(line) in gaps:
                    both+=1
                else:
                    n+=1
                    if args.debug:
                        print(last2+last1+line+"--") 
            last2=last1
            last1=line
    return(both,n)


def count_haps(curH):

    count=0
    scaffs=[]
    with open(curH,'r') as f:
        for line in f:
            if not "GAP" in line:
                x=line.strip().split()
                if x != []:
                    scaff=x[2]
                    if scaff not in scaffs:
                        scaffs.append(scaff)
                        count+=1
    return count            


def count_joins(cur):
    
    with open(cur,'r') as f:
        count=0
        for line in f:
            if "GAP\tTYPE-2\t200" in line:
                count+=1
    joins="\nJoin count:\t"+str(count)
    return joins

def tpf_sanity(tpfout):

    cmd="perl /software/grit/projects/vgp_curation_scripts/test_tpf_sanity.pl -scafflevel "+tpfout
    pyfastaq.utils.syscall(cmd)
    
def fast_check(orig,cur,curH):

    v=orig.replace(".tpf","")+"_plus_"+curH
    cmd="cat "+orig+" "+curH+" > "+v
    pyfastaq.utils.syscall(cmd)
    cmd="python /nfs/team135/yy5/scripts/FastCheck.py -t1 "+orig+" -t2 "+v
    print("\nWhile we're here, I will run Yumi's script to check for accidental component deletion...\n\n"+cmd)
    print("\nIf no output, there were no problems found")
    pyfastaq.utils.syscall(cmd)

def main():

    parser = argparse.ArgumentParser(description='Calculates number of breaks, joins and haplotypic duplicate removals.') 

    #positional args
    parser.add_argument('t1', metavar='tpf1', type=str, help='original tpf')
    parser.add_argument('t2', metavar='tpf2', type=str, help='curated tpf')
    parser.add_argument('t3', metavar='tpf3', type=str, help='haplotigs tpf')

    
    #optional arguments                                                                                                                                                                                                                                              

    parser.add_argument('--rapid', type=str, help='tpf is generated by rapid_split.pl', default=True)
    parser.add_argument("-d","--debug", help="print debug strings", action="store_true")

    #display help when misusage
    if len(sys.argv) <3: 
        parser.print_help()
        
    global args
    args = parser.parse_args()  #gets the arguments
    start_time = datetime.now()

    print("\nChecking "+args.t1+" sanity...\n\n")
    tpf_sanity(args.t1)
    print("\nChecking "+args.t2+" sanity...\n\n")
    tpf_sanity(args.t2)
    print("\nChecking "+args.t3+" sanity...\n\n")
    tpf_sanity(args.t3)
    
    orig=get(args.t1)
    cur=get(args.t2)
    haps=get(args.t3)

    curH=combine(cur,haps)
    changedcomps=compare(curH,orig)
    origgapcount=count_gaps(args.t1) # all old
    curgapcount=count_gaps2(args.t2) # only non-200
    hapgapcount=count_gaps2(args.t3) # only non-200
    hapscount=count_haps(args.t3)
    joins=count_joins(args.t2)

# get new gaps from the curated files   
    gaps=get_gaps(args.t1)
    both,new = compare_gaps(args.t2,gaps)
    hap_both,hap_new = compare_gaps(args.t3,gaps)
    new_total=new+hap_new # combined new gaps

# get paired and single ids (minus gaps)
    curated_pairs,cur_s = get_pairs(args.t2)
    hap_pairs,hap_s = get_pairs(args.t3)
    orig_pairs,orig_s = get_pairs(args.t1)
 
# should be all changed components  
    changes=0
    for i in cur_s+hap_s:
        if i not in orig_s:
            changes+=1


    print("\n\n================================\n")
    
    fast_check(args.t1,args.t2,args.t3)
    if args.debug:
        print("\n-------------------------")
        print("   STATS:")
        print("-------------------------")
        # all original gaps - curated non-200bp gaps - curated non-200bp hap gaps + number of changes/2
        print("\nBreak count:\t"+str(origgapcount-curgapcount-hapgapcount+changedcomps))
        print(joins) # 200bp gaps in the curated 
        print("\nHaps count:\t"+str(hapscount))

    print("\n\n-------------------------")
    print("   STATS (M):")
    print("-------------------------")
    print(f"\nBreak count:\t{round(changes/2)+origgapcount-both}")
    print(f"\nJoin count:\t{new_total}")
    print(f"\nHaps count:\t{hapscount}")

    end_time=datetime.now()
    print("\n-------------------------")
    print('\n\nFINISHED:\t{}'.format(end_time - start_time)+"\n\n")


if __name__ == '__main__':
    main()
