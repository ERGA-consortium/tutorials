import pyfastx
import argparse



def read_file(tpfFile):

	tpfDic = { line.split("\t")[1].strip() : line.split("\t")[1].strip().split(".")[0] for line in open(tpfFile) }

	return tpfDic



parser = argparse.ArgumentParser()
   
parser.add_argument('-f1', type=str, help='first fasta file')
parser.add_argument('-f2', type=str, help='second fasta file')
parser.add_argument('-t1', type=str, help='pre-curation tpf file')
parser.add_argument('-t2', type=str, help='curated tpf file')


args = parser.parse_args()





if args.f1 and args.f2:

	file1 = pyfastx.Fasta(args.f1)
	file2 = pyfastx.Fasta(args.f2)

	print (int(file1.size) - int(file2.size))

if args.t1 and args.t2 :

	td2=read_file(args.t2)
	td1=read_file(args.t1)

	for key, value in td1.items():
		if value not in td2.values():
			print ( "component " + value + " is not in tpf 2")








