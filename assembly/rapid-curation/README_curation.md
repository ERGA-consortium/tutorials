# Rapid curation - Curation tools guide

This section will cover the curation tools and how they work.


## Support 
Slack Channel for curation assembly help:
```
https://join.slack.com/t/assemblycuration/shared_invite/zt-1kx2ww71y-823ruaAxswgQGypgofBaOA
```

PretextView tutorial:

PretextView - Tutorial.pdf
```
https://youtu.be/3IL2Q4f3k3I
https://youtu.be/LWy6pwCQNDU
```

PretextView can be obtained here
```
https://github.com/wtsi-hpag/PretextView/releases
```

Please also see documentation:

RAPID CURATION TRAINING MANUAL.pdf



Telomere identification script (telo_finder.py):

The telo_finder script aims to detect likely telomere motifs in assemblies.
It assumes the assemblies are high quality (hi-fi or similar) and that there are a number of chromosomes (it won't return a result with just one chromosome for example).
It uses a hard-coded set of non-redundant telomere motifs with which it attempts to distinguish true telomere sequences from background noise.


usage: telo_finder.py [-h] [--size size] [--klo klo] [--khi khi] [--ends ends]
                      fasta

Finds most likely telomere motif in Hi-fi or equivalent quality assembly where
telomeres are expected to occur at the ends of multiple scaffolds

positional arguments:
  fasta        assembly fasta

optional arguments:
  -h, --help   show this help message and exit
  --size size  top and tail scf bp (default: 200)
  --klo klo    min kmer (default: 4)
  --khi khi    max kmer (default: 15)
  --ends ends  ends to scan (default: 1000)
  

## Suggested workflow




Split assembly fasta to create a TPF (rapid_split.pl).


Manipulate the assembly in PretextView app - Action any necessary breaks in TPF *contigs* in order to mirror breaks in PretextView. 
This involves editing coordinates and introducing new gaps in the TPF as necessary.Breaks you want to be made require 
that a ">" is inserted as the first character of the gap line in the TPF e.g:

?	scaffold_1:2418601-2556927	scaffold_1	PLUS
GAP	TYPE-2	499
?	scaffold_1:2557427-3079467	scaffold_1	PLUS
>GAP	TYPE-2	23958
?	scaffold_1:3103426-3453292	scaffold_1	PLUS
GAP	TYPE-2	200224
?	scaffold_1:3653517-4134217	scaffold_1	PLUS

No TPF edits are necessary for joins.

Any hap-dups or unlocalised scaffolds should be tagged using meta-data
tag mode in PretextView. Unlocs (linked to chromosome but no specific location found) should be placed at the start/and/or/end 
of a chromosome to be included when the chromosomes are painted in the subsequent step.

Once the map has been rearranged and meta-data tags added paint the chromosomes.

Export the map to AGP. (Some input scaffolds won't feature in the output AGP due to PretextView resolution.  
These will be recovered in the next step).

Run rapid_pretext2tpf_XL.py, pointing it at the edited TPF (where all necessary gaps have been added) and the Pretext AGP. 
This will output a TPF that mirrors the chromosomes that have been built in PretextView, and which will have the same number 
of basepairs as the original TPF. Any errors or warnings should point you to inconsistencies between the TPF and AGP which need 
addressing.The script will produce stats (joins/breaks/hap-dup removals).

Turn the new TPF back into a fasta file using rapid_join.pl.

Produce a new Pretext map from the curated fasta file to check it.

A count of the number of breaks/joins and haplotypic duplicate removals can be produced using rapid_stats.p

An alternative to rapid_pretext2tpf_XL.py is rapid_pretext2tpf.py.  It is highly recommended that the XL version be used.
However, at the user's risk, the non-XL version can be used.  The advantage of this is that the curator doesn't have to
specify which gaps to break the TPF at.  It is possible if a region of high gap density exists on a scaffold that the wrong gap
will be chosen to break the scaffold.  Hence the XL version of the script is strongly recommended.





Scripts and documentation for Rapid curation:

## rapid_split.pl

This script takes a fasta file and produces a TPF on contig level, i.e. it splits at all Ns

Usage:

perl split.pl -fa <fasta>

              -printfa # if you want to print out the split fasta

              -h/help  # this message
	      

##  rapid_pretext2tpf.py	      

usage: rapid_pretext2tpf.py [-h] tpf agp

Designed to take pretext generated AGP and fit your assembly TPF to it.

positional arguments:

  tpf         assembly TPF with gaps as needed to allow rearrangement to match
              the edited PretextView map.

  agp         Pretext agp
  


##	rapid_pretext2tpf_XL.py 

usage: rapid_pretext2tpf_XL.py [-h] tpf agp

Designed to take pretext generated AGP and fit your assembly TPF to it.

positional arguments:
  tpf         assembly TPF with gaps that are needed to allow rearrangement to
              match the edited PretextView map annotated with a ">" as the
              first character of the GAP line.
  agp         Pretext agp

 

	     
## rapid_stats.py

Calculates number of breaks, joins and haplotypic duplicate removals.

positional arguments:
  tpf1           original tpf
  tpf2           curated tpf
  tpf3           haplotigs tpf

optional arguments:
  -h, --help     show this help message and exit
  --rapid RAPID  tpf is generated by rapid_split.pl
usage: rapid_stats.py [-h] [--rapid RAPID] tpf1 tpf2 tpf3

## test_tpf_sanity.pl

Name: test_tpf_sanity.pl

    A script to run some sanity test on a TPF file

Usage: perl tpf_sanity.pl [-scafflevel] FILE

    FILE is a file in TPF format. Either contig or scaffold level.

  Options:

    -scafflevel to use on scaffold level TPFs. Turns onchecks for naming,
    
                order and gaps

## rapid_join.pl

This script takes an original assembly fasta file, a one-line per chromosome pre-csv file and a TPF file generated with split.pl and creates the assembly from the TPF file

Usage:

perl join.pl -fa <fasta>
             -tpf <tpf>
             -csv <pre-csv>
Calculates number of breaks, joins and haplotypic duplicate removals.

positional arguments:

  tpf1           <original tpf>

  tpf2           <curated tpf>

  tpf3           <haplotigs tpf>

optional arguments:

  -h, --help     show this help message and exit

  --rapid RAPID  tpf is generated by rapid_split.pl

usage: rapid_stats.py [-h] [--rapid RAPID] tpf1 tpf2 tpf3

## tpf_cumulative.py

Produces a pseudo-AGP file from a TPF, giving the user cumulative coordinates
for their TPF

positional arguments:

  tpf         <curated tpf>

## dependencies

Dependencies for the python scripts:

rapid_stats.py:
sys, argparse, datetime, pyfastaq

rapid_pretext2tpf.py:
sys, argparse, datetime, pyfastaq, subprocess

tpf_cumulative.py:
sys, argparse, pyfastaq

