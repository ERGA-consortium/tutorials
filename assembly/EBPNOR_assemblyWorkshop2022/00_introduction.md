# Welcome to the Genome assembly, curation and validation workshop!

Today you'll learn how to assemble a whole genome, EBP-Nor style! 

After attending the workshop you should:
- know about the most-used approaches for genome assembly
- be able to assess information inherit in sequencing reads
- be able to validate genome assemblies
- know about manual curation of assemblies


## Our dataset

*Metschnikowia zobellii* is a yeast found in arctic coastal climates. It was first discovered as a water flea parasite, but have since been found in a wide range of organisms, including plants and other arthropods. *Metschnikowia zobellii* has a small genome, with only five chromosomes. You can read more about it [here](https://www.darwintreeoflife.org/news_item/9-species-our-scientists-are-excited-about-sequencing-in-2022/) (scroll down a bit).

[Darwin Tree of Life](https://www.darwintreeoflife.org) has sequenced the yeast, but not published a genome not for it yet. Fortunately for us, they allow everyone to play with the data anyhow, so we will do that. DToL has some interesting webpages where they list several quality measures for the sequencing (some of which we will do in this workshop) [here](https://tolqc.cog.sanger.ac.uk/darwin/fungi/Metschnikowia_zobellii/). It can be worth a look. We downloaded the data from [ENA](https://www.ebi.ac.uk/ena/browser/view/GCA_939531405.1) and subsampled it to get it to the coverages we expect/plan for. 

The genome itself is 14 Mbp, and the PacBio data was almost 17 Gbp, more than 1000x coverage. It was subsampled with [seqtk](https://github.com/lh3/seqtk) like this:
```
seqtk sample ERR9588940.fastq.gz 40000 |gzip > ERR9588940_30x.fastq.gz
```
to get about 30x coverage with PacBio data. For the exact data you are working with, we also added a surprise which you will discover later in the workshop.



## Why do we use a combination of HiFi and Hi-C reads? 

HiFi sequencing creates highly accurate circularized consensus reads. How are these reads generated? By ligating hairpin adapters, the DNA fragment that is being sequenced becomes a circle. This means that the machine can du multiple passes over the same DNA-sequence, to weed out any misread nucleotides. This is how HiFi reads can be so long, while remaining over 99.9% accurate. 

Hi-C sequencing is done to capture how the chromatin is folded within the cell nucleus. By ligating the folded DNA-strands, we can capture which loci are found in close proximity, and thus which parts of the DNA are found within the same chromosomes.

When combining these two, we can create haplotype resolved assemblies, meaning we can separate reads by maternal and paternal origin, without having access to parental data. In diploid, or polyploid organisms, this adds another level of information, and creates more accurate assemblies than a primary and alternate assembly would. 

## Package management

Administrating the different programs that are needed in project can be a hassle. We like conda, especially [miniconda](https://docs.conda.io/en/latest/miniconda.html), and have set up different environments we will use for the different analyses. [Bioconda](https://bioconda.github.io) contain a lot of different packages that are relevant for us, and genomics and bioinformatics in general.

To load conda, do this:
```
eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 
```

There are some of the different programs that are not available through conda. For most of these we use [Singularity containers](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html). 

We mostly set up scripts and arranged data so it is ready to run, but ask you to modify them in some cases. We have backups of everything, but please be careful so you don't delete something you shouldn't.



## Infrastructure

For the different analyses we are doing, we will use [Educloud](https://www.uio.no/english/services/it/research/platforms/edu-research/). To use it, you need an account which you can get here: [https://research.educloud.no/register](https://research.educloud.no/register). The project we are using in this course is ec146, so please ask for access to that one, and we will let you in. 

We will do the work in this course at `/projects/ec146/work` on [Fox](https://www.uio.no/english/services/it/research/platforms/edu-research/help/fox/) which is the HPC part of Educloud. After creating an account, you can log in using `ssh <educloud-username>@fox.educloud.no`. You will be prompted for a One-Time Code for a 2-factor authenticator app (Microsoft Authenticator) and your Fox/Educloud password.

On Fox we will submit jobs/analyses as job scripts. This is for a system called SLURM. Basically, this is instructions to the system for what kind of analysis we are running, or more concretely, how much memory and computing power we need. 

A generic job script might look like this (copied from [Job Scripts on Fox](https://www.uio.no/english/services/it/research/platforms/edu-research/help/fox/jobs/job-scripts.md)):
```
#!/bin/bash

# Job name:
#SBATCH --job-name=YourJobname
#
# Project:
#SBATCH --account=ecXXX
#
# Wall time limit:
#SBATCH --time=DD-HH:MM:SS
#
# Other parameters:
#SBATCH ...

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load SomeProgram/SomeVersion
module list

## Do some work:
YourCommands
```
## Other resources about assembly

We are not the first to create a workshop or tutorials about genome assembly. Here are a couple of good sources for more information.

Vertebrate Genomes Project as a workflow they have implemented in Galaxy, but with enough detail that you can run it yourself outside of Galaxy. They also link to a bit of background material. You can reach it here: 
[https://training.galaxyproject.org/training-material//topics/assembly/tutorials/vgp_genome_assembly/tutorial.html](https://training.galaxyproject.org/training-material//topics/assembly/tutorials/vgp_genome_assembly/tutorial.html)

UC Davis Bioinformatics Core has a genome assembly workshop they held in 2020, which also contains a lot of useful information if you want to deep-dive into this subject. Read it here: [https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/](https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/)



|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/01_GenomeScope2.md)|
|---|
