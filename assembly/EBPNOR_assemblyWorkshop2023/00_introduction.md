# Welcome to the Genome assembly, curation and validation workshop!

Today you'll learn how to assemble a whole genome, EBP-Nor style! 

After attending the workshop you should:
- know about the most-used approaches for genome assembly
- be able to assess information inherit in sequencing reads
- be able to validate genome assemblies
- know about manual curation of assemblies


## Our dataset

The coleseed or turnip sawfly, *Athalia rosae*, is a sawfly found in Europe, Asia, North America and Africa. It is often considered a pest, feeding on plants of the brassica family, such as rapeseed, turnip, mustard and cabbage. 

[Darwin Tree of Life](https://www.darwintreeoflife.org) has sequenced the coleseed sawfly, but not published a genome not for it yet. Fortunately for us, they allow everyone to play with the data anyhow, so we will do that. DToL has some interesting webpages where they list several quality measures for the sequencing (some of which we will do in this workshop) [here](https://tolqc.cog.sanger.ac.uk/darwin/insects/Athalia_rosae/). It can be worth a look. We downloaded the data from [ENA](https://www.ebi.ac.uk/ena/browser/view/GCA_917208135) and subsampled it to get it to the coverages we expect/plan for. 

The genome itself is around 170 Mbp, and the PacBio data was about 18 Gbp, more than 100x coverage. It was subsampled with [seqtk](https://github.com/lh3/seqtk) like this:
```
seqtk sample ERR6548410.fastq.gz  410000 |gzip > ERR6548410_30x.fastq.gz
```
to get about 30x coverage with PacBio data. 


## Why do we use a combination of HiFi and Hi-C reads? 

HiFi sequencing creates highly accurate circularized consensus reads. How are these reads generated? By ligating hairpin adapters, the DNA fragment that is being sequenced becomes a circle. This means that the machine can do multiple passes over the same DNA sequence, to weed out any misread nucleotides. This is how HiFi reads can relatively long, while remaining over 99.9% accurate. 

Hi-C sequencing is done to capture how the chromatin is folded within the cell nucleus. By ligating the folded DNA strands, we can capture which loci are found in close proximity, and thus which parts of the DNA are found within the same chromosomes.

When combining these two, we can create haplotype resolved assemblies, meaning we can separate reads by maternal and paternal origin, without having access to parental data. In diploid, or polyploid organisms, this adds another level of information, and creates more accurate assemblies than a primary and alternate assembly would. 

Testing, by us, but also earlier by Darwin Tree of Life and Vertebrate Genomes Project, among others, has shown that the combination of HiFi and Hi-C, in appropriate coverages, usually generates assemblies that fulfill the Earth Biogenome Project's criteria for [assembly standards](https://www.earthbiogenome.org/assembly-standards). There are other ways to get to these standards, by using combinations of Oxford Nanopore Technologies sequencing data and Illumina, but these are often less straight-forward and involves more steps to a final assembly that the strategy we outline here.

## Package management

Administrating the different programs that are needed in project can be a hassle. We like conda, especially [miniconda](https://docs.conda.io/en/latest/miniconda.html), and have set up different environments we will use for the different analyses. [Bioconda](https://bioconda.github.io) contain a lot of different packages that are relevant for us, and genomics and bioinformatics in general.

To load conda, do this:
```
eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 
```

There are some of the different programs that are not available through conda. For most of these we use [Singularity containers](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html). 

We mostly set up scripts and arranged data so it is ready to run, but ask you to modify them in some cases. We have backups of everything, but please be careful so you don't delete something you shouldn't.

## Infrastructure

For the different analyses we are doing, we will use [Saga](https://documentation.sigma2.no/hpc_machines/saga.html). Everything should be set up properly by now. The project we have at Saga is called nn9984k, and the working folder is `/cluster/projects/nn9984k`. You should set up and do stuff in `/cluster/projects/nn9984k/work/$USERNAME`, but we'll come back to that in the next subject.

On Saga we will submit jobs/analyses as job scripts. This is for a system called SLURM. Basically, this is instructions to the system for what kind of analysis we are running, or more concretely, how much memory and computing power we need. 

A generic job script might look like this (copied from [Saga](https://documentation.sigma2.no/hpc_machines/saga.html)):
```
#SBATCH --account=MyProject
#SBATCH --job-name=MyJob
#SBATCH --time=1-0:0:0
#SBATCH --mem-per-cpu=3G
#SBATCH --ntasks=16

Some commands.
```
## Other resources about assembly

We are not the first to create a workshop or tutorials about genome assembly. Here are a couple of good sources for more information.

Vertebrate Genomes Project as a workflow they have implemented in Galaxy, but with enough detail that you can run it yourself outside of Galaxy. They also link to a bit of background material. You can reach it here: 
[https://training.galaxyproject.org/training-material//topics/assembly/tutorials/vgp_genome_assembly/tutorial.html](https://training.galaxyproject.org/training-material//topics/assembly/tutorials/vgp_genome_assembly/tutorial.html)

UC Davis Bioinformatics Core has a genome assembly workshop they held in 2020, which also contains a lot of useful information if you want to deep-dive into this subject. Read it here: [https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/](https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/)



|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/01_GenomeScope2.md)|
|---|
