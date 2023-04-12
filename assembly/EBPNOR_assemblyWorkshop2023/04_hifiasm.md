# hifiasm tutorial

When choosing an assembler, you need to keep your data in mind. Since we want to create haplotype resolved assemblies, and we have both HiFi and Hi-C reads available, we are going to use **hifiasm**. Hifiasm is fast, easy to use, and creates high quality assemblies with longer contigs. To read more about how this software works, click [*here.*](https://github.com/chhylp123/hifiasm)

## Assembling with hifiasm

```
#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --account=nn9984k
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --ntasks-per-node=10

eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 

conda activate hifiasm

hifiasm -o $1 -t10  \
--h1 $2 \
--h2 $3 \
$4 \
1> hifiasm_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> hifiasm_"`date +\%y\%m\%d_\%H\%M\%S`".err

awk '/^S/{print ">"$2"\n"$3}' $1.hic.hap1.p_ctg.gfa | fold > $1.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' $1.hic.hap2.p_ctg.gfa | fold > $1.hic.hap2.p_ctg.fa
```

We have set up this script for you. What you need to do is to create a run.sh in your working folder (`/cluster/projects/nn9984k/work/<username>/hifiasm`) with this content (with nano for instance): 
 
```
sbatch /cluster/projects/nn9984k/scripts/run_hifiasm.sh iyAthRosa \
/cluster/projects/nn9984k/data/genomic_data/hic/ERR6054981_1_60x.fastq.gz \
/cluster/projects/nn9984k/data/genomic_data/hic/ERR6054981_2_60x.fastq.gz \
/cluster/projects/nn9984k/data/genomic_data/pacbio/iyAthRosa_pacbio.fastq.gz 
```
This script contain the unfiltered HiFi reads. Please replace the reads with the filtered reads you created with HiFiAdapterFilt.

When you have done this, you can submit to the cluster by typing `sh run.sh`.
 
This will run for a while. When testing it ran for 3.2 hours. When next needing the assembly, we will use one which is already done. You can monitor the progress with `squeue -u <username>`.

## Software versions used
```
eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 
conda activate hifiasm
conda list
```
hifiasm version 0.18.5


|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/03_HiFiAdapterFilt.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/05_YaHS.md)|
|---|---|
