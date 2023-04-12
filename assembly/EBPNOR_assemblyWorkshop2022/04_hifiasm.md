# hifiasm tutorial

When choosing an assembler, you need to keep your data in mind. Since we want to create haplotype resolved assemblies, and we have both HiFi and Hi-C reads available, we are going to use **hifiasm**. Hifiasm is fast, easy to use, and creates high quality assemblies with longer contigs. To read more about how this software works, click [*here.*](https://github.com/chhylp123/hifiasm)

## Assembling with hifiasm

```
#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks-per-node=5

eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate hifiasm

hifiasm -o $1 -t5  \
--h1 $2 \
--h2 $3 \
$4 \
1> hifiasm_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> hifiasm_"`date +\%y\%m\%d_\%H\%M\%S`".err

awk '/^S/{print ">"$2"\n"$3}' $1.hic.hap1.p_ctg.gfa | fold > $1.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' $1.hic.hap2.p_ctg.gfa | fold > $1.hic.hap2.p_ctg.fa
```


We have set up this script for you. What you need to do is to create a run.sh in your working folder (`/projects/ec146/work/<username>/hifiasm`) with this content (with nano for instance): 
 
```
sbatch /projects/ec146/scripts/run_hifiasm.sh gsMetZobe \
/fp/projects01/ec146/data/genomic_data/hic/ERR9503460_1_60x.fastq.gz \
/fp/projects01/ec146/data/genomic_data/hic/ERR9503460_2_60x.fastq.gz \
/fp/projects01/ec146/data/genomic_data/pacbio/gsMetZobe_pacbio.fastq.gz
```
This script contain the unfiltered HiFi reads. Please replace the reads with the filtered reads you created with HiFiAdapterFilt.

When you have done this, you can submit to the cluster by typing `sh run.sh`.
 
This should finish in a handful of minutes (when testing it ran for 25 minutes). You can monitor the progress with `squeue -u <username>`.


|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/03_HiFiAdapterFilt.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/05_YaHS.md)|
|---|---|