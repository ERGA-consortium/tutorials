# YaHS tutorial

Congratulations, you have created your sawfly assembly! But now you have to combine your contigs into scaffolds, and for that we use **YaHS**. YaHS stands for “yet another Hi-C scaffolding tool”, and as the name implies, there are a lot of Hi-C scaffolders out there. However, we in EBP-Nor choose to use YaHS because it is fast, creates more contiguous scaffolds, with better genome statistics compared to other widely used scaffolders. To learn more about this, click [*here*](https://github.com/c-zhou/yahs), otherwise scroll down to start scaffolding your haplotype resolved assemblies.

## Scaffolding with YaHS

```
#!/bin/bash
#SBATCH --job-name=yahs
#SBATCH --account=nn9984k
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks-per-node=5

eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 

conda activate yahs

REF=$1

[ -s $REF.bwt ] || bwa index $REF

SAMPLE=$2

mkdir -p outs

[ -s hic_markdup.sort_n.bam ] || bwa mem -t 8 -R '@RG\tSM:$SAMPLE\tID:$SAMPLE' -5SPM $REF \
$3 $4 \
|samtools view -buS - |samtools sort -@1 -n -T tmp_n -O bam - \
|samtools fixmate -mr - -|samtools sort -@1 -T hic_tmp -O bam - |samtools markdup -rs - -  2> hic_markdup.stats |samtools sort -n -@1 -n -T temp_n -O bam\
> hic_markdup.sort_n.bam

#markdup -S is not supported by samtools 1.6

[ -s $REF.fai ] ||samtools faidx $REF

if [ -s $SAMPLE.bin ]; then
        yahs $REF $SAMPLE.bin -o $SAMPLE \
        1> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".err
else
        yahs $REF hic_markdup.sort_n.bam -o $SAMPLE \
        1> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".err
fi

```

As we did with hifiasm, we have set up this script for you. Create a run.sh in your working folder (`/cluster/projects/nn9984k/work/<username>/yahs`) with this content (with `nano` for instance):

```
ln -s ../hifiasm/iyAthRosa.hic.hap1.p_ctg.fa .

#or link from /cluster/projects/nn9984k/data/assemblies/iyAthRosa.hic.hap1.p_ctg.fa if you are not done yet with the assembly

sbatch /cluster/projects/nn9984k/scripts/run_yahs.sh iyAthRosa.hic.hap1.p_ctg.fa \
iyAthRosa \
/cluster/projects/nn9984k/data/genomic_data/hic/ERR6054981_1_60x.fastq.gz \
/cluster/projects/nn9984k/data/genomic_data/hic/ERR6054981_2_60x.fastq.gz
```

When you have done this, you can submit to the cluster by typing `sh run.sh`.

When testing, this ran for 1.6 hours. You can monitor the progress with `squeue -u <username>`.

Here we have simplified matters a bit and only proposes to scaffold one of the assemblies. Usually, we filter the Hi-C data based on unique k-mers in the two assemblies. This would lead to the reads used to scaffold hap1 would not contain reads with k-mers that are unique to hap2, and vice versa. 

## Software versions used
```
eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 
conda activate yahs
conda list
```
bwa version 0.7.17 

samtools version 1.6

yahs version 1.2a.2

|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/04_hifiasm.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/06_gfastats.md)|
|---|---|
