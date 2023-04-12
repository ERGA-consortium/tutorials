# YaHS tutorial

Congratulations, you have created your yeast assembly! But now you have to combine your contigs into scaffolds, and for that we use **YaHS**. YaHS stands for “yet another Hi-C scaffolding tool”, and as the name implies, there are a lot of Hi-C scaffolders out there. However, we in EBP-Nor choose to use YaHS because it is fast, creates more contiguous scaffolds, with better genome statistics compared to other widely used scaffolders. To learn more about this, click [*here*](https://github.com/c-zhou/yahs), otherwise scroll down to start scaffolding your haplotype resolved assemblies.

## Scaffolding with YaHS

```
#!/bin/bash
#SBATCH --job-name=yahs
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=10

eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate yahs

REF=$1

[ -s $REF.bwt ] || bwa index $REF

SAMPLE=$2

mkdir -p outs

[ -s hic_markdup.sort_n.bam ] || bwa mem -t 10 -R '@RG\tSM:$SAMPLE\tID:$SAMPLE' -5SPM $REF \
$3 $4 \
|samtools view -buS - |samtools sort -@1 -n -T tmp_n -O bam - \
|samtools fixmate -mr - -|samtools sort -@1 -T hic_tmp -O bam - |samtools markdup -rs - -  2> hic_markdup.stats |samtools sort -n -@1 -T temp_n -O bam \
> hic_markdup.sort_n.bam

[ -s $REF.fai ] ||samtools faidx $REF

if [ -s $SAMPLE.bin ]; then
        yahs $REF $SAMPLE.bin -o $SAMPLE \
        1> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".err
else
        yahs $REF hic_markdup.sort_n.bam -o $SAMPLE \
        1> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> yahs_"`date +\%y\%m\%d_\%H\%M\%S`".err
fi

```

As we did with hifiasm, we have set up this script for you. Create a run.sh in your working folder (`/projects/ec146/work/<username>/yahs`) with this content (with `nano` for instance):

```
ln -s ../hifiasm/gsMetZobe.hic.hap1.p_ctg.fa .

sbatch /projects/ec146/scripts/run_yahs.sh gsMetZobe.hic.hap1.p_ctg.fa \
gsMetZobe \
/fp/projects01/ec146/data/genomic_data/hic/ERR9503460_1_60x.fastq.gz \
/fp/projects01/ec146/data/genomic_data/hic/ERR9503460_2_60x.fastq.gz 
```

When you have done this, you can submit to the cluster by typing `sh run.sh`.

This should finish in a handful of minutes (when testing it ran for 6 minutes). You can monitor the progress with `squeue -u <username>`.

Here we have simplified matters a bit and only proposes to scaffold one of the assemblies. Usually, we filter the Hi-C based on unique k-mers in the two assemblies. So that the reads used to scaffold hap1 would not contain reads with k-mers that are unique to hap2, and vice versa. 


|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/04_hifiasm.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/06_gfastats.md)|
|---|---|
