# Rapid curation tutorial

Although hifiasm is a great assembler, and YaHS can create chromosome length scaffolds, assembly errors do occur. Whether there are contigs that are misassembled, or scaffolds that are harder for the software to place, sometimes we need to manually curate the assemblies in order to reach our EBP-Nor assembly standards (you can read more about that here). To do this, we use the **Rapid curation** suite, developed by the GRIT-team at the Wellcome Sanger Institute, and **PretextView**, which you´ll learn more about in the last tutorial. If you want to read more about the code you´ll be using today, click [here](https://gitlab.com/wtsi-grit/rapid-curation/-/blob/main/README_software.md), and if you want to read more about why curation is so important for good quality reference genomes, click [here.](https://academic.oup.com/gigascience/article/10/1/giaa153/6072294) 

## Running the Rapid curation suite

### Run the suite

As with the other programs in the assembly pipeline, we have set up a script for you (see the code chunk below):

```
#!/bin/bash
#SBATCH --job-name=curation
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=5

eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate curation

WORKDIR=$PWD/data
DESTDIR=$PWD/out
HICDIR=$2

export SINGULARITY_BIND="
$WORKDIR:/data,\
$HICDIR:/hic,\
$DESTDIR:/output,\
$TMP_DIR:/tmp
"

#hic
#singularity run /fp/projects01/ec146/opt/rapid-curation/rapid_hic_software/runHiC.sif -q 0 -s $1
#rm $HOME/hic_done
#this did not work for some reason, and we were unable to figure it out.

bwa index data/ref.fa

bwa mem -t 8 -5SPM data/ref.fa \
$4 $5 \
|samtools view -buS - |samtools sort -@1 -n -T tmp_n -O bam - \
|samtools fixmate -mr - -|samtools sort -@1 -T hic_tmp -O bam - |samtools markdup -rsS - -  2> hic_markdup.stats |samtools sort -@1 -T temp_n -O bam\
> hic_markdup.sort.bam

#coverage
minimap2 -ax map-hifi \
         -t 5 data/ref.fa \
	$3 \
| samtools sort -@16 -O BAM -o coverage.bam

samtools view -b -F 256 coverage.bam > coverage_pri.bam

samtools index coverage_pri.bam

bamCoverage -b coverage_pri.bam -o coverage.bw

#gaps
singularity run /fp/projects01/ec146/opt/rapid-curation/rapid_hic_software/runGap.sif -t $1

#repeats
singularity run /fp/projects01/ec146/opt/rapid-curation/rapid_hic_software/runRepeat.sif -t $1  -s 10000

#for some reason, all reads had mapq == 0, so we'll cheat:
samtools view -h hic_markdup.sort.bam | PretextMap -o $1.pretext --sortby length --sortorder descend --mapq 0

#telomers
#singularity run /fp/projects01/ec146/opt/rapid-curation/rapid_hic_software/runTelo.sif -t $1 -s $3
#skipping telomers since they are not regular in budding yeast, at least not to our knowledge

#put it together
bigWigToBedGraph coverage.bw  /dev/stdout |PretextGraph -i $1.pretext -n "PB coverage"

cat out/*_gap.bedgraph  | PretextGraph -i $1.pretext -n "gaps"

#cat out/*_telomere.bedgraph |awk -v OFS="\t" '{$4 *= 1000; print}' | PretextGraph -i $1.pretext -n "telomers"

bigWigToBedGraph  out/*_repeat_density.bw /dev/stdout | PretextGraph -i $1.pretext -n "repeat density"
```

This script creates both the Hi-C contact map that we´ll use in PretextView, and the overlays we´ll use to inform our edits during curation. 

### Starting the script

To run the script above, create a `run.sh` file in a new curation-directory, and run the code using `sh run.sh`. 

```
mkdir -p data
mkdir -p out

cat  /fp/projects01/ec146/data/fcsgx/gsMetZobe_clean.fa > data/ref.fa 

echo "/hic/hic_yeast.bam" > data/cram.fofn

sbatch /projects/ec146/scripts/run_rapidcuration.sh gsMetZobe /fp/projects01/ec146/data/genomic_data/hic/  /fp/projects01/ec146/data/genomic_data/pacbio/gsMetZobe_pacbio.fastq.gz /fp/projects01/ec146/data/genomic_data/hic/ERR9503460_1_60x.fastq.gz /fp/projects01/ec146/data/genomic_data/hic/ERR9503460_2_60x.fastq.gz 
```

After this is finished, you should be left with an out.pretext file, and this can be used for manual curation. However, since the yeast species you were working on had some ambiguous Hi-C contact signals, we´re gonna let you try your hand at curating the EBP-Nor brook lamprey instead! To download this file to your local computer, open a new terminal window, navigate to where you want to place your file, and use the code below:

```
scp -r <username>@fox.educloud.no:/projects/ec146/data/pretext_data/lampetra_planeri_hap1.pretext .

scp -r <username>@fox.educloud.no:/projects/ec146/data/pretext_data/lampetra_planeri_hap1.tpf .
```

### For information: converting fastq files to BAM
The rapid curation suite requires Hi-C reads to be in a BAM format. To create that, we did this:
```
#!/bin/bash
#SBATCH --job-name=convert_bam
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=48G
#SBATCH --ntasks-per-node=10


module load picard/2.24.0-Java-11

java -Xms6g -Xmx6g -jar $EBROOTPICARD/picard.jar FastqToSam \
F1=ERR9503460_1_60x.fastq.gz \
F2=ERR9503460_1_60x.fastq.gz \
O=hic_yeast.bam \
SM=hic_yeast \
TMP_DIR=.
```


|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/09_FCS_GX.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/11_PretextView.md)|
|---|---|
