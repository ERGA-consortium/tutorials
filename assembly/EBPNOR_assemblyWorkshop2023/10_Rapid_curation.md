# Rapid curation tutorial

Although hifiasm is a great assembler, and YaHS can create chromosome length scaffolds, assembly errors do occur. Whether there are contigs that are misassembled, or scaffolds that are harder for the software to place, sometimes we need to manually curate the assemblies in order to reach our EBP-Nor assembly standards (you can read more about that here). To do this, we use the **Rapid curation** suite, developed by the GRIT-team at the Wellcome Sanger Institute, and **PretextView**, which you´ll learn more about in the last tutorial. If you want to read more about the code you´ll be using today, click [here](https://gitlab.com/wtsi-grit/rapid-curation/-/blob/main/README_software.md), and if you want to read more about why curation is so important for good quality reference genomes, click [here.](https://academic.oup.com/gigascience/article/10/1/giaa153/6072294) 

## Running the Rapid curation suite

### Run the suite

As with the other programs in the assembly pipeline, we have set up a script for you (see the code chunk below):

```
#!/bin/bash
#SBATCH --job-name=curation
#SBATCH --account=nn9984k
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=5

eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 

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
singularity run /cluster/projects/nn9984k/opt/rapid-curation/rapid_hic_software/runHiC.sif -q 0 -s $1
#rm $HOME/hic_done

#coverage
minimap2 -ax map-hifi \
         -t 5 data/ref.fa \
	$3 \
| samtools sort -@16 -O BAM -o coverage.bam

samtools view -b -F 256 coverage.bam > coverage_pri.bam

samtools index coverage_pri.bam

bamCoverage -b coverage_pri.bam -o coverage.bw

#gaps
singularity run /cluster/projects/nn9984k/opt/rapid-curation/rapid_hic_software/runGap.sif -t $1

#repeats
singularity run /cluster/projects/nn9984k/opt/rapid-curation/rapid_hic_software/runRepeat.sif -t $1

#telomers
singularity run /cluster/projects/nn9984k/opt/rapid-curation/rapid_hic_software/runTelo.sif -t $1 -s TTAGG

cp out/out.pretext $1.pretext

#put it together
bigWigToBedGraph coverage.bw  /dev/stdout |PretextGraph -i $1.pretext -n "PB coverage"

cat out/*_gap.bedgraph  | PretextGraph -i $1.pretext -n "gaps"

cat out/*_telomere.bedgraph |awk -v OFS="\t" '{$4 *= 1000; print}' | PretextGraph -i $1.pretext -n "telomers"

bigWigToBedGraph  out/*_repeat_density.bw /dev/stdout | PretextGraph -i $1.pretext -n "repeat density"
```

This script creates both the Hi-C contact map that we'll use in PretextView, and the overlays we'll use to inform our edits during curation. 

### Starting the script

To run the script above, create a `run.sh` file in a new curation-directory, and run the code using `sh run.sh`. 

```
mkdir -p data
mkdir -p out

#or use your own
cat  /cluster/projects/nn9984k/data/fcsgx/iyAthRosa_clean.fa > data/ref.fa 

echo "/hic/ERR6054981_60x.bam" > data/cram.fofn

sbatch /cluster/projects/nn9984k/scripts/run_rapidcuration.sh iyAthRosa /cluster/projects/nn9984k/data/genomic_data/hic/  /cluster/projects/nn9984k/data/genomic_data/pacbio/iyAthRosa_pacbio.fasta.gz 

```

After this is finished, you should be left with an iyAthRosa.pretext file, and this can be used for manual curation. 

If you don´t want to wait for your scripts to finish, and you want to start curating right away, we have provided both the files you need to do so. To download the PRETEXT file to your local computer, open a new terminal window, navigate to where you want to place the file, and use the code below:

```
scp -r <username>@saga.sigma2.no:/cluster/projects/nn9984k/data/pretext/iyAthRosa.pretext .

```

You also need a TPF file to curate the sawfly assembly, and this is created from your decontaminated fasta. To create the TPF file, run the code below in your curation directory:

```
eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)"
conda activate curation

#or use your own
ln -s /cluster/projects/nn9984k/data/fcsgx/iyAthRosa_clean.fa

perl /cluster/projects/nn9984k/opt/rapid-curation/rapid_split.pl -fa iyAthRosa_clean.fa
```


### For information: converting fastq files to BAM
The rapid curation suite requires Hi-C reads to be in a BAM format. To create that, we did this:
```
#!/bin/bash
#SBATCH --job-name=convert_bam
#SBATCH --account=nn9984k
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=48G
#SBATCH --ntasks-per-node=10


module load picard/2.24.0-Java-11

java -Xms6g -Xmx6g -jar $EBROOTPICARD/picard.jar FastqToSam \
F1=ERR6054981_1_60x.fastq.gz \
F2=ERR6054981_1_60x.fastq.gz \
O=ERR6054981_60x.bam \
SM=hic_sawfly \
TMP_DIR=.
```


|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/09_FCS_GX.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/11_PretextView.md)|
|---|---|
