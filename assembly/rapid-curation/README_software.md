# Rapid curation - Software guide
  
Slack Channel for curation assembly help:
```
https://join.slack.com/t/assemblycuration/shared_invite/zt-yezxfd4w-P0xJdV1TJg47OaQKJaqlhA
```

This README will cover the software side of things and how it works.

## Set up

First directory variables must be set up, this will tell the singularity image
where the working directory and the data are located.
```
WORKDIR=/Where/you/will/be/working
DESTDIR=/Where/the/output/will/go
```

For Example:
```
WORKDIR=/user/grit/curation/organism/
DESTDIR=/user/grit/curation/organism/output
```

Next, we must group these variables together and bind them to a location inside
the singularity image.
```
export SINGULARITY_BIND="
/nfs:/nfs,\ <-- This will need changing depending on your working environment.
/lustre:/lustre,\ <-- This will need chaniging depending on where data is stored.
$WORKDIR:/data,\ <-- This does not need changing
$DESTDIR:/output,\ <-- This does not need changing
"
``` 

$SINGULARITY_BIND mounts the specified paths to paths inside the singularity image.
Files on paths not noted here will not be be visible and cause the pipeline to fail.
Only the text before the ':' should be changed to the base directory e.g. `/User, /nfs, /lustre`


For Example:
```
export SINGULARITY_BIND="
/user:/nfs,\
/user:/lustre,\
$WORKDIR:/data,\
$DESTDIR:/output,\
"
```

## Running the individual steps
### Command Flags
EXPLAIN

EXPLAIN bsub is LSF specific, if you use SLURM you need to change your job command

### Extra Notes
*_done files found in the output directories are empty files which aid in the control of
the underlying snakemake pipeline.

*.conf files in the same folders detail the commands used by the pipeline inside of the Singularity environment.

### HiC - Generates the HiGlass files and Pretext Map
<details open>
<summary> Open </summary>

This script generates the HiGlass and Pretext files required for curation of you genome.

The output of this program includes:
- ref.fa		<-- Original fasta file
- ref.fa.amb
- ref.fa.ann
- ref.fa.bwt		<-- Burrows Wheeler Transform file
- ref.fa.fai		<-- Fasta index file
- ref.fa.pac
- ref.fa.sa
- MORE

Example, used in GRIT:
```
echo '/software/singularity-v3.6.4/bin/singularity run /lustre/scratch123/tol/teams/grit/yy5/sigs/runHiC.sif -q 0 -s ilWatBina1_dr'|bsub -J test_sig -q basement -o test_sig.o -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'
```

##### Scripts and indepth usage
<details>
run-HiC uses a number of pieces of software:

- SamTools

- Python2.7

- Java8

- SALSA

- PretextMap

- bwa

- bamToBed

- juicer_tools 1.8.9_jcuda.0.8

- picard 2.18.11

- bammarkduplicates2

- arima_mapping_pipeline

##### Requires:
- fofn - file of input BAM or CRAM
- ref-fa - genome sequence in FASTA format

##### Commands:
These have been provided as the most basic command line commands

1 - Indexing
<details>
`bwa index {assembly}.fasta` & `samtools faidx {assembly}.fasta`
</details>

2 - Mapping
<details>
This step splits cram/bam into two files into the two sets of reads.
The top line may need modifying dependant on your own files. This needs to be performed for each cram file.

```
rgline=$(samtools view -H {bam1} | grep "@RG"| perl -spe 's/\t/\\t/g')
samtools view -hf 0x40 {bam1} | samtools fastq - | bwa mem -t15 -B8 -H'$rgline' {assembly}.fa - | perl arima_mapping_pipeline/filter_five_end.pl | samtools view -@4 -b - > bam1.mem.filt.1
samtools view -hf 0x80 {bam1} | samtools fastq - | bwa mem -t15 -B8 -H'$rgline' {assembly}.fa - | perl arima_mapping_pipeline/filter_five_end.pl | samtools view -@4 -b - > bam1.mem.filt.2
```
</details>

3 - Combine
<details>
These reads should now be re-combined. 

```
perl arima_mapping_pipeline/two_read_bam_combiner.pl mem.filt.1 mem.filt.2 /software/grit/bin/samtools 0 | samtools sort -@12 -T outfile_1.mem.filt.paired.sort.tmp -o outfile_1.mem.filt.paired.sort - 
```
</details>

4 - Merge
<details>
Multiple re-combined bams must be combined, an index is also produced.

```
samtools merge -@12 - {list of bams} | tee merge.bam | samtools index -c -@4 - {assembly}.csi
```
</details>

5 - Mark Duplicates
<details>
Mark Duplicates and create index.

```
bammarkdup2 I=merge.bam O=merge.mkdup.bam M=merge.mkdup.bam.metrics.txt tmpfile={outdir}/bammkdup2 markthreads=16
samtools index -c merge.mkdup.bam
```
</details>

6 - Bam2Bed
<details>

```
samtools view -@4 -u -F0x400 merge.mkdup.bam | bamToBed | sort -k4 --parallel=8 -S50G >merge.mkdup.bed
```
</details>

7 - HiGlass Start
<details>
Produce a .genome file and pre.bed to produce HiGlass cool and mcool files.

```
cut -f1,2 {assembly}.fa.fai | sed 's/-/_/g'|sort -k2,2 -nr > {assembly}.genome

paste -d '\t' - - < merge.mkdup.bed | sed 's/-/_/g' | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($1 > $7) {print substr($4,1,length($4)-2),$12,$7,$8,"16",$6,$1,$2,"8",$11,$5} else { print substr($4,1,length($4)-2),$6,$1,$2,"8",$12,$7,$8,"16",$5,$11} }' | tr '\-+' '01'  | sort --parallel=8 -S10G -k3,3d -k7,7d > pre.bed
```
</details>

8 - HiGlass End
<details>

```
cooler cload pairs -0 -c1 3 -p1 4 -c2 7 -p2 8 {assembly}.genome:1000 pre.bed {assembly}.cool
cooler zoomify -o {assembly}.mcool {assembly}.cool
```
</details>

9 - Pretext
<details>

```
samtools view -h merge.bam | pretextMap -o {assembly}.pretext --sortby length --mapq 0
```
</details>

</details>

</details>

### Coverage Track - Coverage plot used by HiGlass and Pretext
<details open>
<summary> Open </summary>

EXPLANATION

The output of this program includes:
- coverage.bed			<--
- geval.bed			<--
- {sample}.bw			<-- HiGlass and Pretext
- {input Pacbio fasta}.bam	<-- Input fasta converted to bam
- merged.bam			<-- Multiple bams are merged
- merged_sort.bam		<-- Merged bam is then sorted
- my.genome			<--
- pre.bam			<--
- ref.mmi			<-- Minimap Index File
- sort_tmp/			<--

Example, used in GRIT:
```
echo "/software/singularity-v3.6.4/bin/singularity run /nfs/team135/yy5/sif/runCoverage.sif -t ilWatBina1_dr" | bsub -J coverage -q basement -o coverage.o -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'
```
</details>

### Repeat Density Track - Needed for Pretext
<details open>
<summary> Open </summary>

EXPLANATION

The output of this program includes:
- bin.bed		<--
- density_nodot.bed	<--
- ref.fa.density.bed	<-- 
- ref.fa.genome		<--
- ref.fa.intersect.fa	<--
- ref.fa.repeat.bed	<--
- ref.fa.stage1		<--
- ref.fa.wm		<--
- sorted.genome		<--
- sorted_intersect.bed	<--
- {sample}_repeat_density.bw <-- HiGlass and Pretext

Example, used in GRIT:
```
echo "/software/singularity-v3.6.4/bin/singularity run /nfs/team135/yy5/sif/runRepeat.sif -t ilWatBina1_dr" | bsub -J repeat -q basement -o repeat.o -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'
```
</details>

### Gap Track - Needed for HiGlass and Pretext
<details open>
<summary> Open </summary>

The Gap Track utilises SeqTK to identity strings of N's wich indicate a gap in the sequence. This can then be used to inform curators on possible locations they can break scaffolds.
The output of this program is:

- {sample}_gap.bed	<-- HiGlass
- {sample}_gap.bedgraph <-- Pretext

seqtk can be found on GitHub [here](https://github.com/lh3/seqtk).

Example, used in GRIT:
```
echo "/software/singularity-v3.6.4/bin/singularity run /nfs/team135/yy5/sif/runGap.sif -t ilWatBina1_dr" | bsub -J gap -q basement -o gap.o -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'
```
</details>

### Telomere Track - Needed for HiGlass and Pretext
<details open>
<summary> Open </summary>

The Telomere Track uses the find_telomere program created by the VGP and modified to accept a custom telomere motif to search a given genome.
The output of this program include:

- {sample}_telomere.bed      <-- HiGlass
- {sample}_telomere.bedgraph <-- Pretext
- ref.telomere		     <-- 
- ref.windows		     <-- 

The original VGP's find_telomere files are found [here](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere).	

Example, used in GRIT:
```
echo "/software/singularity-v3.6.4/bin/singularity run /nfs/team135/yy5/sif/runTelo.sif -t ilWatBina1 -s ilWatBina1_dr" | bsub -J telo -q basement -o telo.o -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'
```
</details>

