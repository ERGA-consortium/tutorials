# HiFiAdapterFilt tutorial

Now that you have learned a bit more about your dataset you are *almost* ready to start the assembly process. Before you can start, you need to remove any remaining adapter sequences that may still be attached to your HiFi reads. To do this, we in EBP-Nor use **HiFiAdapterFilt**. To learn more about how this software works, click [*here*](https://github.com/sheinasim/HiFiAdapterFilt), and to do it yourself, follow the tutorial below.

## Filtering adapter sequences with HiFiAdapterFilt

```
#!/bin/bash
#SBATCH --job-name=hifiadaptfilt
#SBATCH --account=nn9984k
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=5

eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 

conda activate hifiadapterfilt

export PATH=/cluster/projects/nn9984k/opt/HiFiAdapterFilt/:$PATH
export PATH=/cluster/projects/nn9984k/HiFiAdapterFilt/DB:$PATH

pbadapterfilt.sh -t 5
```

We have set up this script for you. What you need to do is to create a run.sh in your working folder (`cluster/projects/nn9984k/work/<username>/hifiadaptfilt`) with this content (with nano for instance):

```
ln -s /cluster/projects/nn9984k/data/genomic_data/pacbio/iyAthRosa_pacbio.fastq.gz  . 
sbatch /cluster/projects/nn9984k/scripts/run_hifiadaptfilt.sh
```  
When you have done this, you can submit to the cluster by typing `sh run.sh`.

This should finish in a handful of minutes (when testing it ran for about 20 minutes). You can monitor the progress with `squeue -u <username>`.

HiFiAdapterFilt creates several files, for instance the filtered file and a statistics file: 

```
Started on Sat Feb  4 08:58:25 CET 2023
For the iyAthRosa_pacbio dataset:
Removing reads containing adapters a minimum of 44 bp in length and 97% match.

Number of ccs reads: 410000
Number of adapter contaminated ccs reads: 205 (0.05% of total)
Number of ccs reads retained: 409795 (99.95% of total)

Finished on Sat Feb  4 09:18:25 CET 2023
```

You should have similar content.

## Software versions used
```
eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 
conda activate hifiadapterfilt
conda list
```
HiFiAdapterFilt version Second release

blast version 2.13.0

bamtools version 2.5.2

|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/02_Smudgeplot.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/04_hifiasm.md)|
|---|---|
