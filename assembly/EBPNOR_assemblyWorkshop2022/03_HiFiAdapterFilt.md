# HiFiAdapterFilt tutorial

Now that you have learned a bit more about your dataset you are *almost* ready to start the assembly process. Before you can start, you need to remove any remaining adapter sequences that may still be attached to your HiFi reads. To do this, we in EBP-Nor use **HiFiAdapterFilt**. To learn more about how this software works, click [*here*](https://github.com/sheinasim/HiFiAdapterFilt), and to do it yourself, follow the tutorial below.

## Filtering adapter sequences with HiFiAdapterFilt

```
#!/bin/bash
#SBATCH --job-name=hifiadaptfilt
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=5

eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate hifiadapterfilt

export PATH=/fp/projects01/ec146/opt/HiFiAdapterFilt/:$PATH
export PATH=/fp/projects01/ec146/opt/HiFiAdapterFilt/DB:$PATH

pbadapterfilt.sh -t 5
```

We have set up this script for you. What you need to do is to create a run.sh in your working folder (`/projects/ec146/work/<username>/hifiadaptfilt`) with this content (with nano for instance):

```
ln -s /fp/projects01/ec146/data/genomic_data/pacbio/gsMetZobe_pacbio.fastq.gz . 
sbatch /projects/ec146/scripts/run_hifiadaptfilt.sh
```  
When you have done this, you can submit to the cluster by typing `sh run.sh`.

This should finish in a handful of minutes (when testing it ran for 1 minute). You can monitor the progress with `squeue -u <username>`.

HiFiAdapterFilt creates several files, for instance the filtered file and a statistics file: 

```
Started on Wed Dec  7 11:31:02 CET 2022
For the gsMetZobe_pacbio dataset:
Removing reads containing adapters a minimum of 44 bp in length and 97% match.

Number of ccs reads: 41700
Number of adapter contaminated ccs reads: 12 (0.028777% of total)
Number of ccs reads retained: 41688 (99.9712% of total)

Finished on Wed Dec  7 11:31:48 CET 2022
```

You should have similar content.



|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/02_Smudgeplot.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/04_hifiasm.md)|
|---|---|