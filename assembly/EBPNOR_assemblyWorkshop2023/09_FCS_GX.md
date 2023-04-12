# FCS-GX tutorial

Contaminants can end up in your assemblies in various different ways. Maybe someone touched samples without gloves? Maybe there were symbionts living on the organism when it was sampled? Or maybe the sample was contaminated during the sequencing run? Luckily, there are several genomic decontamination tools available, and the one we use in EBP-Nor is the **NCBI Foreign Contamination Screen (FCS)** tool suite (click [here](https://github.com/ncbi/fcs) to read more). This program suite can identify and remove contaminant sequences, whether it is adaptor sequences, vector contamination or foreign organisms. In today's workshop, we are going to focus on the latter, and below you can find the code to run your own decontamination script. 

## Decontaminating the sawfly assembly

```
#!/bin/bash
#SBATCH --job-name=fcsgx
#SBATCH --account=nn9984k
#SBATCH --partition=bigmem
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --ntasks-per-node=10

eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 
conda activate base
#The base system seems to have python2, this sets up python3 which is needed for fcs

#this copies the database to a shared memory. This works, but not sure how safe it is.
mkdir /dev/shm/fcs   #make directory. 
rsync -av /cluster/projects/nn9984k/opt/fcs/gxdb /dev/shm/fcs

export SHM_LOC=/dev/shm/fcs/gxdb

#export SHM_LOC=/cluster/projects/nn9984k/opt/fcs/gxdb

echo "GX_NUM_CORES=10" > env.txt

python3 /cluster/projects/nn9984k/opt/fcs/run_fcsgx.py --fasta $1 \
--gx-db  "${SHM_LOC}/all" --split-fasta --tax-id $2 \
--gx-db-disk "${SHM_LOC}/all.gxi" \
--container-engine singularity --image /cluster/projects/nn9984k/opt/fcs/fcsgx.sif
```

As we have done earlier, we have set up this script for you. Create a run.sh in your working folder (`/projects/ec146/work/<username>/fcsgx`) with this content (with `nano` for instance):

```
sbatch /cluster/projects/nn9984k/scripts/run_fcsgx.sh assembly.fasta \
taxonomy_id
```
You have to modify the run.sh script based on your assembly file and you have to find the NCBI taxonomy ID for *Athalia rosae* and input that.

Unfortunately this program requires a lot of memory to run (["approximately 470 GiB"](https://github.com/ncbi/fcs/wiki/FCS-GX)). If it is given unsufficient memory, the running time can increase by a factor of 10000x. On Saga, there are not that [many nodes](https://documentation.sigma2.no/hpc_machines/saga.html) with a lot of memory. There are 8 so-called bigmem nodes, which should be abel to handle multiple jobs each of the script above. However, these are quite heavily used, so it is not certain that we will be able to run our jobs here. If configured properly, it is quite quick (1-30 minutes when testing). 

We should coordinate this, so only a couple people submit to the cluster. Let us know when you are at this point.

After running the decontamination script, which foreign contaminants did you find?

If you were unable to run it, you can take a look at the results of a previous run at ` /cluster/projects/nn9984k/data/fcsgx/iyAthRosa_scaffolds_final.37344.taxonomy.rpt`, `/cluster/projects/nn9984k/data/fcsgx/fcs.log ` and `/cluster/projects/nn9984k/data/fcsgx/iyAthRosa_scaffolds_final.37344.fcs_gx_report.txt`. 

To remove contamination, you can do something like this:
```
eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 

conda activate seqtk

grep ">" iyAthRosa_scaffolds_final.fa  |tr -d ">" |sort > all_sequences
grep EXCLUDE *fcs_gx_report.txt |cut -f 1 |sort > exclude_sequences

comm -23 all_sequences exclude_sequences > keep_sequences

seqtk subseq iyAthRosa_scaffolds_final.fa keep_sequences > iyAthRosa_clean.fa 
```
You can do it on the command line, or put it in a small script. Adjust for possible different file names.




|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/08_Merqury.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/10_Rapid_curation.md)|
|---|---|
