# FCS-GX tutorial

Contaminants can end up in your assemblies in various different ways. Maybe someone touched samples without gloves? Maybe there were symbionts living on the organism when it was sampled? Or maybe the sample was contaminated during the sequencing run? Luckily, there are several genomic decontamination tools available, and the one we use in EBP-Nor is the **NCBI Foreign Contamination Screen (FCS)** tool suite (click [here](https://github.com/ncbi/fcs) to read more). This program suite can identify and remove contaminant sequences, whether it is adaptor sequences, vector contamination or foreign organisms. In today´s workshop, we are going to focus on the latter, and below you can find the code to run your own decontamination script. 

## Decontaminating the yeast assembly

```
#!/bin/bash
#SBATCH --job-name=fcsgx
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=48G
#SBATCH --ntasks-per-node=10

export SHM_LOC=/fp/projects01/ec146/opt/fcs

echo "GX_NUM_CORES=10" > env.txt

python3 /fp/projects01/ec146/opt/fcs/run_fcsgx.py --fasta $1 \
--gx-db  "${SHM_LOC}/gxdb/all" --split-fasta --tax-id $2 \
--gx-db-disk "${SHM_LOC}/gxdb/all.gxi" \
--container-engine singularity --image /fp/projects01/ec146/opt/fcs/fcsgx.sif
```

As we have done earlier, we have set up this script for you. Create a run.sh in your working folder (`/projects/ec146/work/<username>/fcsgx`) with this content (with `nano` for instance):

```
sbatch /projects/ec146/scripts/run_gcsgx.sh assembly.fasta \
taxonomy_id
```
You have to modify the run.sh script based on your assembly file and you have to find the taxonomy ID for *Metschnikowia zobellii* and input that.

Unfortunately this program requires a lot of memory to run (["approximately 470 GiB"](https://github.com/ncbi/fcs/wiki/FCS-GX)). If it is given unsufficient memory, the running time can increase by a factor of 10000x. On Fox, there are not that [many nodes](https://www.uio.no/english/services/it/research/platforms/edu-research/help/fox/system-overview.md) with a lot of memory. The normal nodes have 501 GiB RAM, while the GPU accelerated nodes have up to 1006 GiB. Ideally, the job should have been allocated a bit more memory than what it strictly needs, but that is not easy here. Luckily, it should run in a handful of minutes if configured properly (1-30 minutes when tested). 

We should coordinate this, so only a couple people submit to the cluster. Let us know when you are at this point, and we can coordinate this.

After running the decontamination script, which foreign contaminants did you find?

If you were unable to run it, you can take a look at the results of a previous run at `/projects/ec146/data/fcsgx/gsMetZobe_scaffolds_final.27328.taxonomy.rpt`, `/projects/ec146/data/fcsgx/fcsgx.log` and `/projects/ec146/data/fcsgx/gsMetZobe_scaffolds_final.27328.fcs_gx_report.txt`. 

To remove contamination, you can do something like this:
```
eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate seqtk

grep ">" gsMetZobe_scaffolds_final.fa  |tr -d ">" |sort > all_sequences
grep EXCLUDE *fcs_gx_report.txt |cut -f 1 |sort > exclude_sequences

comm -23 all_sequences exclude_sequences > keep_sequences

seqtk subseq gsMetZobe_scaffolds_final.fa keep_sequences > gsMetZobe_clean.fa 
```
You can do it on the command line, or put it in a small script. Adjust for possible different file names.




|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/08_Merqury.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/10_Rapid_curation.md)|
|---|---|
