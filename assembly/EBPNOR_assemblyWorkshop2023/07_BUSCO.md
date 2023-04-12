# BUSCO tutorial

How do you measure how “complete” your assembly is? Since we know our sawfly assemblies are about the expected size (which we estimated using GenomeScope2), and that our assembly statistics are good (which we found out by using gfastats), we now want to know that the genes we expect to find are there. **BUSCO** is a tool that can be used to compare your assemblies to lists of near universal orthologous genes (i.e. genes with common ancestry and function), to find out if you have managed to assemble them correctly. The assumption is that if you use the correct lineage dataset (based on what organism you are assembling), you´ll be able to find most of the genes within your assembly. If you want to read more about BUSCO, click [*here.*](https://busco.ezlab.org/busco_userguide.html)


## Running BUSCO

This is the script we use to run BUSCO:

```
#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --account=nn9984k
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks-per-node=5

eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 

conda activate busco

prefix=${1%.*}

mkdir -p busco_${prefix}
origdir=$PWD
cd busco_${prefix}

busco -c 5 -i ${origdir}/$1 -l /cluster/projects/nn9984k/opt/busco_downloads/lineages/hymenoptera_odb10 -o assembly -m genome  --offline > busco.out 2> busco.err
``` 

Create a new directory in your work area named `busco`. Make a new `run.sh` file with `nano run.sh`, and copy the code below into that file:

```
ln -s ../yahs/iyAthRosa_scaffolds_final.fa .
sbatch /cluster/projects/nn9984k/scripts/run_busco.sh iyAthRosa_scaffolds_final.fa
```

You can also run the script in the `hifiasm` and/or `yahs` folder directly. The script will make a subfolder with the results (look at the file called `busco.out` in the subfolder).

When you have done this, you can submit to the cluster by typing `sh run.sh`. If you scaffolded both haplotypes, repeat the process for both assemblies.


## Reviewing the BUSCO results

Before you start reviewing your results, look at the [list of lineages](https://busco-data.ezlab.org/v5/data/lineages/) found within the OrthoDB database. Which lineage would you pick for the yeast assemblies? Is it the same as the one we specified in the script above?

**Try to answer these questions by reviewing your BUSCO results:**

1. Which of the haplotypes has the highest percentage of complete BUSCO genes?

2. What do you think it means for the assembly quality if the number of fragmented, missing or duplicated BUSCO genes are high? Discuss with the people on your table. 

We have also downloaded more lineages than what is used in the script above. Take a look in `/cluster/projects/nn9984k/opt/busco_downloads/lineages/`. What do these different lineages mean? Can it be useful to use a less specific lineage in BUSCO?

(Hint: `busco  --list-datasets` shows all possible datasets/lineages and their hierarchy.) 

## Software versions used
```
eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 
conda activate busco
conda list
```
metaeuk version 6.a5d39d9 

busco version 5.4.5

|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/06_gfastats.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/08_Merqury.md)|
|---|---|
