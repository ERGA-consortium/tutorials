# gfastats tutorial

Now it´s time to evaluate the assemblies you have created, starting with **gfastats**. This tool summarises some important assembly metrics such as number of scaffolds and contigs, scaffold and contig N50, and the number of gaps in the assembly. If you want to learn more about gfastats, click [here.](https://github.com/vgl-hub/gfastats) When you are ready, run the code and try to answer the questions below.

## Evaluating the gfastats results

Run this code in the terminal command line, and fill out the table below.


```
eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 

conda activate gfastats

gfastats your_assembly.fa
```
The two first lines enable access to the conda system, and the gfastats environment in particular. 


Metric | Value
-------|-------
Number of scaffolds |
Total scaffold length |
Average scaffold length |
Scaffold N50 |
Largest scaffold |
Number of contigs |
Total contig length |
Average contig length |
Contig N50 |
Number of gaps |
Total gap length | 
Average gap length |
Gap N50 |


Why do you think each of these metrics are important for evaluating the quality of your assemblies? Discuss with the person sitting next to you.

When everything is done (hifiasm + YaHS), you can run gfastats on all the different assemblies (iyAthRosa.hic.hap1.p_ctg.fa, iyAthRosa.hic.hap2.p_ctg.fa and iyAthRosa_scaffolds_final) and compare then against each other. (Hint: These can also be found at `/cluster/projects/nn9984k/data/assemblies` if you are missing one or several.) 

The sequecinging data and approach is chosen to be able to fulfill the criteria/standards Earth Biogenome Project has for genome assemblies. You can review these [here](https://www.earthbiogenome.org/assembly-standards). Does this assembly or these assemblies fulfill the N50 contig length criterium?

## Software versions used
```
eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 
conda activate gfastats
conda list
```
gfastats version 1.3.6

|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/05_YaHS.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/07_BUSCO.md)|
|---|---|
