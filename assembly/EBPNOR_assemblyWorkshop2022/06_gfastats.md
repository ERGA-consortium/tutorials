# gfastats tutorial

Now it´s time to evaluate the assemblies you have created, starting with **gfastats**. This tool summarises some important assembly metrics such as number of scaffolds and contigs, scaffold and contig N50, and the number of gaps in the assembly. If you want to learn more about gfastats, click [here.](https://github.com/vgl-hub/gfastats) When you are ready, run the code and try to answer the questions below.

## Evaluating the gfastats results

Run this code in the terminal command line, and fill out the table below.


```
eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 
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


|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/05_YaHS.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/07_BUSCO.md)|
|---|---|