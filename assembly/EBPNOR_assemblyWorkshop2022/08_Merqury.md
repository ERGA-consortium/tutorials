# Merqury tutorial

Another way to validate your assemblies is by using **Merqury**. This k-mer based tool compares k-mers from the unassembled reads to the assemblies you have created, to find the degree of k-mer completeness, i.e. the percentage of the k-mers from the unassembled reads that are found within your assemblies. The advantage of this genome validation tool is that it can be used without any references, which is optimal when assessing de novo assemblies. If you want to learn more about Merqury, click [*here.*](https://github.com/marbl/merqury)

## Running Merqury

To run Merqury, create a new directory named `Merqury` and copy the code in the chunk below into a `run.sh` file:

```
#!/bin/bash
#SBATCH --job-name=merqury
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=5

eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate merqury

#$1 reads,  $2 first asm, $3 second asm, $4 output prefix 

ln -s $1 .

f=$(basename $1)

j=${f%.fastq.gz}
meryl k=21 threads=10 memory=8g count output $j.meryl $f

merqury.sh $j.meryl $2 $3 > $4_merqury.out 2> $4_merqury.err
```

Create a run.sh script with the following content (modify so it corresponds to what you have):

```
sbatch /projects/ec146/scripts/run_merqury.sh reads.fastq.gz first.fasta second.fasta prefix
```

## Interpreting a Merqury assembly spectrum plot

![Merqury plot 1](https://user-images.githubusercontent.com/110542053/206440295-74db4b51-d5f8-43d4-961f-c027c0f080af.png)

![Merqury plot 2](https://user-images.githubusercontent.com/110542053/206440369-ef889d7f-08ee-4f77-8a87-c0b4042fdb9d.png)

These are two of the Merqury plots generated for the EBP-Nor assembled Svalbard reindeer (*Rangifer tarandus*), specifically the spectra_asm plots. These are especially useful for evaluating haploid species assemblies, like *Metschnikowia zobellii*. This is because they can detect k-mers that are shared between both haplotypes (shown in green), and k-mers that are unique to each haplotype (shown in blue and red). The grey shaded area indicates the reads that are not found in either assemblies, and this gives us an indication about the degree of k-mer completeness, as mentioned in the introduction. 

Open a new terminal window, navigate to a folder where you want to place your Merqury png-files, and copy them using this line of code:

```
scp -r <username>@fox.educloud.no:/projects/ec146/work/<username>/merqury/"*.png" .
```

Look at the k-mer spectrum plots that you generated for *Metschnikowia zobellii* and answer these questions:

1. Does one of the haplotype have more unique reads than the other? Or are the haplotypes pretty similar?

2. Are a lot of the k-mers shared between the assemblies, or are most of the kmers within the assemblies unique to their own haplotype?


|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/07_BUSCO.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/09_FCS_GX.md)|
|---|---|
