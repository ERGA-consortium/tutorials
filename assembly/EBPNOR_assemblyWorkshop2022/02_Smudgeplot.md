# Smudgeplot tutorial

The creators of GenomeScope2 created another way to visualize and estimate the ploidy and genome structure, using heterozygous k-mer pairs instead of k-mer frequency distribution. This software, called **Smudgeplot**, creates a “heatgraph” where a gradient from blue to yellow indicates the relative frequency of k-mer pairs. If the read coverage is good enough, you will have clear “smudges” of yellow which indicates the ploidy of your sequenced organism. The tutorial on how to create a smudgeplot can be found below, but if you want to read more about the software, you can click [*here.*](https://github.com/KamilSJaron/smudgeplot) 

## Creating a smudgeplot

```
#!/bin/bash
#SBATCH --job-name=smudgeplot
#SBATCH --account=ec146
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=5

eval "$(/fp/projects01/ec146/miniconda3/bin/conda shell.bash hook)" 

conda activate smudgescope

#reads as fastq
reads=$1
k=21
ploidy=2

mkdir -p tmp
echo $reads > FILES

[ -s reads.kmc_suf ] || kmc -k$k -t5 -m38 -ci1 -cs10000 @FILES reads tmp/

[ -s reads.histo ] || kmc_tools transform reads histogram reads.histo -cx10000

L=$(smudgeplot.py cutoff reads.histo L)
U=$(smudgeplot.py cutoff reads.histo U)
echo $L $U # these need to be sane values
# L should be like 20 - 200
# U should be like 500 - 3000

kmc_tools transform reads -ci"$L" -cx"$U" dump -s kmcdb_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o kmcdb_L"$L"_U"$U" < kmcdb_L"$L"_U"$U".dump

smudgeplot.py plot kmcdb_L"$L"_U"$U"_coverages.tsv
```

We have set up this script for you. What you need to do is to create a run.sh in your working folder (`/projects/ec146/work/<username>/smudgeplot`) with this content (with nano for instance): 
 
```
sbatch /projects/ec146/scripts/run_smudgeplot.sh /fp/projects01/ec146/data/genomic_data/pacbio/gsMetZobe_pacbio.fastq.gz
```

When you have done this, you can submit to the cluster by typing `sh run.sh`.
 
This should finish in a handful of minutes (when testing it ran for 1.5 minutes). You can monitor the progress with `squeue -u <username>`.

Smudgeplot produces several files in addition to the plot itself. You can for instance look at `smudgeplot_verbose_summary.txt` which contain the same information as the plot, but in text.

  
## Interpreting your smudgeplot

This is a smudgeplot generated for the EBP-Nor river lamprey (*Lampetra fluviatilis*):

![smudgeplot_riverlamprey](https://user-images.githubusercontent.com/110542053/206215771-1649b262-b685-4946-a869-397ff69ce533.png)

Here we see that based on the coverage of the heterozygous kmer pairs, the most represented haplotype structure for the river lamprey is diploid. 

Which is the most representet haplotype for *Metschnikowia zobellii*? To look at the plots you created, open a new terminal window, and navigate to a directory where you want to place your files. When you have found the place where you want to save them, use this code to copy them to your local computer:

```
scp -r <username>@fox.educloud.no:/projects/ec146/work/<username>/smudgeplot/"*.png" .
```


|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/01_GenomeScope2.md)|[Next](https://github.com/ebp-nor/genome-assembly-workshop-2022/blob/main/03_HiFiAdapterFilt.md)|
|---|---|