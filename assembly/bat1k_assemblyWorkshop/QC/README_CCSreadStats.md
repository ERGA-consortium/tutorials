# 1. Get initial statistics from the CCS reads

[seqkit](https://bioinf.shenwei.me/seqkit/usage/#stats) can be used to get some initial read statistics from the CCS reads 

```bash 
# 1. store singularity command (with all arguments) into a variable to facilitate the command line call
SING_CMD="singularity exec --no-home --cleanenv -B /projects /projects/dazzler/pippel/prog/assembly-workshop/assembly-workshop_v0.6.3.sif"
# 2. lets assume all PacBio ccs fastq files are in following /projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data
CCS_DIR=/projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data
ls ${CCS_DIR}/*.ccs.fastq.gz > pb.fofn 
# 3. check your pb.fofn
cat pb.fofn 
# this should give you somthing like
#/projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data/m64015_191126_012040.ccs.fastq.gz
#/projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data/m64015_191127_073410.ccs.fastq.gz
#/projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data/m64046_190805_101354.ccs.fastq.gz
#/projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data/m64046_190807_094011.ccs.fastq.gz
#/projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data/m64046_190813_152807.ccs.fastq.gz
#/projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data/m64046_190819_115653.ccs.fastq.gz
#/projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data/m64046_191007_114257.ccs.fastq.gz
#/projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data/m64046_191022_094945.ccs.fastq.gz 
# 4.1 run seqkit stats, on a local machine 
for x in $(cat pb.fofn)
do 
    out=$(basename ${x%.*.*}).readStats.txt
    ${SING_CMD} seqkit stats -a ${x} > ${out}
    ## 4.2 if you have enough cores and IO-bandwidth you can run the seqkit call in parallel in the background  
    ## seqkit usually uses 2-cores 
    # ${SING_CMD} seqkit stats -a ${x} > ${out} &
done 
# OR 4.3 submit it to a compute cluster, e.g. on SLURM cluster trhe submission could look like:
##  sbatch -c 2 -n 1 -p batch -a 1-$(wc -l < pb.fofn) -J seqkit --wrap="${SING_CMD} seqkit stats -a \$(sed -n \${SLURM_ARRAY_TASK_ID}p pb.fofn) > \$(basename \$(sed -n \${SLURM_ARRAY_TASK_ID}p pb.fofn) | sed -e 's:.fastq.gz:.readStats.txt:')" 
```

check the results: `cat *.readStats.txt`

|file|format|type|num_seqs|sum_len|min_len|avg_len|max_len|Q1|Q2|Q3|sum_gap|N50|Q20(%)|Q30(%)|
|:---|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|m64015_191126_012040.ccs.fastq.gz|FASTQ|DNA|782,884|14,945,850,409|44|19,090.8|38,463|18,350|19,295|20,143|0|19,392|98.36|96.18|
|m64015_191127_073410.ccs.fastq.gz|FASTQ|DNA|730,066|13,946,895,644|45|19,103.6|38,718|18,362|19,304|20,152|0|19,402|98.4|96.27|
|m64046_190805_101354.ccs.fastq.gz|FASTQ|DNA|1,437,262|18,907,954,820|45|13,155.5|42,537|11,642|12,730|14,357|0|13,078|98.76|97.16|
|m64046_190807_094011.ccs.fastq.gz|FASTQ|DNA|1,367,406|18,080,287,584|50|13,222.3|39,874|11,683|12,793|14,436|0|13,151|98.74|97.16|
|m64046_190813_152807.ccs.fastq.gz|FASTQ|DNA|1,068,654|16,992,227,735|48|15,900.6|49,720|13,213|15,024|17,659|0|15,818|98.44|96.45|
|m64046_190819_115653.ccs.fastq.gz|FASTQ|DNA|933,213|14,980,358,470|47|16,052.5|47,286|13,327|15,185|17,842|0|15,984|98.37|96.28|
|m64046_191007_114257.ccs.fastq.gz|FASTQ|DNA|77,436|1,429,015,292|43|18,454.1|48,026|15,059|17,785|21,356|0|19,154|98.05|95.53|
|m64046_191022_094945.ccs.fastq.gz|FASTQ|DNA|276,770|5,189,355,592|45|18,749.7|49,356|15,168|17,956|21,644|0|19,334|97.58|94.48|

back to [main QC page](./README.md)