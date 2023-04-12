In case your sequencing provider supports the raw subreads from the machine you need to create the circular consensus sequences from those subreads first. 
The PacBio tool [ccs](https://ccs.how/), which is contained in the assembly_workshop conda environemt and in the singularity container, can be used this. 

```bash 
    ccs [options] <IN.subreads.bam|xml> <OUT.ccs.bam|fastq.gz|xml>
```

As this step can be quite time consuming we recommend, to use the `--chunk` option and distribute the jobs on a compute cluster.  

The following code snippet shows how to run a ccs pipeline on a SLURM compute cluster with the given singularity container. We partition a SMRT cell into 140 chunks and use 24 cores. This ensures that each ccs run usually finishes within 20 minutes. 

```bash 
SING_CMD="singularity exec --no-home --cleanenv -B /projects /projects/dazzler/pippel/prog/assembly-workshop/assembly-workshop_v0.6.3.sif"
SMRT=m64015_191127_073410.subreads.bam 
outDir=./
outFilePrefix=$(basename ${SMRT%.subreads.bam}).ccs
# chunk data into N parts: 
CHUNKS=140
# threads per ccs job 
THREAD=24
# log level  TRACE, DEBUG, INFO, WARN, FATAL
LOGLEVEL=INFO
# limit concurrent number of jobs that run on a partition at the same time 
# if you don't want any limitation set it to CHUNKS
#CONCURRENT_JOBS=${CHUNKS}
CONCURRENT_JOBS=140
# MEM per job 
mem=16G
# job time limit
tim=04:00:00
# partition 
partition="batch"

# submit ccs jobs to the cluster
jid1=$(sbatch --array=1-${CHUNKS}%${CONCURRENT_JOBS} -c ${THREAD} -n 1 -p ${partition} --mem=${mem} -e ${outDir}/"${outFilePrefix}".%A.%a.err -o "${outDir}"/"${outFilePrefix}".%A.%a.out --time=${tim} -J ccs --wrap="${SING_CMD} ccs ${SMRT} ${outDir}/${outFilePrefix}.\${SLURM_ARRAY_TASK_ID}.bam --chunk \${SLURM_ARRAY_TASK_ID}/${CHUNKS} -j ${THREAD} --log-level ${LOGLEVEL} --log-file ${outDir}/${outFilePrefix}.\${SLURM_ARRAY_TASK_ID}.log --report-file ${outDir}/${outFilePrefix}.\${SLURM_ARRAY_TASK_ID}.report.txt && ccs --version | head -n 1 | tr -d \")(\" |  awk '{print \$NF}' > ${outDir}/${outFilePrefix}.\${SLURM_ARRAY_TASK_ID}.version")
echo "submit ccs jobs: ${jid1##* }"
# submit merge job to the cluster
jid1=$(sbatch -c ${THREAD} -n 1 --dependency=afterok:"${jid1##* }" -p ${partition} -J merge --mem=${mem} --time=${tim} --wrap="${SING_CMD} samtools merge -c -p -f -@${THREAD} ${outDir}/${outFilePrefix}.bam ${outDir}/${outFilePrefix}.*.bam")
echo "submit bam merge job: ${jid1##* }"
# create index file
jid1=$(sbatch -c ${THREAD} -n 1 --dependency=afterok:"${jid1##* }" -p ${partition} -J pbindex --mem=${mem} --time=${tim} --wrap="${SING_CMD} pbindex ${outDir}/${outFilePrefix}.bam")
echo "submit bam index job: ${jid1##* }"
# create fastq file
jid1=$(sbatch -c ${THREAD} -n 1 --dependency=afterok:"${jid1##* }" -p ${partition} -J bam2fastq --mem=${mem} --time=${tim} --wrap="${SING_CMD} bam2fasta -o ${outDir}/${outFilePrefix} ${outDir}/${outFilePrefix}.bam")
echo "submit bam2fastq job: ${jid1##* }"
```
back to [main QC page](./README.md)