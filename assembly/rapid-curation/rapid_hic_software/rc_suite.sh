#!/bin/bash

SINGULARITY_DIR=/nfs/team135/yy5/sif
SINGULARITY_BIN=/software/singularity-v3.6.4/bin
baseWORKDIR=/lustre/scratch123/tol/teams/grit/geval_pipeline/geval_runs/rc
baseDESTDIR=/lustre/scratch123/tol/teams/grit/geval_pipeline/geval_runs/rc


while getopts "p:s:t:" arg; do
  case $arg in
    h)
      echo "usage" 
      ;;
    p)
	  myprocess=$OPTARG
	  ;;
    s)
      teloseq=$OPTARG
      ;;
    t)
      sample=$OPTARG
      ;; 
  esac
done


# Check an assembly prefix has been specified
if [ -z "$myprocess" ]; then
  echo "ERROR: you must specify process that you want to run" >&2
  exit 1
fi

if [ -z "$sample" ]; then
  echo "ERROR: you must specify sample name to run" >&2
  exit 1
fi

WORKDIR=$baseWORKDIR/${sample}/data
DESTDIR=$baseDESTDIR/${sample}
logdir=$baseDESTDIR/${sample}/logs
mkdir -p ${logdir}

IFS=', ' read -ra process2run <<<${myprocess}; declare -p process2run;


for each in "${process2run[@]}"
do
  if [ "$each" = "hic" ]; then 
  	DESTDIR=$baseDESTDIR/${sample}/hic
  	mkdir -p ${DESTDIR}
    export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output," 
    echo "${SINGULARITY_BIN}/singularity run ${SINGULARITY_DIR}/runHiC.sif -q 0 -s ${sample}"|bsub -J runHic -q basement -o ${logdir}/hic.o -e ${logdir}/hic.e -n 10 -M80000 -R'select[mem>80000] rusage[mem=80000] span[hosts=1]'

  elif [ "$each" = "coverage" ]; then 
  	DESTDIR=$baseDESTDIR/${sample}/coverage
  	mkdir -p ${DESTDIR}
    export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output," 
    echo "${SINGULARITY_BIN}/singularity run ${SINGULARITY_DIR}/runCoverage.sif -t ${sample}"|bsub -J runCoverage -q long -o ${logdir}/coverage.o -e ${logdir}/coverage.e -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'

  elif [ "$each" = "gap" ]; then 
  	DESTDIR=$baseDESTDIR/${sample}/gap
  	mkdir -p ${DESTDIR}
    export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output," 
    echo "${SINGULARITY_BIN}/singularity run ${SINGULARITY_DIR}/runGap.sif -t ${sample}"|bsub -J runGap -q normal -o gap.o -n 2 -M5000 -R'select[mem>5000] rusage[mem=5000] span[hosts=1]'

  elif [ "$each" = "repeat" ]; then 
  	DESTDIR=$baseDESTDIR/${sample}/repeat
  	mkdir -p ${DESTDIR}
    export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output," 
    echo "${SINGULARITY_BIN}/singularity run ${SINGULARITY_DIR}/runRepeat.sif -t ${sample} -s 10000"|bsub -J runRepeat -q normal -o ${logdir}/repeat.o -e ${logdir}/repeat.e -n 5 -M10000 -R'select[mem>10000] rusage[mem=10000] span[hosts=1]'

  elif [ "$each" = "telomere" ]; then 
  	DESTDIR=$baseDESTDIR/${sample}/telomere
  	mkdir -p ${DESTDIR}
    export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output," 
    echo "${SINGULARITY_BIN}/singularity run ${SINGULARITY_DIR}/runTelo.sif -t ${sample} -s ${teloseq}"|bsub -J runTelomere -q normal -o ${logdir}/telomere.o -e ${logdir}/telomere.e -n 5 -M10000 -R'select[mem>10000] rusage[mem=10000] span[hosts=1]'
  else
    echo "process name is not valid"  
  fi;
done

