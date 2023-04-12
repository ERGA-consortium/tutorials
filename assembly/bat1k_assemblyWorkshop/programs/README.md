
Most tools that are used in the assembly workshop can be installed with the following 3 options. 

There are two external packages DeepVariant and Higlass which need to be installed separately. 
We recommend to use the available Docker/Singularity conatainer (see E) External packages)

# A) Use a Singularity container (recommended)

Install Singularity from [sylabs.io](https://sylabs.io/guides/3.0/user-guide/installation.html). 
The singularity container `assembly-workshop_v0.6.3.sif` can be downloaded from the MPI Owncloud service: [assembly_workshop](https://cloud.mpi-cbg.de/index.php/s/ErYzv6Eae1lZfOa) and can be set up in the following way: 

```bash
# 1. create a directory for the container (full path)
SING_CONT_DIR=/projects/dazzler/pippel/prog/assembly-workshop
mkdir -p ${SING_CONT_DIR}
cd ${SING_CONT_DIR}

# 2. download container file
curl 'https://cloud.mpi-cbg.de/index.php/s/ErYzv6Eae1lZfOa/download?path=%2F&files=assembly-workshop_v0.6.3.sif'  -o assembly-workshop_v0.5.sif

# 3. test container 
singularity exec --no-home --cleanenv ${SING_CONT_DIR}/assembly-workshop_v0.6.3.sif bash /assembly_workshop.sh
singularity exec --no-home --cleanenv ${SING_CONT_DIR}/assembly-workshop_v0.6.3.sif hifiasm --version

# 4. to access files on your host system, you need to bind the host path(s) into the container (for more info: https://sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html) 
singularity exec --no-home --cleanenv -B /projects ${SING_CONT_DIR}/assembly-workshop_v0.6.3.sif ls ${SING_CONT_DIR}/assembly-workshop_v0.5.sif
```
back to [main assembly_workshop page](https://git.mpi-cbg.de/assembly/assembly_workshop)

# B) Use Anaconda ( ... and from source)
Install Anaconda from [https://docs.anaconda.com](https://docs.anaconda.com/anaconda/install/index.html#), create a new conda environment, and install the the software packages from the file `conda-linux-64.lock`.
```bash
# 1. create new environment from package file 
conda create --name assembly_workshop --file conda-linux-64.lock
# 2. load conda environment 
conda activate assembly_workshop
# 3. most of the assembly tools can be used now from the command line, e.g.  
hifiasm -help
```

But the following 3 software packages must be installed from source. (We are working on a anaconda version).
* [FastK](https://github.com/thegenemyers/FASTK)
* [MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK)
* [GENESCOPE.FK](https://github.com/thegenemyers/GENESCOPE.FK)

back to [main assembly_workshop page](https://git.mpi-cbg.de/assembly/assembly_workshop)

# C) Build from source 

* [hifiasm](https://github.com/chhylp123/hifiasm)
* [purge_dups](https://github.com/dfguan/purge_dups)
* [FastK](https://github.com/thegenemyers/FASTK)
* [MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK)
* [GENESCOPE.FK](https://github.com/thegenemyers/GENESCOPE.FK)
* [seqkit](https://github.com/shenwei356/seqkit)
* [pbmm2](https://github.com/PacificBiosciences/pbmm2)
* [gcpp](https://github.com/PacificBiosciences/gcpp)
* [ccs](https://github.com/PacificBiosciences/ccs)
* [bam2fastx](https://github.com/PacificBiosciences/bam2fastx)
* [samtools](https://github.com/samtools/samtools)
* [bcftools](https://github.com/samtools/bcftools)
* [freebayes](https://github.com/freebayes/freebayes)
* [bwa](https://github.com/lh3/bwa)
* [pairtools](https://github.com/open2c/pairtools)
* [busco](https://gitlab.com/ezlab/busco)

back to [main assembly_workshop page](https://git.mpi-cbg.de/assembly/assembly_workshop)
___________________________

## [D) build container from scratch]

Just for completeness (and documentation) here are all steps to recreate the containers (docker and singularity): 

1. The initial assembly tools were installed in a clean conda environment. 
```bash 
conda create --name assembly-workshop
conda activate assembly-workshop
conda config --add channels conda-forge
conda config --add channels bioconda  
conda install purge_dups pbmm2 pbccs pbgcpp samtools bcftools seqkit freebayes hifiasm bam2fastx bwa pairtools r-base r-argparse r-minpack.lm busco
```

2. Export (and adapt) conda environment  
```bash 
conda env export --from-history > environment_asm.yml
```

3. Use [conda-lock](https://github.com/conda-incubator/conda-lock) to render the requirements into locked pinnings for different architectures. If you need conda packages for other architectures you need to modify this step appropriatelty. More details can found on [Uwe's blog post](https://uwekorn.com/2021/03/01/deploying-conda-environments-in-docker-how-to-do-it-right.html).
```bash
# the following command creates the file conda-linux-64.lock
conda lock -p linux-64 -f environment_asm.yml
```

4. run the build.sh script to create docker and singularity containers.
```bash 
bash build.sh
```

5. **Troubleshooting**

In case you are building the singularity container on a compute cluster (as we did) you might get the following error: 
```bash 
INFO:    Starting build...
FATAL:   While performing build: conveyor failed to get: Error writing blob: write /tmp/bundle-temp-817373221/oci-put-blob352501850: no space left on device
``` 
If this is the case you can prepend the envionment variable SINGULARITY_TMPDIR=/path/onHost/WithEnough/Storage/tmp to the singularity build call in the build.sh script. E.g.:
``` bash  
SINGULARITY_TMPDIR=/path/onHost/WithEnough/Storage/tmp SINGULARITY_DISABLE_CACHE=true singularity build "$SINGULARITY_IMAGE" "docker-daemon://$DOCKER_TAG"
```

back to [main assembly_workshop page](https://git.mpi-cbg.de/assembly/assembly_workshop)

# E) External packages

## [DeepVariant](https://github.com/google/deepvariant)

DeepVariant is used in error correction of the assembly with PacBio HiFi data. We strongly recommend to do this.

```bash 
# easiest installation is via pulling the docker image 
BIN_VERSION=1.2.0
SINGULARITY_DISABLE_CACHE=true singularity pull docker://gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}"
```

## [higlass](https://github.com/higlass/higlass)

Higlass is used to visualize the HiC density plots of the assembly, especially in the manual curation.
We think easiest installation is via the available [HiGlass Docker container](https://github.com/higlass/higlass-docker). 
```bash  
# get the latest image 
docker pull higlass/higlass-docker # Ensure that you have the latest.
# run the docker image NOTE: you need to modify the paths that are bound to the container 
docker run --detach \
           --publish 8888:80 \
           --volume ~/hg-data:/data \
           --volume ~/hg-tmp:/tmp \
           --name higlass-container \
           higlass/higlass-docker
```

back to [main assembly_workshop page](https://git.mpi-cbg.de/assembly/assembly_workshop)
