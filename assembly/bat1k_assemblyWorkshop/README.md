# Bat1K assembly workshop

The Bat1K assembly workshop provides an introduction into _de novo_ genome assemblies based on [PacBio HiFi](https://www.pacb.com/smrt-science/smrt-sequencing/hifi-reads-for-highly-accurate-long-read-sequencing/) and [Illumina HiC](https://arimagenomics.com/) reads. 
The assembly pipeline will be presented step-by-step beginning with a QC of the input reads, contig assembly, haplotype purging, HiC-scaffolding and ending with the manual curation step. Our sequencing data is from one individual of the bat species [_Tadarida brasiliensis_](https://en.wikipedia.org/wiki/Mexican_free-tailed_bat). Other PacBio HiFi data sets that could be used for training can be found here: [www.pacb.com](https://www.pacb.com/blog/hifi-data-release/).

The slides for the workshop are available as [google doc](https://docs.google.com/presentation/d/1R-ZdZCxhrK8-WGVK8SiRMcozJN2u9yr8XrHYozcNDiA/edit?usp=sharing). 

## Species overview

* Project: Bat1K
* Species: Tadarida brasiliensis (Mexican free-tailed bat)
* Estimated genome size: ~2.9Gb &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ([http://genomesize.com](http://genomesize.com))
* Expected chromosomes: 2n = 48  &nbsp;&nbsp;&nbsp;&nbsp; ([chromosome_size](https://git.mpi-cbg.de/dibrov/chromosome_size))
* Classification: [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9438)
* NCBI taxonomy ID: 9438

![Tadarida brasiliensis](images/tadarida_brasiliensis.jpeg)


&nbsp;**0. Software and Data Prerequisites** 

An overview of the used data sets is given here: [data](data)

An installation tutorial for the required software of this workshop can be found in [programs](programs). 


**&nbsp;1. Preprocessing and Quality Control (QC) of input reads**   
&nbsp;&nbsp;&nbsp;[Overview](QC/README.md)     
**&nbsp;&nbsp;&nbsp;1.1 create PacBio CCS reads (optional)**    
&nbsp;&nbsp;&nbsp;[HowTo](QC/README_CCS.md)  
**&nbsp;&nbsp;&nbsp;1.2 Initial QC of CCS reads with seqkit**      
&nbsp;&nbsp;&nbsp;[HowTo](QC/README_CCSreadStats.md)    
**&nbsp;&nbsp;&nbsp;1.3 Create kmer table with FastK**      
&nbsp;&nbsp;&nbsp;[HowTo](QC/README_FastK.md)    
**&nbsp;&nbsp;&nbsp;1.4 Create genome size estimate with GeneScopeFK**      
&nbsp;&nbsp;&nbsp;[HowTo](QC/README_GeneScopeFK.md)    
**&nbsp;&nbsp;&nbsp;1.5 MerquryFK**      
&nbsp;&nbsp;&nbsp;[HowTo](QC/README_MerquryFK.md)

**&nbsp;2. Genome assembly with hifiasm**     
&nbsp;presented in the [google doc](https://docs.google.com/presentation/d/1R-ZdZCxhrK8-WGVK8SiRMcozJN2u9yr8XrHYozcNDiA/edit?usp=sharing/#slide=id.g10458f6996b_0_0)

**&nbsp;3. Purge_Dups**     
&nbsp;presented in the [google doc](https://docs.google.com/presentation/d/1R-ZdZCxhrK8-WGVK8SiRMcozJN2u9yr8XrHYozcNDiA/edit?usp=sharing/#slide=id.g10458f6996b_0_146)

**&nbsp;4. HiC scaffolding**        
&nbsp;presented in the [google doc](https://docs.google.com/presentation/d/1R-ZdZCxhrK8-WGVK8SiRMcozJN2u9yr8XrHYozcNDiA/edit?usp=sharing/#slide=id.g10458f6996b_0_242)

**&nbsp;5. Manual curation**   
&nbsp;The tool is included in the container (>=v0.6.3) but can also be downloaded from the MPI gitlab: [ManualCurationHiC](https://git.mpi-cbg.de/assembly/programs/manualcurationhic) 

&nbsp;presented in the [google doc](https://docs.google.com/presentation/d/1R-ZdZCxhrK8-WGVK8SiRMcozJN2u9yr8XrHYozcNDiA/edit?usp=sharing/#slide=id.g10458f6996b_0_218)

**&nbsp;6. Gap Closing**  
&nbsp; A pipeline, still adapted to MPI-CBG compute infrastructure, can be found on our MPI gitlab: 
[gap_closing](https://git.mpi-cbg.de/assembly/programs/gap_closing). We are in the progress to intergrate this into the container as well.

&nbsp;presented in the [google doc](https://docs.google.com/presentation/d/1R-ZdZCxhrK8-WGVK8SiRMcozJN2u9yr8XrHYozcNDiA/edit?usp=sharing/#slide=id.g1049e928914_5_13)


**&nbsp;7. Assembly Polishing**      
&nbsp; A pipeline, still adapted to MPI-CBG compute infrastructure, can be found on our MPI gitlab: 
[deepvariant_polishing](https://git.mpi-cbg.de/assembly/programs/polishing). We are in the progress to intergrate this into the container as well.

&nbsp;presented in the [google doc](https://docs.google.com/presentation/d/1R-ZdZCxhrK8-WGVK8SiRMcozJN2u9yr8XrHYozcNDiA/edit?usp=sharing/#slide=id.g10458f6996b_0_227)

**&nbsp;8. QV analysis**        
&nbsp;&nbsp;&nbsp;[HowTo](QV/README_MerqueryFK.md)  