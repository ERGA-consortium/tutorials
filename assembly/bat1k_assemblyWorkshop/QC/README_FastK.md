# Run FastK on the CCS data 

```bash 
# 1. store singularity command (with all arguments) into a variable to facilitate the command line call
SING_CMD="singularity exec --no-home --cleanenv -B /projects /projects/dazzler/pippel/prog/assembly-workshop/assembly-workshop_v0.6.3.sif"
``` 
To get an overview about the FastK arguments just run `${SING_CMD} FastK` on the command line. A more in depth introduction of all FastK tools can be found on Gene Myer's github repo [FastK](https://github.com/thegenemyers/FASTK).

```
Usage: FastK [-k<int(40)>] -t[<int(4)>]] [-p[:<table>[.ktab]]] [-c] [-bc<int(0)>]
               [-v] [-N<path_name>] [-P<dir(/tmp)>] [-M<int(12)>] [-T<int(4)>]
                 <source>[.cram|.[bs]am|.db|.dam|.f[ast][aq][.gz] ...

      -v: Verbose mode, output statistics as proceed.
      -T: Use -T threads.
      -N: Use given path for output directory and root name prefix.
      -P: Place block level sorts in directory -P.
      -M: Use -M GB of memory in downstream sorting steps of KMcount.

      -k: k-mer size.
      -t: Produce table of sorted k-mers & counts >= level specified
      -p: Produce sequence count profiles (w.r.t. table if given)
     -bc: Ignore prefix of each read of given length (e.g. bar code)
      -c: Homopolymer compress every sequence

```

In our test scenario we will use the arguments `-v -t1 -T24 -k40 -NmTadBra1 -P./tmp`. By default intermediate files (e.g. for sorting) are written into `/tmp`. For some data sets these can get quite big. Therefore we suggest to explicitely set `-P/path`. Note that between the argument switch and the following parameter e.g. `-T24` **no space is allowed**. 

```bash
# 2. create a local tmp folder
mkdir -p tmp
# 3. lets assume all PacBio ccs fastq files are in following /projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data
# 4. run FastK
$SING_CMD FastK -v -t1 -T24 -k40 -NmTadBra1 -P./tmp /projects/dazzlerAssembly/asm_mTadBra/Bat1K_workshop/data/*.ccs.fastq.gz
```

The above FastK command finished in less then 17 minutes on our test data set and created two files `mTadBra1.ktab`, `mTadBra1.hist`, and multiple hidden files `.mTadBra1.ktab.1 .mTadBra1.ktab.2 ... .mTadBra1.ktab.38`. 


back to [main QC page](./README.md)