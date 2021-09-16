
## Module to compute Differentially Accessibile Regions (DARs) and Similarly Accessibile Regions (SARs) genome-wide.

![Fig3_4Au](https://user-images.githubusercontent.com/61786710/133560437-1959f2d9-287e-44ca-af7f-dd17058cb9ee.png)


here one can find a set of computatitonal tools used to identify differentially accessible regions (DARs) and similarly accessible regions (SARs) genome-wide from the GC-methylation data of bulked scNOME_seq. The output is in bed-file format. 

The detailed description of the method in in original manuscript:
> 'Integrative analysis of transcriptomic and epigenomic data reveals distinct patterns for developmental and housekeeping genes'.

#### Main scripts and functions
``compute_dar_sar_main.m`` : 


#### A. default option

**usage** 

To run it, a user should submit minimal amount of input fils and the names/path of ouput files:

`[DAR,SARH]=compute_dar_sar_main(chrN,name_chr,file_win,file_countEct,file_countEnd,file_countMes,outDAR,outSAR);`

where
- chrN is an interger (here, ``chrN=13``)
- name_chr is a string (here, ``'chr13'``)
- `file_win,file_countEct,file_countEnd,file_countMes - names and paths of input files`
- ``outDAR,outSAR``- names and paths of output files

In this case the default values of the five thresholds below:
```
   thr_cov=25;
   thr_dif=0.51;
   thr_sim=0.36; 
   thr_accH=0.35;
   thr_accL=0.20;
   ```
   
will be generated within the function.

In the test folders the GC counts for chromosome 13 are stored, and the pathes/names of input and output files are generated for the
har-coded function for compactness: ``[file_win,file_countEct,file_countEnd,file_countMes,outDAR,outSAR]=test_file_names(name_chr);`` 

Therefore, to run it from a command line one should run the following lines:
```
   name_chr='chr13';chrN=13;[file_win,file_countEct,file_countEnd,file_countMes,outDAR,outSAR]=test_file_names(name_chr);   
   [DAR,SARH]=compute_dar_sar_main(chrN,name_chr,file_win,file_countEct,file_countEnd,file_countMes,outDAR,outSAR);
```

 #### B. User-defined thresholds
 
 If the user wants to vary parameters,e.g. to filter with smaller coverage, then they should  submit 13 inputs, where the last five
 are user defined thresholds
 
 e.g. **usage** is:
 ```
    thr_cov=15;
    thr_dif=0.53;
    thr_sim=0.35; 
    thr_accH=0.35;
    thr_accL=0.20; 
    name_chr='chr13';chrN=13;[file_win,file_countEct,file_countEnd,file_countMes,outDAR,outSAR]=test_file_names(name_chr);
           [DAR,SARH]=compute_dar_sar_main(chrN,name_chr,file_win,file_countEct,file_countEnd,file_countMes,outDAR,outSAR,thr_cov,thr_dif,thr_sim,thr_accH,thr_accL);
```

#### Description of the test example 

The main function used for reading the GC-data (arranged per chromosome), computing it in a fixed size window, filtering it by coverage, separately for each lineage. It further computes a dominance index for each window across the three lineages. It uses this dominance index to infer DARs and SAR as most dominant and even windows. The user should specify the input variables: GC-files/path, window file/path and output path. Currently, it is hard-coded in the beinning of the script.

#### Test data

``test_data`` folder contains examples of GC-methylation data for each lineage (pseudo-bulked) chr 13, precomputed coordinates of each 100bp window for chromosome 13 and expected outputs of DARs and SARs for chr13.

- the function test_file_names.m
generates names of the given test files arranged in the test_data folder: 
```
[file_win,file_countEct,file_countEnd,file_countMes,outDAR,outSAR]=test_file_names(name_chr); 
```

#### Dependencies
- ``prefilter_count_chr.m`` : function used to prefilter windows based on specified coverage threhold.

- ``dar_sarHL_levGC_dominance.m ``: function used to compute a dominance index for each of prefiltered windows across lineages, and select
DARs and SARs based on it.

 - ``save_chromatin_layers.m`` : function used to write sets of DARs and SARs in bed format into specified output file/path


*Requires MatlabR2014b or octave*
