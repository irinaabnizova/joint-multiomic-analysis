## Module to compute Differentially methylated Regions (DMRs) and Similarly Methylated Regions (SMRs) genome-wide.

![Fig4new](https://user-images.githubusercontent.com/61786710/133770572-8e2bcab1-8897-4371-b7b3-5f20617decf0.png)

here one can find a set of computatitonal tools used to identify differentially methylated regions (DMRs) and similarly methylated regions (SMRs) genome-wide from the CpG-methylation data (BS_seq) of bulked scNOME_seq. The output is in bed-file format. 

The detailed description of the method in in original manuscript:
> 'Integrative analysis of transcriptomic and epigenomic data reveals distinct patterns for developmental and housekeeping genes'.

#### Main scripts and functions
``compute_dmr_smr_main_2.m`` : 


#### A. default option

**usage** 

To run it, use the code:
`run_main_chr13.m`

user should submit minimal amount of input fils and the names/path of ouput files:

`[DMR,SMRH]=compute_dmr_smr_main_2(chrN,name_chr,file_win,file_countEct,file_countEnd,file_countMes,outDMR,outSMR);`

where
- chrN is an interger (here, ``chrN=13``)
- name_chr is a string (here, ``'chr13'``)
- `file_win,file_countEct,file_countEnd,file_countMes - names and paths of input files`
- ``outDMR,outSMR``- names and paths of output files

In this case the default values of the five thresholds below:
```
   thr_cov=5;
   thr_dif=0.5;
   thr_sim=0.35; 
   thr_meH=0.7;
   thr_meL=0.35;
   ```
   
will be generated within the function.

In the test folders theCpG counts for chromosome 13 are stored, and the pathes/names of input and output files are generated within the example code.


 #### B. User-defined thresholds
 
 If the user wants to vary parameters,e.g. to filter with smaller coverage, then they should  submit/insert in the body of the test scripts  13 inputs, where the last five  are user defined thresholds
 
 e.g. **usage** is:
 ```
   thr_cov=5;
   thr_dif=0.51;
   thr_sim=0.36; 
   thr_meH=0.6;
   thr_meL=0.3;
   
    name_chr='chr13';chrN=13;[file_win,file_countEct,file_countEnd,file_countMes,outDMR,outSMR]=test_file_names(name_chr);
           [DMR,SMRH]=compute_dar_sar_main_2(chrN,name_chr,file_win,file_countEct,file_countEnd,file_countMes,outDMR,outSMR,thr_cov,thr_dif,thr_sim,thr_meH,thr_meL);
           
```

#### Description of the test example 

The main function used for reading the CpG (arranged per chromosome), computing it in a fixed size window, filtering it by coverage, separately for each lineage. It further computes a dominance index for each window across the three lineages. It uses this dominance index to infer DARs and SAR as most dominant and even windows. The user should specify the input variables: CpG-files/path, window file/path and output path. Currently, it is hard-coded in the beinning of the script.

#### Test data

``test_data`` folder contains examples of CpG-methylation data for each lineage (pseudo-bulked) chr 13, precomputed coordinates of each 100bp window for chromosome 13 and expected outputs of DMRs and SMRs for chr13.

- the function test_file_names.m
generates names of the given test files arranged in the test_data folder: 
```
[file_win,file_countEct,file_countEnd,file_countMes,outDMR,outSMR]=test_file_names(name_chr); 
```

#### Dependencies
- ``prefilter_count_chr.m`` : function used to prefilter windows based on specified coverage threhold.

- ``DMR_SMR_lev_cov_levCG_PuHL.m, just_DMR_thrSept.m, get_pu_ch_meth_3 ``: functions used to compute a dominance index for each of prefiltered windows across lineages, and select DMRs and SMRs based on it.

 - ``save_DMR_SMR_EctEndMes_H07_L05.m`` : function used to write sets of DMRs and SMRs in bed format into specified output file/path


*Requires MatlabR2014b or octave*
