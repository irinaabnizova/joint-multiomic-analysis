
## Specificity module: codes to compute how specific are two sets of TF motifs:

![Fig5new](https://user-images.githubusercontent.com/61786710/134467015-f0e34627-799f-410f-be7c-726bc5154601.png)

It consists of the following steps:

    - compute feature-vector/matrix for each TF motif list
    - computing pairwise similarities between lists of enriched TFBS motifs between and within two given sets: e.g. DEGs core=promoters and SEGs core prooters.
    - compute median similarity values within  and between sets of TF lists
    - performs permutation test for significance of specificity: suffles ranks and motifs within lists
    
  

The detailed description of the method in in the original manuscript:
> 'Integrative analysis of transcriptomic and epigenomic data reveals distinct patterns for developmental and housekeeping genes'.

**This set of tools is used for producing:**
 - similarity matrix 2x2 for median values of within and between sets similarity values
 - together with specificity value (as relative proportion of difference between similarities)  
 - probability  of obtaing the same specificity by chance given similar local nucleotide content


#### Main scripts and functions
``MAIN_compute_promoter_specificity.m:`` 


**usage** 

## A. default option

to run it, user should submit minimal amount of input files and the names/path of ouput files:

``[matrix_sim_Original,vectors_sim_Shuffled,specificity,prob]=MAIN_compute_promoter_specificity(Feature_List,list_SEG,list_Ect,list_End,list_Mes);``


where

Input files are ``Feature_List, list_SEG,list_Ect,list_End,list_Mes``:

- ``Feature_List`` - list of all TF motifs in the given mice databases, alphabetical order

- ``list_SEG,list_Ect,list_End,list_Mes`` are TFBS lists within SARs and three DARs

#### data_specificity data

In the data_specificity folders the motif lists for promoters and enhancers (DEGs and SEGs) are stored, 
and the pathes/names of input and output files are generated for the hard-coded function for compactness: 

``[FLCC,list_SEG,list_Ect,list_End,list_Mes,outOriginal,outShuffled]=test_enhancer_names()`` 

Therefore, to run it from a command line one should run the following lines:

```  
[FLCC,list_SEG,list_Ect,list_End,list_Mes,outOriginal,outShuffled]=test_enhancer_names();  
[matrix_sim_Original,vectors_sim_Shuffled,specificity,prob]=MAIN_compute_specificity(FLCC,list_SEG,list_Ect,list_End,list_Mes);

```

### B. User defined parameters

If the user wants to vary parameters,e.g. then they should  submit 7 inputs, where the last two
 are user defined thresholds for the length of lists (in proportions, should be less than 1) and user defined quartile weight values:

e.g. usage is:

```
    fS=1;quart_weights=[0.8,0.4,0.2,0.1];
    [FLCC,list_SEG,list_Ect,list_End,list_Mes,outOriginal,outShuffled]=test_enhancer_names(); 
    [matrix_sim_Original,vectors_sim_Shuffled,specificity,prob]=MAIN_compute_specificity(FLCC,list_SEG,list_Ect,list_End,list_Mes,fS,quart_weights);
```


#### Dependencies

The functions reading data, computing stats and visualisation:
  
``permute_original_loop.m``          : function used to compute original and shuffled matrices of similarity
             
``M2_fromM4_noDiag.m``               : function used to transform 4x4 into 2x2 metrix     

``put_Zero_rank_uniform_bw_a_b.m ``  : function used to shuffle ranks and positions of TF motifs          
   
``quartileControl_feature_permute.m:`` function used to permute ranks      
    
``sim_matrix_sameType_1.m``           : function used to compute original matrix             
 
``test_enhancer_names.m``             : function used to generate lists name for test

``uniform_bw_a_b.m ``                 : sample uniformly between a and b
        
   

*Requires Mtalab, RStudio, tidyverse and upsetR packages, and R*
