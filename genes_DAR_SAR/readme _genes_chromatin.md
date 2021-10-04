### genes DAR SAR module contains functions which link differentially and similarly accessible chromatin and genes
![TSS_dar_sar](https://user-images.githubusercontent.com/61786710/135826676-31351141-2882-4991-8c58-4c34ff4b3ee0.png)
The detailed description of the method in in the original manuscript:
> 'Integrative analysis of transcriptomic and epigenomic data reveals distinct patterns for developmental and housekeeping genes'.

#### Main script and functions
``dar_sar_around_DEG_SEG_chr1.m`` : 


**usage** 

To run it, open the code in Matlab/octave and click run command:
``dar_sar_around_DEG_SEG_chr1.m``

**dependencies**
```
around_target_dis2start.m  
enhs_targets_stats_chr_sorted.m  
histo5K_normVicNov.m             
nearest2targets_dis2start.m      
separate_chr.m                   
smooth_histo_TSS_fig.m           
smooth_win3_float.m              
smooth_win5_float.m              
test_data
```

#### test data

``test_data`` folder contains sets of Differentially Expressed Genes (DEGs) for each lineage chr1, and of Similarly Expressed Genes across the three lineages (SEGs) chr1,
together with DARs and sARs chr1:

```
DAR_chr1.txt  
DEG_chr1.txt  
SAR_chr1.txt  
SEG_chr1.txt  
```

*Requires Matlab or octave*
