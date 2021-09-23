## Module to compute distances between the nearest (expressed ) genes

![Fig1](https://user-images.githubusercontent.com/61786710/134468637-48e09cd7-ef96-4a58-a393-3fc257dd2ba0.png)

Here one can find a set of computatitonal tools used to to compute distances between the nearest (expressed ) genes

The detailed description of the method in in the original manuscript:
> 'Integrative analysis of transcriptomic and epigenomic data reveals distinct patterns for developmental and housekeeping genes'.

#### Main scripts and functions
``DEG_SEG_distances_1.Rmd`` : 


**usage** 

To run it, open the code in RStudio and click on running RStudio options:
``DEG_SEG_distances_1.Rmd``

The standard  input files and the names/path of ouput files are supplied in the gene data folder:

#### gene data

``gene_data`` folder contains sets of Differentially Expressed Genes (DEGs) for each lineage, and of Similarly Expressed Genes across the three lineages (SEGs)
as well as pseudo-bulk data for the transcriptomics:
```
DifEct_high.txt                                 
DifEnd_high.txt                                 
DifMes_high.txt                                 
E7_5_pseudobulk_log2CPM_200331.txt              
RPKM_Ect_Mes_19117.txt                          
log2RPKM_filtered_for_gene_body_annotation.txt  
sim_par_thrn12.txt                              
tab_DifChEct_3.txt                              
tab_DifChEnd_3.txt                              
tab_DifChMes_3.txt                              
tab_SimCh_3.txt                                 
tab_SimCh_3_add.txt   

```

#### Dependencies

There are functions reading data, computing stats and visualisation:
  
```
DEGs3_SEGs_Ch_High_data.R      
DEGs_SEGs_Chadj.R              
DEGs_SEGs_IrChadj.R            
data_genes                     
density_DifEct.R               
density_DifEnd.R               
density_DifMes.R               
density_Sim.R                  
distanceTSS_ect.R              
distanceTSS_end.R              
distanceTSS_mes.R              
distanceTSS_sim.R              
len_DEG_SEG.R                  
readCh_adj_Ect.R               
readCh_adj_End.R               
readCh_adj_Mes.R               
readCh_adj_Sim.R               
readCh_adj_Sim_other.R         
read_NPSB.R                    
read_filter_expr_Chr_recent.R  
summ_med_dist_count.R          
vEct_Len.R                     
vMes_Len.R                     
violin5_dist.R                 

```
   

*Requires RStudio, tidyverse package, and R*
