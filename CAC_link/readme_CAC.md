
### This directory contains files computing CAC(Chromatin Abundance Coefficient), which is used to link putative enhancers with their target genes![Fig3n_2](https://user-images.githubusercontent.com/61786710/151939059-cd2d04e3-2a4f-4dde-9dc0-f4723eaf2d41.png)
![Fig3n_2](https://user-images.githubusercontent.com/61786710/151939109-88dba9d3-551c-4119-8013-7e4208a2dc75.png)

- To compute CAC, the number of open chromatin windows (DARs or SARs) was calculated across the TSS vicinity Vi for each gene gi in a set first, N(Vi). Then the number is normalised by the count of expressed genes in this vicinity, Ng(Vi) (normalised accessible region frequency):
                                   
                                   normalised AR frequency in Vi = N(Vi)  / Ng(Vi) 

- It is then compared with the average gene expression (within the same normalised frequency counts) across genes having GEi  , by computing a Pearson correlation coefficient:
                      	
                        CAC = Pearson correlation (N(Vi)  / Ng(Vi) , mean(GEi ))

- If the CAC value is positive and high  (R2>0.7) we define the corresponding sets of genes and DARs as linked. 
