# Joint multiomic analysis of mouse gastrulation data
![Fig6new_1](https://user-images.githubusercontent.com/61786710/132846059-046dbf8e-15ee-4c8d-bbe7-06892457ac45.png)


## Joint_multiomic_analysis: computational methods used to analyse jointly the transcriptome, chromatin accessibility and DNA methylation data from scNOME_seq data 



This set of methods allows to define two contrasting sets of genes: differentially expressed genes (DEGs) and similarly expressed genes (SEGs) from the transcriptomics data per lineage (pseudobulk).
It further identifies epigenetic patterns genome-wide, for both chromatin and DNA methylation layers.
It puts together gene expression and epigentic patterns, and establishes features separating developental (DEGs) vs housekeeping (SEGs) gene regulatory programs.
It analyses TF binding per lineage , depending on all three layers. It compares TF motif repertoires and infers pioneering and cooperative lineage-specific TF sets.

Check the original manuscript where this method was applied for context details:

'Integrative analysis of transcriptomic and epigenomic data reveals distinct patterns for developmental and housekeeping gene regulation'

Tools described here can potentially be extended to disentangle gene regulatory programs in other tissues.


 The tools are organised in several modules found in separate folders:
 
 - compute CG and GC counts in a 100 bp window
 
 - compute DARs, SARs: 
Computational methods used to identify differentially accessible regions (DARs) and similarly accessible regions (SARs) genome-wide
from the GC-methylation data of bulked scNOME_seq. This set of tools is used for producing of these two sets of chromatin accessible regions in bed-file format. 
 
 - compute DMRs, SMRs
 - compute DEG and SEGs
 - link chromatin and genes
 - link DNA methylation and genes
 - define lineage specific pioneering and cooperative TFs

- The data set is from (Argelaguet 2020). scNOME_seq tachnology as descibed in (Clark 2029, Kelly 2012)
Requirements: 
*Matlab'octave, C (for CG/GC methylation windows precomputing),RStudio, R*
