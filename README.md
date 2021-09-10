# joint_multiomic_analysis 
C:\Users\ia1\Documents\irina_docs\EpiGeVar\JointDifPap_currSept\final_figures_text\rene_figs12Aug\Fig6new_1.png 

joint_multiomic_analysis :: Computational methods used to analyse jointly the transcriptome, chromatin accessibility and DNA methylation data from scNOME_seq data (Clark 2018).

This set of methods allows to define two contrasting sets of genes: differentially expressed genes (DEGs) and similarly expressed genes (SEGs) from the transcriptomics data per lineage (pseudobulk).
It further identifies epigenetic patterns genome-wide, for both chromatin and DNA methylation layers.
It puts together gene expression and epigentic patterns, and establishes features separating developental (DEGs) vs housekeeping (SEGs) gene regulatory programs.
It analyses TF binding per lineage , depending on all three layers. It compares TF motif repertoires and infers pioneering and cooperative lineage-specific TF sets.

Check the original manuscript where this method was applied for context details:

'Integrative analysis of transcriptomic and epigenomic data reveals distinct patterns for developmental and housekeeping genes'

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

Requirements: 
Matlab,
C (for CG/GC methylation windows precomputing),
R tidyverse
