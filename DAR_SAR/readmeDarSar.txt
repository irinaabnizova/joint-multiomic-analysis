Compute_DARSAR

Computatitonal methods used to identify differentially accessible regions (DARs) and similarly accessible regions (SARs) geno-wide
from the GC-methylation data of bulked scNOME_seq.
This set of tools is used for producing of these two sets of chromatin accessible regions in bed-file format. 

Check the original manuscript where this method was applied for context details:

'Integrative analysis of transcriptomic and epigenomic data reveals distinct patterns for developmental and housekeeping genes'


Main scripts
compute_dar_sar_dominance_chrN.m : main script used for reading the GC-data (arranged per chromosome), computing it in a fixed size window,
filtering it by coverage, separately for each lineage. It further computes a dominance index for each window across the three 
lineages. It uses this dominance index to infer DARs and SAR as most dominant and even windows.
The user should specify the input variables: GC-files/path, window file/path and output path. Currently, it is hard-coded in the
beinning of the script.

Dependencies
test_data : contains examples of GC-methylation data for each lineage (pseudo-bulked) chr 13, precomputed coordinates of each 100bp window for 
chromosome 13 and expected outputs of DARs and SARs for chr13.

prefilter_count_chrN.m : function used to prefilter windows based on specified coverage threhold.

dar_sarHL_levGC_dominance.m : function used to compute a dominance index for each of prefiltered windows across lineages, and select
DARs and SARs based on it.

save_chromatin_EctEndMes.m : function used to write sets of DARs and SARs in bed format into specified output file/path

Requires MatlabR2014b
