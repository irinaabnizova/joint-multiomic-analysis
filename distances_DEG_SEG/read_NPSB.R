read_NPSB<-function(param){
  
  #--------------  _other point to data in the main location
  #C:/Users/ia1/Documents/R_work/Ect_Mes_DEG/E7_5_pseudobulk_log2CPM_200331.txt

  
  #--------------------read and filter new pseudobulk
  print(param)
  
  #param3='Christel new PSP'
  
  
  
  read_tsv("data_genes/E7_5_pseudobulk_log2CPM_200331.txt", col_types = cols(Chromosome = col_character())) -> expressionRaw
  #expressionRaw
  
  expressionRaw %>% 
    filter(Chromosome != "ERCC" & Chromosome != "MT" & Chromosome != "Y") %>% 
    select(Probe, Chromosome, Start, End, `Probe Strand`, Ectoderm_counts, Endoderm_counts, Mesoderm_counts) %>% 
    distinct(Probe, .keep_all = T) %>% 
    rename(gene = Probe, Strand = 'Probe Strand', Ecto = Ectoderm_counts, Endo = Endoderm_counts, Meso = Mesoderm_counts) -> main_raw
  #main_raw
  
  #dim (main_raw)->size_mr
  #size_mr
  
  
}