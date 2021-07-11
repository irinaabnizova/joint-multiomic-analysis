readCh_adj_Sim<-function(){
read_tsv('data_genes/tab_SimCh_3.txt', col_types = cols(chr = col_character())) -> DifEct_Ir
#tab_SimCh_3
DifEct_Ir %>%
  rename(gene=Row)->DifEct_Ir
DifEct_Ir
}
