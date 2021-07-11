DEGs3_SEGs_Ch_High_data<-function(){
  
  #------------------ output with strand, ...thrE, set
  #-------------from other dir=data_genes/
  
  
  
  #-----------inputs internal hard-coded for each set
  #read Ect Mes Sim separately, merge them
  
  
  #-----------------output is unified
  #gg = splitted in three (DifEc DifM Sim) matrix, ect end mes

#print(param)
#param="Irina DEGs, Chastity adj"
  
  #-------------------New Diff High Genes: unique 
  #param5='Unique DEGs, high': already got 'meat'
  
  #-------------from other dir=data_genes/
  
  read_tsv("data_genes/DifEct_high.txt", col_types = cols(Chromosome = col_character())) -> DEct
  
  DEct%>%
    rename(chr=Chromosome, starts=Start, ends=End, strand=Strand)-> DEct
  
  read_tsv("data_genes/DifEnd_high.txt", col_types = cols(Chromosome = col_character())) -> DEnd
  DEnd%>%
    rename(chr=Chromosome, starts=Start, ends=End, strand=Strand)-> DEnd
  
  
  read_tsv("data_genes/DifMes_high.txt", col_types = cols(Chromosome = col_character())) -> DMes
  DMes%>%
  rename(chr=Chromosome, starts=Start, ends=End, strand=Strand)->DMes
  
  
  #--------------------------all genes Diff expressed, upregulated in one layer 
  
  #SEGs--------------sim (from Irina Chastity adjasted)
  source("readCh_adj_Sim_other.R")
  Sim<-readCh_adj_Sim_other()
  Sim%>%
    select(-Var8)%>%
    mutate(set="Sim")->Sim
  
  Sim%>%
    #mutate(strand1=(ifelse(strand==1)),'+','-')->Sim1
    mutate(strand1=case_when(strand==1 ~ "+",
                             strand==0 ~ "-")
    )%>%
    select(-strand)%>%
    select(gene,chr,starts,ends,strand1,ecto,endo,meso,set)%>%
    rename(strand=strand1)->Sim
    
 
  
  #source("densityGE_any.R")
  #pec<-densityGE_any(DEct)
  #pen<-densityGE_any(DEnd)
  #pm<-densityGE_any(DMes)
  #ps<-densityGE_any(Sim)
  
  #-------------------------
  #pc<-densityGE_any(Comm)
  
  #p<- grid.arrange(pc,pec,pen,pm,ps, ncol=1)
  #p
  
  #pga <-
   # ggarrange( pc,pec,pen,pm,ps, labels = c("Common","upEctoderm","upEndoderm", "upMesoderm","Similar"),
    #           common.legend = TRUE, legend = "bottom"
  #  )
  #pga
  #-------------------
  
  
  #pga <-
    #ggarrange( pec,pen,pm,ps, labels = c("upEctoderm","upEndoderm", "upMesoderm","Similar"),
    #           common.legend = TRUE, legend = "bottom"
   # )
  #pga
  


#---------------------------get lists
DEct %>% 
  distinct(gene) %>% 
  as_vector() -> DEct_vector

DEnd %>% 
  distinct(gene) %>% 
  as_vector() -> DEnd_vector

#len_DEn<-length(DEnd_vector)
#len_DEn

DMes %>% 
  distinct(gene) %>% 
  as_vector() -> DMes_vector

#------------------Sim
Sim %>% 
  distinct(gene) %>% 
  as_vector() -> Sim_vector


#concatenate DEGs SEGs as two- column matrix: gene set
#sh be EXACTLY same columns!!!
# A tibble: 515 x 9
#gene    chr     starts     ends strand   ecto   endo  meso  
#<chr>   <chr>    <dbl>    <dbl>  <dbl>  <dbl>  <dbl> <dbl> 
 # 1 Fam178b 1     36562694 36683183      0  1.87  -0.966 -2.14 
#2 Aff3    1     38177326 38664955      0  0.794  0.569 -2.04 
#3 Pgap1   1     54473000 54557684      0  2.29   0.252 -1.08 



EcEnMS=rbind(DEct,DEnd,DMes,Sim)

EcEnMS %>%
  mutate(thrE=-12)%>%
 select(-set)%>%
    mutate(set=case_when(gene %in% DMes_vector ~"DifM",
                          gene %in% DEct_vector ~"DifEc",
                          gene %in% DEnd_vector ~"DifEn",
                          gene %in% Sim_vector ~"Sim")
    )->separated_DEGs
#separated_DEGs
#------------------------split into matrices

diff<-separated_DEGs %>%
  group_by(set)

group_split(diff)->gg
#gg


#------------------------split into matrices

#gg

#--------------------------------then later on

#DifEct_ <-gg[[1]]
#DifEct_

#DifMes_ <-gg[[2]]
#DifMes_

#DifEnd_ <-gg[[3]]

#Sim_ <-gg[[4]]

}