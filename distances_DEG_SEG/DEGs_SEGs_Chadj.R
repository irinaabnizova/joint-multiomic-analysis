DEGs_SEGs_Chadj<-function(param){
  
  #-----------------output is
  #gg = splitted in three (DifEc DifM Sim) matrix

#2.1. DeGs Unique Ect Mes Irina
print(param)
#param2="Irina DEGs, Chastity adj"

#2.1--------------------------MZ unique Diff genes

#2.2-----------------Z unique Diff genes

#2.4-----------------MZ Ect vs Mes Diff_gene_MesEctZF

#2.3-----------------Z Ect vs Mes Diff_gene_MesEctZF

#2.5---------------------Ch_3 (by shifting up depending on GE, unique)
source("readCh_adj_Ect.R")
DifEct_Ir<-readCh_adj_Ect()
source("readCh_adj_End.R")
DifEnd_Ir<-readCh_adj_End()
source("readCh_adj_Mes.R")
DifMes_Ir<-readCh_adj_Mes()


#--------------sim
source("readCh_adj_Sim.R")
Sim_Ir<-readCh_adj_Sim()




#2.0 #-----------------visualise param2

source("density_DifEct.R")
pe<-density_DifEct(DifEct_Ir)
#pe
source("density_DifEnd.R")
pen<-density_DifEct(DifEnd_Ir)

source("density_DifMes.R")
pm<-density_DifMes(DifMes_Ir)
#pm

source("density_Sim.R")
ps<-density_Sim(Sim_Ir)
#ps

#ggarrange(pe,pm, ncols=1)
grid.arrange(pe,pen,pm,ps)


#get lists
DifEct_Ir %>% 
  distinct(gene) %>% 
  as_vector() -> DEct_vector

DifEnd_Ir %>% 
  distinct(gene) %>% 
  as_vector() -> DEnd_vector

DifMes_Ir %>% 
  distinct(gene) %>% 
  as_vector() -> DMes_vector

#------------------Sim
Sim_Ir %>% 
  distinct(gene) %>% 
  as_vector() -> Sim_vector


#concatenate DEGs SEGs as two- column matrix: gene set
EcEnMS=rbind(DifEct_Ir,DifEnd_Ir,DifMes_Ir,Sim_Ir)
EcEnMS %>%
 select(-Var8)%>%
    mutate(type=case_when(gene %in% DMes_vector ~"DifM",
                          gene %in% DEct_vector ~"DifEc",
                          gene %in% DEnd_vector ~"DifEn",
                          gene %in% Sim_vector ~"Sim")
    )->separated_DEGs

#------------------------split into matrices

diff<-separated_DEGs %>%
  group_by(type)

group_split(diff)->gg
#gg

#--------------------------------then later on

#DifEct_ <-gg[[1]]
#DifEct_

#DifMes_ <-gg[[2]]
#DifMes_

#Sim_ <-gg[[3]]

}