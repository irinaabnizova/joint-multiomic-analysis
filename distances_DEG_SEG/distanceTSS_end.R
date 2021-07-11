distanceTSS_end <- function(main_raw, thrGE){

  #input is a tibble
  #gene
  #<chr>
  #  Chromosome
  #<chr>
  #  Start
 # <dbl>
  #  End
  #<dbl>
  #  Strand
  #<chr>
  #GEect
  #GEend
  #GEmes
  
  #Xkr4	1	3205901	3671498	-	-4.350450000	-2.398244100	-3.207620000
  #Rp1	1	4343507	4360314	-	-9.973417000	-8.495143000	-3.278361600
  #Sox17	1	4490928	4496413	-	-4.482705600	6.171816300	-4.565535500
  #Mrpl15	1	4773206	4785739	-	4.971182000	4.610790000	4.919617700
  #Lypla1	1	4807788	4886770	+	4.222963300	4.634207700	4.020493000
  #Tcea1	1	4857814	4897909	+	5.305976400	4.969354000	5.186199700
  
  
  #0------filter due to ect >thrGE
  #1.1  flip REV strand
  #create swapped  columns for REV strand, StartSW EndSW
  
  
  main_raw %>%
    filter(Endo > thrGE) %>%
    mutate(StartSW = ifelse(Strand == "-", End, Start), 
    EndSW = ifelse(Strand =="-", Start, End)) %>% 
    select(gene, Chromosome,StartSW, EndSW, Strand,Ecto,Endo,Meso)-> main_sw
   main_sw

#2. introduce addiational columns by shift functions

#source("shift_next.R")
#source("shift_prev.R")
shift_next <- function(col1){
  #next
  col1 %>% lead()->col2
  n=length(col2)
  # n
  col2[n]=col2[n-2]#last gene just repeats n-2 start (to compute dist further)
  return(col2)
}

shift_prev <- function(col1){
  #previous
  col1 %>% lag()->col2
  col2[1]=col1[1]
  # col2[n]=col2[n-1]#last gene just repeats n-2 start (to compute dist further)
  return(col2)
}


#2.1  select only needed columns: gene,Chromosome,StartSW,EndSW, Strand 
#AND
#- compute distance bw shifted starts(just righ to left distance) and right neighbour gene(gene_dist)
#AND
#create new columns to compute and vizualise new tible with min-dist

main_sw %>%
  select(gene,Chromosome,StartSW,EndSW, Strand,Ecto,Endo,Meso) %>%
  group_by(Chromosome)%>%
  mutate(StartN=shift_next(StartSW)) %>%
  mutate(geneN=shift_next(gene)) %>%
  mutate(distN=abs(StartN - StartSW)) %>%
  mutate(distPr=shift_prev(distN)) %>%
  mutate(genePr=shift_prev(gene))%>%
  mutate(dif_dist=(distN-distPr))%>%
  mutate(min_dist=if_else(dif_dist >= 0, distPr,distN))%>%
  mutate(min_gene=if_else(dif_dist >= 0, genePr,geneN))%>%
  mutate(gene1=gene)->mi_dist_gene
mi_dist_gene 

#2.2-correlate min and just distance
#pcd<- ggplot(mi_dist_gene,aes(distN,min_dist))+geom_point(colour="magenta", size=1)+ ggtitle("min distance bw genes vs just distance")
#pcd

#3 - compute Length column, and rename template file

mi_dist_gene %>%
  select(gene,Chromosome,StartSW,EndSW, Strand, min_dist,min_gene,Ecto,Endo,Meso) %>%
  rename(chr=Chromosome, Start=StartSW, End=EndSW)%>%
  mutate(Len=abs(End-Start))%>%
  ungroup()-> gene_start_dist
gene_start_dist



}