summ_med_dist_count<-function(separated_DS_any){
  
  #sh be column 'set' and min_dist; any= any_layer
  
  separated_DS_any %>% 
    group_by(set) %>% 
    summarise(median_lin = median(min_dist),
              count=n()) -> set_med_n
  #medians_lin.mes
  
  #separated_DS_mes %>% 
   # group_by(set) %>% 
   # summarise(count=n()) ->set_counts_mes
  #set_counts_mes    
  
  #medians_lin.mes %>%
    #inner_join(set_counts_mes) -> set_med_n_mes
  #set_med_n_mes
  
}