violin5_dist<- function(separated_CDS, med_Dif_EcMSm,param2,thrGE){
  #input: long format for set type: DifEct vs CommEct
  #param1="Christel or Carine DEGs"
  
 
    #ceiling(med_dist_DifEct/1000)->med_dist_DifEct
  
  ceiling(med_Dif_EcMSm$median_lin/1000)->med_dist_DifEct
  #med_dist_DifEct
  med_Dif_EcMSm$count->count_DifEct
  #count_DifMes
  
  separated_CDS %>% 
    ggplot(aes(set, log2(min_dist), fill = set)) +
    geom_violin(draw_quantiles = 0.5) +
    theme_bw(base_size = 16) +
    #scale_fill_manual(values = c("green","pink","yellow","yellow")) +
    scale_fill_manual(values = c("green","cyan","magenta","orange","white","yellow","cyan","magenta","yellow","black")) +
    ylab("dist2nearest expr gene (log2 bp)\n") +
    xlab(element_blank()) +
    theme(axis.title.y = element_text(size = 14), plot.title = element_text(size = 14, hjust = 0.5)) + 
    annotate(geom = "text", x = 1, y = 17, label = paste(med_dist_DifEct[1], "kb")) +
    annotate(geom = "text", x = 2, y = 17, label = paste(med_dist_DifEct[2], "kb")) +
    annotate(geom = "text", x = 3, y = 17, label = paste(med_dist_DifEct[3], "kb")) +
    annotate(geom = "text", x = 4, y = 17, label = paste(med_dist_DifEct[4], "kb")) +
    annotate(geom = "text", x = 5, y = 17, label = paste(med_dist_DifEct[5], "kb")) +
    
    annotate(geom = "text", x = 2, y = 25, label = paste(param2)) +
    
    annotate(geom = "text", x = 1, y = 2, label = paste("n=",count_DifEct[1])) +
    annotate(geom = "text", x = 2, y = 2, label = paste("n=",count_DifEct[2])) +
    annotate(geom = "text", x = 3, y = 2, label = paste("n=",count_DifEct[3])) +
    annotate(geom = "text", x = 4, y = 2, label = paste("n=",count_DifEct[4])) +
    annotate(geom = "text", x = 5, y = 2, label = paste("n=",count_DifEct[5])) +
    
    #annotate(geom = "text", x = 2.5, y = 25, label = paste(thrGE)) +
    #annotate(geom = "text", x = 2, y = 23, label = paste("thrGE=",thrGE)) +
    ggtitle("dist for DEGs SEGs expressed in one layer") -> violin1.any
  
   #violin1.ect
 
}