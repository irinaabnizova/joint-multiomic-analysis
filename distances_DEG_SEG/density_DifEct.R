density_DifEct <- function(DifEct){
  
  DifEct %>%
      gather(key=layer, value=GE,ecto,endo,meso) -> dfect
      tail(dfect)
      dim(dfect)
                           pe1<-ggplot(dfect,aes (GE))+
                             geom_density(aes (fill=layer), color ="black", alpha=0.5)+
                             scale_fill_manual(values=c("green", "blue", "magenta"))+
                             #ggtitle("Gene expression distribution each layer")+
                             #ylim(0, 0.3)+
                             theme_minimal(base_size = 14)+
                             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                             xlab("Gene expression, upEct")+
                             ylab("density")
                           pe1
}