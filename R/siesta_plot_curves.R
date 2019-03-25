siesta_plot_curves <- function(Curves = curves,
           Results = results,
           cores = n_cores){
    library("foreach")
    library("doSNOW")
    library("dplyr")
    
    no_cores <- cores
    cl <- makeCluster(no_cores)
    registerDoSNOW(cl)  
    message(getDoParWorkers())
    
    id <- unique(Curves$id)
    
    foreach(i = 1:length(id), .packages = c("ggplot2", "gridExtra", "dplyr", "tidyr")) %dopar% {
      t.curves <- Curves[which(Curves$id == id[i]),]
      t.results <- Results[which(Results$id == id[i]),]
      
      if(nrow(t.curves) > 0 & nrow(t.results) > 0){
        gg1 <-
          ggplot()  +
          geom_line(data = t.curves, aes(x = Temperature, color = Type, group = Sample, y = prValue), alpha = 0.8) +
          geom_point(data = t.curves, aes(x = Temperature, color = Type, group = Sample, y = Value), alpha = 0.3, size = 1) +
          theme_classic()+
          labs(x = "Temperature",
               y = "Value",
               title = t.curves$`Gene names`) +
          theme(legend.position = "bottom")
        
        suppressMessages(gg2 <- tableGrob(
          t.results %>% mutate(
            V = sprintf("%.2f", estimate),
            se = sprintf("%.2f", std.error),
            R2 = sprintf("%.2f", rSquared)
          )  %>% filter(term == 'Tm') %>%
            select(Sample, V, se) %>%
            unite(Val, V, se, sep="±") %>%
            spread(Sample, Val), theme = ttheme_minimal(base_size = 5), rows=""))
        
        m1 <- arrangeGrob(gg1, gg2, heights=c(13, 1))
        
        ggsave(paste("fit_", id[i], ".pdf", sep = ""), plot = m1, device = "pdf", width = 7, height = 7)
      }
    }
    stopCluster(cl)
  }
