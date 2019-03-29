siesta_analysis <- function(results, treatment, vehicle){
 require(broom)
  require(tidyverse)

  results.t <- results %>%
    filter(rSquared >= 0.95) %>%
    filter(term == "Tm") %>%
    dplyr::select(-c(term, std.error:rSquared)) %>%
    separate(Sample, sep = "_", remove = F, c("Cell_line", "Treatment", "Rep"))
  
  g <- ggplot(results.t) +
    geom_density(aes(x = estimate, fill = Sample), alpha = 0.2) +
    facet_wrap(~Treatment, ncol = 2) +
    theme_minimal() +
    theme(legend.position = "bottom")
 ggsave(g, "Tm_distribution.pdf")
  
  #CoFaEnzy vs. (Enzy, CoFa, Cntrl)
  rES <- results.t %>%
    group_by(id) %>%
    do({
      t <- .
      pid <- unique(t$id)
      
      ES <- t %>% filter(Treatment == "CoFaEnzy")
      nES <- t %>% filter(Treatment != "CoFaEnzy")
            
      if(nrow(ES) >= 2 & nrow(nES) >= 2){
        t.test <- try(t.test(ES$estimate, nES$estimate, var.equal = T)$p.value, silent = T)
        
        if(is.numeric(t.test)){
          res <- data.frame("id" = pid,
                            mean.CoFaEnzy = mean(ES$estimate, na.rm = T),
                            mean.nCoFaEnzy = mean(nES$estimate, na.rm = T),
                            p.value = t.test)
        }
       else{
       res <- data.frame()
       }
      }
     else{
       res <- data.frame()
       }
      res
    })
  write_tsv(rES, "CoFaEnzy_vs_nCoFaEnzy.tsv")
  
  #Enzy vs. (CoFa, Cntrl)
  rE <- results.t %>%
    filter(Treatment != "SubEnz") %>%
    group_by(id) %>%
    do({
      t <- .
      pid <- unique(t$id)
      
      E <- t %>% filter(Treatment == "Enzy")
      nE <- t %>% filter(Treatment != "Enzy")
      
      if(nrow(E) >= 2 & nrow(nE) >= 2){
        t.test <- try(t.test(E$estimate, nE$estimate, var.equal = T)$p.value, silent = T)
        
        if(is.numeric(t.test)){
          res <- data.frame("id" = pid,
                            mean.Enzy = mean(E$estimate, na.rm = T),
                            mean.nEnzy = mean(nE$estimate, na.rm = T),
                            p.value = t.test)
        }
       else{
       res <- data.frame()
       }
      }
     else{
       res <- data.frame()
       }
      res
    })
  write_tsv(rE, "Enzy_vs_nEnzy.tsv")
  
  #CoFa vs. (Enzy, Cntrl)
  rS <- results.t %>%
    filter(Treatment != "SubEnz") %>%
    group_by(id) %>%
    do({
      t <- .
      pid <- unique(t$id)
      
      S <- t %>% filter(Treatment == "CoFa")
      nS <- t %>% filter(Treatment != "CoFa")
      
      if(nrow(S) >= 2 & nrow(nS) >= 2){
        t.test <- try(t.test(S$estimate, nS$estimate, var.equal = T)$p.value, silent = T)
        
        if(is.numeric(t.test)){
          res <- data.frame("id" = pid,
                            mean.CoFa = mean(S$estimate, na.rm = T),
                            mean.nCoFa = mean(nS$estimate, na.rm = T),
                            p.value = t.test)
        }
       else{
       res <- data.frame()
       }
      }
     else{
       res <- data.frame()
       }
      res
    })
  write_tsv(rS, "CoFa_vs_nCoFa.tsv")
  
  # plot CoFa vs Cntrl
  p <- results.t %>%
    group_by(id, Treatment) %>%
    summarise(mean.estimate = mean(estimate, na.rm = T))
  
  t <- p %>%
    filter(Treatment == "CoFa" | Treatment == "Cntrl" ) %>%
    spread(Treatment, mean.estimate) %>%
    drop_na()
  
  g <- ggplot(t) +
    geom_point(aes(x = Cntrl, y = CoFa), alpha = 0.5) +
    theme_minimal()
  ggsave(g, "CoFa_vs_Cntrl.pdf")
  
  # plot ES-E vs ES-S
  t <- p %>%
    filter(Treatment != "Cntrl") %>%
    spread(Treatment, mean.estimate) %>%
    mutate("CoFaEnzy-Enzy" = CoFaEnzy - Enzy) %>%
    mutate("CoFaEnzy-CoFa" = CoFaEnzy - CoFa) %>%
    drop_na()
  
  g <- ggplot(t) +
    geom_point(aes(x = `CoFaEnzy-Enzy`, y = `CoFaEnzy-CoFa`), alpha = 0.5) +
    theme_minimal() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0)
 ggsave(g, "siesta_plot.pdf")
}
