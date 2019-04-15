siesta_analysis <- function(results, treatment, vehicle, rSquared.filter =  0.95){
 require(broom)
 require(tidyverse)

 results_export <- results %>%
    filter(term == "Tm") %>%
    dplyr::select(-c(term, std.error:rSquared)) %>%
    mutate("Treatment" = gsub("\\_.*", "", Sample))
  write_tsv(results_export, "all_results_export.tsv")
 
  results_export_opls <- results_export %>%
    select(-Treatment) %>%
    spread(Sample, estimate)
   write_tsv(results_export_opls, "all_results_export_opls.tsv")
 
  results.t <- results %>%
    filter(rSquared >= rSquared.filter) %>%
    filter(term == "Tm") %>%
    dplyr::select(-c(term, std.error:rSquared)) %>%
    separate(Sample, sep = "_", remove = F, c("Cell_line", "Treatment", "Rep"))
 
  ggplot(results.t) +
    geom_density(aes(x = estimate, fill = Sample), alpha = 0.2) +
    facet_wrap(~Treatment, ncol = 2) +
    theme_minimal() +
    theme(legend.position = "bottom")
 ggsave(paste("Tm_distribution_", rSquared.filter, ".pdf", sep = "_"))
  
  #CoFaEnzy vs. (Enzy, CoFa, CTRL)
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
  write_tsv(rES,  paste("CoFaEnzy_vs_nCoFaEnzy", rSquared.filter, ".tsv", sep = "_"))
 
  #Enzy vs. (CoFa, CTRL)
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
  write_tsv(rE,  paste("Enzy_vs_nEnzy", rSquared.filter, ".tsv", sep = "_"))
 
  #CoFa vs. (Enzy, CTRL)
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
  write_tsv(rS,  paste("CoFa_vs_nCoFa", rSquared.filter, ".tsv", sep = "_"))
  # plot CoFa vs CTRL
  p <- results.t %>%
    group_by(id, Treatment) %>%
    summarise(mean.estimate = mean(estimate, na.rm = T))
  
  t <- p %>%
    filter(Treatment == "CoFa" | Treatment == "CTRL" ) %>%
    spread(Treatment, mean.estimate) %>%
    drop_na()
  
  ggplot(t) +
    geom_point(aes(x = CTRL, y = CoFa), alpha = 0.5) +
    theme_minimal()
  ggsave(paste("CoFa_vs_CTRL", rSquared.filter, ".pdf", sep = "_"))
 
  ggplot(rS) +
    geom_point(aes(x = log2(mean.CoFa/mean.nCoFa), y = -log10(p.value), id = id), alpha = 0.5) +
    theme_minimal()
   ggsave(paste("CoFa_vs_CTRL_vulcano", rSquared.filter, ".pdf", sep = "_"))
 
  # plot ES-E vs ES-S
  t <- p %>%
    filter(Treatment != "CTRL") %>%
    spread(Treatment, mean.estimate) %>%
    mutate("CoFaEnzy-Enzy" = CoFaEnzy - Enzy) %>%
    mutate("CoFaEnzy-CoFa" = CoFaEnzy - CoFa) %>%
    drop_na()
  
  ggplot(t) +
    geom_point(aes(x = `CoFaEnzy-Enzy`, y = `CoFaEnzy-CoFa`), alpha = 0.5) +
    theme_minimal() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0)
 ggsave(paste("siesta_plot", rSquared.filter, ".pdf", sep = "_"))
}
