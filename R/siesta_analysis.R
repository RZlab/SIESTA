siesta_analysis <- function(results, treatment, vehicle){
  require(broom)
  require(tidyverse)

  results.t <- results %>%
    filter(term == "Tm") %>%
    dplyr::select(-c(term, std.error:rSquared)) %>%
    separate(Sample, sep = "_", remove = F, c("Cell_line", "Treatment", "Rep"))
  
  ggplot(results.t) +
    geom_density(aes(x = estimate, fill = Sample), alpha = 0.2) +
    theme_minimal()
  
  #ES vs. (E, S, C)
  rES <- results.t %>%
    group_by(id) %>%
    do({
      t <- .
      pid <- unique(t$id)
      
      ES <- t %>%
        filter(Treatment == "ATPEnz")
      
      nES <- t %>%
        filter(Treatment != "ATPEnz")
      
      
      if(nrow(ES) >= 2 & nrow(nES) >= 2){
        t.test <- try(t.test(ES$estimate, nES$estimate, var.equal = T)$p.value, silent = T)
        
        if(is.numeric(t.test)){
          res <- data.frame("id" = pid,
                            mean.ES = mean(ES$estimate, na.rm = T),
                            mean.nES = mean(nES$estimate, na.rm = T),
                            p.value = t.test)
        }
      }
      res
    })
  write_tsv(rES, "ES_vs_nES.tsv")
  
  #E vs. (S, C)
  rE <- results.t %>%
    filter(Treatment != "ATPEnz") %>%
    group_by(id) %>%
    do({
      t <- .
      pid <- unique(t$id)
      
      E <- t %>%
        filter(Treatment == "Enz")
      
      nE <- t %>%
        filter(Treatment != "Enz")
      
      
      if(nrow(E) >= 2 & nrow(nE) >= 2){
        t.test <- try(t.test(E$estimate, nE$estimate, var.equal = T)$p.value, silent = T)
        
        if(is.numeric(t.test)){
          res <- data.frame("id" = pid,
                            mean.E = mean(E$estimate, na.rm = T),
                            mean.nE = mean(nE$estimate, na.rm = T),
                            p.value = t.test)
        }
      }
      res
    })
  write_tsv(rE, "E_vs_nE.tsv")
  
  #S vs. (C)
  rS <- results.t %>%
    filter(Treatment == "ATP" | Treatment == "C" ) %>%
    group_by(id) %>%
    do({
      t <- .
      pid <- unique(t$id)
      
      S <- t %>%
        filter(Treatment == "ATP")
      
      nS <- t %>%
        filter(Treatment != "ATP")
      
      
      if(nrow(S) >= 2 & nrow(nS) >= 2){
        t.test <- try(t.test(S$estimate, nS$estimate, var.equal = T)$p.value, silent = T)
        
        if(is.numeric(t.test)){
          res <- data.frame("id" = pid,
                            mean.S = mean(S$estimate, na.rm = T),
                            mean.nS = mean(nS$estimate, na.rm = T),
                            p.value = t.test)
        }
      }
      res
    })
  write_tsv(rS, "E_vs_nE.tsv")
  
  # plot S vs C
  p <- results.t %>%
    group_by(id, Treatment) %>%
    summarise(mean.estimate = mean(estimate, na.rm = T))
  
  t <- p %>%
    filter(Treatment == "ATP" | Treatment == "C" ) %>%
    spread(Treatment, mean.estimate) %>%
    drop_na()
  
  ggplot(t) +
    geom_point(aes(x = C, y = ATP), alpha = 0.5) +
    theme_minimal()
  
  
  # plot ES-E vs ES-S
  t <- p %>%
    filter(Treatment != "C") %>%
    spread(Treatment, mean.estimate) %>%
    mutate("ES-E" = ATPEnz - Enz) %>%
    mutate("ES-S" = ATPEnz - ATP) %>%
    drop_na()
  
  ggplot(t) +
    geom_point(aes(x = `ES-E`, y = `ES-S`), alpha = 0.5) +
    theme_minimal() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0)

}
