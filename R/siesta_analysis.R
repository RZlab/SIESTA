siesta_analysis <- function(results = results, rSquared.filter =  0.7){
  library(tidyverse)
  library(plotly)
  library(ggrepel)
  library(htmlwidgets)
  library(broom)
  
  results <- results %>%
    filter(rSquared > rSquared.filter) %>%
    filter(term == "Tm") %>%
    filter(!grepl("CTRL", Sample))
  
  results.t <- results %>%
    filter(term == "Tm") %>%
    dplyr::select(-c(term:rSquared)) %>%
    separate(Sample, sep = "_", remove = F, c("Cell_line", "Treatment", "Rep"))
  
  n <- results.t %>%
    group_by(id) %>%
    summarise(n = n()) %>%
    filter(n == 6)
  
  t <- results.t %>%
    filter(id %in% n$id) %>%
    filter(Treatment != "CTRL") %>%
    group_by(id, `Gene names`, Treatment) %>%
    summarise(mean.estimate = mean(Tm_root, na.rm = T),
              diff.estimate = diff(Tm_root),
              sd.estimate = sd(Tm_root, na.rm = T)) %>%
    filter(abs(sd.estimate) < 2.5)
  
  tpp <- t %>%
    group_by(id, `Gene names`) %>%
    summarise(n = n()) %>%
    filter(n == 3)
  
  t <- t %>%
    filter(id %in% tpp$id)
  
  oplsda <- results %>%
    filter(id %in% t$id) %>%
    select(-c(term:rSquared)) %>%
    spread(Sample, Tm_root)
  write_tsv(oplsda, "export_opls.tsv")
  
  d <- results.t %>%
    filter(id %in% t$id)

  rCoFa <- d %>%
    filter(Treatment != 'Enzy') %>%
    group_by(id) %>%
    do({
      t <- .
      pid <- unique(t$id)
      
      ES <- t %>% filter(Treatment == "CoFaEnzy")
      S <- t %>% filter(Treatment == "CoFa")
      
      if(nrow(ES) >= 2 & nrow(S) >= 2){
        t.test <- try(t.test(ES$Tm_root, S$Tm_root, var.equal = T)$p.value, silent = T)
        
        if(is.numeric(t.test)){
          res <- data.frame("id" = pid,
                            `Gene names` = t$`Gene names`[1],
                            mean.CoFaEnzy = mean(ES$Tm_root, na.rm = T),
                            mean.CoFa = mean(S$Tm_root, na.rm = T),
                            p.value_CoFa = t.test)
        }
        else{
          res <- data.frame()
        }
      }
      else{
        res <- data.frame()
      }
      res
    }) %>%
    mutate(Tm_diff_CoFa = (mean.CoFaEnzy - mean.CoFa)) %>%
    mutate(sig_vs_CoFa = ifelse(p.value_CoFa < 0.05 & abs(Tm_diff_CoFa) > 1, 1, 0))

  rEnzy <- d %>%
    filter(Treatment != 'CoFa') %>%
    group_by(id) %>%
    do({
      t <- .
      pid <- unique(t$id)
      
      ES <- t %>% filter(Treatment == "CoFaEnzy")
      E <- t %>% filter(Treatment == "Enzy")
      
      if(nrow(ES) >= 2 & nrow(E) >= 2){
        t.test <- try(t.test(ES$Tm_root, E$Tm_root, var.equal = T)$p.value, silent = T)
        
        if(is.numeric(t.test)){
          res <- data.frame("id" = pid,
                            `Gene names` = t$`Gene names`[1],
                            mean.CoFaEnzy = mean(ES$Tm_root, na.rm = T),
                            mean.Enzy = mean(E$Tm_root, na.rm = T),
                            p.value_Enzy = t.test)
        }
        else{
          res <- data.frame()
        }
      }
      else{
        res <- data.frame()
      }
      res
    }) %>%
    mutate(Tm_diff_Enzy = (mean.CoFaEnzy - mean.Enzy)) %>%
    mutate(sig_vs_Enzy = ifelse(p.value_Enzy < 0.05 & abs(Tm_diff_Enzy) > 1, 1, 0))

  dat <- rCoFa %>%
    full_join(rEnzy, by = c('id', 'Gene.names')) %>%
    mutate(SIESTA = ifelse((p.value_CoFa < 0.05 & p.value_Enzy < 0.1 & sign(Tm_diff_CoFa) == sign(Tm_diff_Enzy) & abs(Tm_diff_CoFa) > 1 & abs(Tm_diff_Enzy) > 1) |
                             (p.value_CoFa < 0.1 & p.value_Enzy < 0.05 & sign(Tm_diff_CoFa) == sign(Tm_diff_Enzy) & abs(Tm_diff_CoFa) > 1 & abs(Tm_diff_Enzy) > 1),
                           1, 0))
  write_tsv(dat, '2D_SIESTA.txt')
  
  
  g <- ggplot(dat, aes(x = Tm_diff_CoFa, y = Tm_diff_Enzy, colour = as.factor(SIESTA))) +
    geom_point(alpha = 0.5, size = 2, shape = 16) +
    theme_minimal() +
    scale_color_manual(values=c("#999999", "#00BA38")) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    theme(legend.position = 'bottom') +
    xlim(-10, 10) +
    ylim(-8, 8) +
    geom_text_repel(data = subset(dat, SIESTA == 1),
                    aes(label = paste(`Gene.names`)), size = 3, box.padding = unit(0.35, 'lines'),
                    point.padding = unit(0.3, 'lines'))
  ggsave('2D_SIESTA.pdf', width = 8, height = 8, units = "in")
  
  gp <- ggplotly(ggplot(dat, aes(x = Tm_diff_CoFa, y = Tm_diff_Enzy, colour = as.factor(SIESTA), id = Gene.names)) +
                   geom_point(alpha = 0.5) +
                   theme_minimal() +
                   geom_hline(yintercept = 0, linetype = 'dashed') +
                   geom_vline(xintercept = 0, linetype = 'dashed') +
                   theme(legend.position = 'bottom') +
                   geom_text_repel(data = subset(dat, SIESTA == 1),
                                   aes(label = paste(`Gene.names`)), size = 3, box.padding = unit(0.35, 'lines'),
                                   point.padding = unit(0.3, 'lines')))
  saveWidget(gp,'2D_SIESTA.html')
}
