siesta_clean_proteins <- function(data = data){ 
    data <- data %>%
      filter(is.na(`Only identified by site`)) %>%
      filter(is.na(`Potential contaminant`)) %>%
      filter(is.na(Reverse)) %>%
      filter(Peptides >= 2)
    
    s_completeness <- data %>%
      group_by(id) %>%
      summarise(n = n()) %>%
      mutate(completeness = n/length(temperatures)) %>%
      group_by(completeness) %>%
      summarise(nn = n()) %>%
      mutate(prop.nn = nn/sum(nn) * 100) %>%
      mutate(txt = paste0(round(prop.nn, 2), "%", '\n', 'n=', nn)) %>%
      mutate(completeness = paste0('found in ', completeness, ' replicate(s)'))

    s_completeness$pos = (cumsum(c(0, s_completeness$nn)) + c(s_completeness$nn / 2, .01))[1:nrow(s_completeness)]

    ggplot(s_completeness, aes(1, nn, fill = completeness)) +
      geom_col(color = 'black', position = position_stack(reverse = T), show.legend = T, alpha = 0.7) +
      geom_text_repel(aes(x = 1.4, y = pos, label = txt), nudge_x = 0.3, segment.size = 0.7, show.legend = F) +
      coord_polar('y') +
      theme_void() +
      theme(legend.position = 'bottom') +
      scale_fill_manual(values = c("#F8766D", "#999999", "#4287f5", "#60f542", "#150b3b", "#cf0c0c", "#1ad6b1", "#d6d01a"))
    ggsave("histogram_completeness.pdf")
           
    s_replicates <- data %>%
        group_by(Sample) %>%
        summarise(n = n()) %>%
        mutate(n.proteins = n/length(temperatures)) %>%
        mutate(prop.n.proteins = n.proteins/length(unique(data$id)))

     ggplot(s_replicates) +
        geom_bar(aes(x = Sample, y = prop.n.proteins), position = "dodge", stat = "identity") +
        ylim(0, 1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave("histogram_replicates.pdf")
    
     ggplot(s_replicates) +
        geom_bar(aes(x = Sample, y = n.proteins), position = "dodge", stat = "identity") +
         theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave("relative_histogram_replicates.pdf")
    
    return(data)
  }
