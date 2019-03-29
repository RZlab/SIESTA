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
      mutate(prop.nn = nn/sum(nn) * 100)

    g <- ggplot(s_completeness, aes(x = "", y = nn, fill = as.factor(completeness)))+
        geom_bar(width = 1, stat = "identity") +
        coord_polar("y", start=0) +
        theme_void() +
        geom_text(aes(label = paste0(round(prop.nn, 1), "%")), position = position_stack(vjust = 0.5)) +
        theme(legend.position = "bottom")
    ggsave(g, "histogram_completeness.pdf"
           
    s_replicates <- data %>%
        group_by(Sample) %>%
        summarise(n = n()) %>%
        mutate(n.proteins = n/length(temperatures)) %>%
        mutate(prop.n.proteins = n.proteins/length(unique(data$id)))

     g <- ggplot(s_replicates) +
        geom_bar(aes(x = Sample, y = prop.n.proteins), position = "dodge", stat = "identity") +
        ylim(0, 1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(g, "histogram_replicates.pdf"
    
    return(data)
  }
