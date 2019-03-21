siesta_clean_proteins <- function(data = data){
    data <- data %>%
      filter(is.na(`Only identified by site`)) %>%
      filter(is.na(`Potential contaminant`)) %>%
      filter(is.na(Reverse)) %>%
      filter(Peptides >= 2)
    
    return(data)
  }
