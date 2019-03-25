siesta_results <-
  function(data = normdata,
           startPars = c("Tm" = 50, "Pl" = 0, "b" = 0.05),
           vehicle=c(),
           treatment=c(),
           cores = n_cores) {

    cl <- create_cluster(cores = cores)
    cluster_library(cl, "CETSA")
    cluster_library(cl, "ggplot2")
    cluster_library(cl, "gridExtra")
    cluster_library(cl, "broom")
    
    cluster_copy(cl, startPars)
    cluster_copy(cl, fitPeptide)
    cluster_copy(cl, vehicle)
    cluster_copy(cl, treatment)
    set_default_cluster(cl)
    
    data %>%
      partition(id) %>%
      do({
        pepdata <- .
        pid <- unique(pepdata$id)
        models <- try(fitPeptide(pepdata, startPars), silent = T)
        result <- data_frame()
        if (class(models) != 'try-error') {
          models %>%
            group_by(Sample) %>%
            do({
              m = .$model[[1]]
              res <- try(tidy(m), silent = TRUE)
              if (class(res) == 'try-error')
                res <- data.frame()
              else{
                res[, 'sigma'] = .$sigma
                res[, 'rSquared'] = .$rSquared
              }
              res
            }) %>%
            ungroup() %>%
            mutate(id = pid) -> result
        }
        result
      }) -> fitted
    message("analyzed")
    fitted %>% collect()
    
  }
