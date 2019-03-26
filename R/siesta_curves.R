siesta_curves <- function(data = normdata,
           startPars = c("Tm" = 50, "Pl" = 0, "b" = 0.05),
           plotCurves = F,
           resultPath = '.',
           vehicle=c(),
           treatment=c(),
           cores = n_cores) {
    cl <- create_cluster(cores = cores)
    cluster_library(cl, "CETSA")
    cluster_library(cl, "ggplot2")
    cluster_library(cl, "gridExtra")
    cluster_library(cl, "broom")
    
    cluster_copy(cl, startPars)
    cluster_copy(cl, plotCurves)
    cluster_copy(cl, resultPath)
    cluster_copy(cl, fitPeptide)
    cluster_copy(cl, vehicle)
    cluster_copy(cl, treatment)
    set_default_cluster(cl)
    if (!dir.exists(file.path(resultPath, "plots"))){dir.create(file.path(resultPath, "plots"), recursive = T)}
    
    data %>%
      partition(id) %>%
      do({
        pepdata <- .
        pid <- unique(pepdata$id)
        models <- try(fitPeptide(pepdata, startPars), silent = T)
        pepdata_m <- data_frame()
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
          
          if (nrow(result) >= 1) {
            temps <- unique(pepdata$Temperature)
            xtemps <- seq(min(temps), max(temps), length.out = 100)
            pepdata_m =  models %>%
              group_by(Sample) %>%
              do(data.frame(
                Sample = .$Sample,
                Temperature = xtemps,
                prValue = predict(.$model[[1]], list(x = xtemps))
              )) %>%
              full_join(pepdata, by = c("Sample", "Temperature")) %>%
              arrange(Sample, Temperature) %>%
              mutate(Type = Sample)
            
            for(i in 1:nrow(pepdata_m)){
              pepdata_m$Type[i] <- substr(pepdata_m$Type[i], 1, (nchar(pepdata_m$Type[i])-1))
            }
            
            pepdata_m$id <- pid
          }
        }
        pepdata_m
      }) -> fitted
    message("fitted")
    fitted %>% collect()
    
  }
