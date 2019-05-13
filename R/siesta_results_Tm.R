siesta_results_Tm <- function (data = normdata,
                               startPars = c("Pl" = 0, "a" = 550, "b" = 10),
                               vehicle = c(),
                               treatment = c(),
                               cores = n_cores,
                                temperatures = temperatures) 
{
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
  fitted_Tm <- data %>% partition(id) %>% do({
    pepdata <- .
    pid <- unique(pepdata$id)
    models <- try(fitPeptide(pepdata, startPars), silent = T)
    result <- data_frame()
    if (class(models)[1] != "try-error") {
      result <- models %>% group_by(Sample) %>% do({
        meltPStr  <- paste("(1 - Pl) * 1 / (1+exp((x-Tm)/b/x)) + Pl", "-0.5")
        meltPExpr <- parse(text = meltPStr)
        
        m = .$model[[1]]
        res <- try(tidy(m), silent = TRUE)
        if (class(res)[1] == "try-error"){
          res <- data.frame()
        }
        else {
          res[, "sigma"] = .$sigma
          res[, "rSquared"] = .$rSquared
          
          ur <- try(uniroot(f = function(fExpr, Tm, Pl, b ,x){eval(fExpr)},
                        fExpr = meltPExpr, Tm = coef(m)[['Tm']], b = coef(m)[['b']], Pl = coef(m)[['Pl']],
                        interval = c(min(temperatures, max(temperatures)), tol = 0.0001),
                            silent = T)
          if(class(ur) != "try-error"){
            res[, "Tm"] = ur$root
          }
          else{
            res[, "Tm"] = NA
          }
          
        }
        res
      }) %>% ungroup() %>% mutate(id = pid)
    }
    result
  })
  message("analyzed")
  fitted_Tm %>% collect()
}
