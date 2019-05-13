fitPeptide <- function(pepdata, startPars = c("Pl" = 0, "a" = 550, "b" = 10)){
  pepdata %>%
    group_by(Sample) %>%
    do(
      model=fitSigmoid(.[,c("Temperature","Value")],startPars),
      yVec = .$Value
    ) %>%
    filter(class(model)=='nls') %>%
    rowwise() %>%
    summarise(
      Sample=Sample,
      sigma=sigma(model),
      rSquared = rSquared(model, yVec),
      model=list(model)) %>%
    arrange(sigma)
}
