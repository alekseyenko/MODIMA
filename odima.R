ODIMAstat = function(exposure, mediator, response){
  EM = bcdcor(exposure, mediator)
  MRE = pdcor(mediator, response, exposure)
  ODIMAstat <- EM*MRE
  attr(ODIMAstat, "names") <- "ODIMA stat"
  return(ODIMAstat)
}

permuteDist = function(d){
  d=as.matrix(d)
  n=nrow(d)
  p=sample(n)
  as.dist(d[p,p])
}

odima = function(exposure, mediator, response, nrep=999){
  method <- "WARNING:Specify the number of replicates nrep>0 to perform the test"
  if (!is.null(nrep)){
    nrep <- floor(nrep)
    if (nrep <1)
      nrep <- 0
    if (nrep > 0)
      method <- cat(paste("ODIMA: A Method for Omnibus Distance Mediation Analysis \n sample estimates are bias-corrected distance correlation (bcdcor) of indicated pairs and \n partial distance correlation (pdcor) of",
                      deparse(substitute(exposure)), "and", 
                      deparse(substitute(mediator)), "removing",
                      deparse(substitute(response))))
  }
  else {
    nrep <- 0
  }
  if(bcdcor(exposure, mediator)< pdcor(mediator, response, exposure)){
    S <- replicate(nrep,expr = bcdcor(permuteDist(exposure),mediator)*pdcor(mediator, response, exposure))
  }
  else{
    S <- replicate(nrep, expr = bcdcor(exposure,mediator)*pdcor(mediator, permuteDist(response), exposure))
  }
  p.value <- ifelse( nrep>0, ((1+sum(S>ODIMAstat(exposure, mediator, response)))/(nrep+1)), 1)
  bcdcorEM <- bcdcor(exposure, mediator)
  bcdcorER <- bcdcor(exposure, response)
  bcdcorMR <- bcdcor(mediator, response)
  pdcorMRE <- pdcor(mediator, response, exposure)
  attr(bcdcorEM, "names") <- "Exposure-Mediator bcdcor"
  attr(bcdcorER, "names") <- "Exposure-Response bcdcor"
  attr(bcdcorMR, "names") <- "Mediator-Response bcdcor"
  attr(pdcorMRE, "names") <- "Mediator-Response-Exposure pdcor"
  
  e <- list(method = method,
            data.name = base::paste("number of permutations + 1:", nrep+1, "\n"),
            statistic = ODIMAstat(exposure, mediator, response), 
            p.value = p.value,
            estimates = c(bcdcorEM, bcdcorER, bcdcorMR, pdcorMRE)
            #other options for getAnywhere(htest): parameter, alternative, null.value, conf.int
            )
  class(e) <- "htest"
  return(e)
  }

subset.distance = function(d, sub) as.dist(as.matrix(d)[sub, sub])
