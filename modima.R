MODIMAstat = function(exposure, mediator, response){
  EM = bcdcor(exposure, mediator)
  MRE = pdcor(mediator, response, exposure)
  MODIMAstat <- EM*MRE
  attr(MODIMAstat, "names") <- "MODIMA stat"
  return(MODIMAstat)
}

permuteDist = function(d){
  d=as.matrix(d)
  n=nrow(d)
  p=sample(n)
  as.dist(d[p,p])
}

modima = function(exposure, mediator, response, nrep=999){
  method <- "WARNING:Specify the number of replicates nrep>0 to perform the test"
  if (!is.null(nrep)){
    nrep <- floor(nrep)
    if (nrep <1)
      nrep <- 0
    if (nrep > 0)
      method <- paste("MODIMA: a Method for Multivariate Omnibus Distance Mediation Analysis")
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
  p.value <- ifelse( nrep>0, ((1+sum(S>MODIMAstat(exposure, mediator, response)))/(nrep+1)), 1)
  bcdcorEM <- bcdcor(exposure, mediator)
  bcdcorER <- bcdcor(exposure, response)
  bcdcorMR <- bcdcor(mediator, response)
  pdcorMRE <- pdcor(mediator, response, exposure)
  attr(bcdcorEM, "names") <- paste0("exposure", "\u2013", "mediator", " bcdcor")
  attr(bcdcorER, "names") <- paste0("exposure", "\u2013", "response", " bcdcor")
  attr(bcdcorMR, "names") <- paste0("mediator", "\u2013", "response", " bcdcor")
  attr(pdcorMRE, "names") <- paste0("mediator", "\u2013", "response", "\u2013", "exposure", " pdcor")
  e <- list(method = method,
            data.name = base::paste("exposure = ", deparse(substitute(exposure)),
                                    "\n       mediator = ", deparse(substitute(mediator)),
                                    "\n       response = ", deparse(substitute(response)),
                                    "\nnumber of permutations + 1:", nrep+1,
                                    "\n sample estimates are \n \t-bias-corrected distance correlation (energy::bcdcor) of indicated pairs and \n \t-partial distance correlation (energy::pdcor) of indicated triple"
                                    #deparse(substitute(exposure)), "and", deparse(substitute(mediator)), "removing", deparse(substitute(response))
                                    ),
            statistic = MODIMAstat(exposure, mediator, response), 
            p.value = p.value,
            estimates = c(bcdcorEM, bcdcorER, bcdcorMR, pdcorMRE)
            #other options for getAnywhere(htest): parameter, alternative, null.value, conf.int
            )
  class(e) <- "htest"
  return(e)
  }

subset.distance = function(d, sub) as.dist(as.matrix(d)[sub, sub])
