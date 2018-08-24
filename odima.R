mediationStat = function(exposure, mediator, response){
  EM = bcdcor(exposure, mediator)
  MRE = pdcor(mediator, response, exposure)
  mediationStat <- EM*MRE
  attr(mediationStat, "names") <- "MediationStat"
  return(mediationStat)
}

permuteDist = function(d){
  d=as.matrix(d)
  n=nrow(d)
  p=sample(n)
  as.dist(d[p,p])
}

pdCorMediationTest = function(exposure, mediator, response, nrep=99999){
  if (is.null(nrep))
    nrep <- 0
  if(bcdcor(exposure, mediator)< pdcor(mediator, response, exposure)){
    S <- replicate(nrep,expr = bcdcor(permuteDist(exposure),mediator)*pdcor(mediator, response, exposure))
  }
  else{
    S <- replicate(nrep, expr = bcdcor(exposure,mediator)*pdcor(mediator, permuteDist(response), exposure))
  }
  mediationTest <- (1+sum(S>mediationStat(exposure, mediator, response)))/(nrep+1)
  dataname <- paste("replicates ", nrep, sep = "")
  method <- ifelse(nrep > 0, "pdCor mediation test of independence of mediator and response removing exposure", "WARNING:Specify the number of replicates nrep>0 to perform the test")
  bcdcor_estimates <- c(bcdcor(exposure, mediator)[1], bcdcor(exposure, response), bcdcor(mediator, response)#, pdcor(distjsd_early, dR_early, dE_early)
                        )
    e <- list(method = method, 
            statistic = mediationStat(exposure, mediator, response), 
            estimates = bcdcor_estimates,
            p.value = mediationTest,
            replicates = nrep,
            data.name = dataname)
  class(e) <- "htest"
  return(e)
  }

subset.distance = function(d, sub) as.dist(as.matrix(d)[sub, sub])
