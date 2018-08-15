mediationStat = function(exposure, mediator, response){
  EM = bcdcor(exposure, mediator)
  MRE = pdcor(mediator, response, exposure)
  EM*MRE
}

permuteDist = function(d){
  d=as.matrix(d)
  n=nrow(d)
  p=sample(n)
  as.dist(d[p,p])
}

pdCorMediationTest = function(exposure, mediator, response, nrep=99999){
  if(bcdcor(exposure, mediator)< pdcor(mediator, response, exposure)){
    S = replicate(nrep,expr = bcdcor(permuteDist(exposure),mediator)*pdcor(mediator, response, exposure))
  }
  else{
    S = replicate(nrep, expr = bcdcor(exposure,mediator)*pdcor(mediator, permuteDist(response), exposure))
  }
  (1+sum(S>mediationStat(exposure, mediator, response)))/(nrep+1)
}

subset.distance = function(d, sub) as.dist(as.matrix(d)[sub, sub])
