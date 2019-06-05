#used to install Jun Chen MedTest::MedOmniTest
# if(!require(MedTest))
   install.packages("/Users/bashir/Box Sync/mediation/Zhang Chen MedTest MedOmniTest/MedTest", repos = NULL, type="source")

#used to install Joshua Sampson MultiMed::medTest
# if(!require(MultiMed))
#   BiocManager::install("MultiMed", version = "3.8")

# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("MultiMed", version = "3.8")

library(permute); packageVersion("permute")
library(lattice); packageVersion("lattice")
library(matrixStats); packageVersion("matrixStats")
library(vegan); packageVersion("vegan")
library(MedTest); packageVersion("MedTest") #Jun Chen packageVersion("MedTest")'1.0'
library(MultiMed); packageVersion("MultiMed") #Sampson cluster version 1.12.0
library(energy); packageVersion("energy") #energy packageVersion on cluster 1.7.2
#library(purrr); packageVersion("purrr") #purrr cluster version 0.2.4
library(iterators); packageVersion("iterators")
library(foreach); packageVersion("foreach")
library(doParallel); packageVersion("doParallel") #doParallel cluster version 1.0.11
doParallel::registerDoParallel(parallel::makeCluster(10))

#c("permute", "lattice", "matrixStats", "vegan", "MedTest", "MultiMed", "energy", "iterators", "foreach", "doParallel") 
#doParallel::registerDoParallel(cores=4)
#cores <- detectCores(logical = FALSE)
#cl <- makeCluster(cores)
#registerDoParallel(cl, cores=cores)

#library(gpuR); packageVersion("gpuR") 
#library("RevoUtilsMath"); packageVersion("RevoUtilsMath")
#getMKLthreads()

#library(snow); packageVersion("snow")
#library(doSNOW); packageVersion("doSNOW")
#cluster <- makeCluster(3)
#cluster <- makeCluster(as.numeric(Sys.getenv('LSB_MAX_NUM_PROCESSORS')))
#registerDoSNOW(cluster)

#Simulate standard normal
simulate_1mediator = function(n, alpha, beta, gamma){
  sigma2E <- 1
  sigma2M <- 1 #- alpha^2*sigma2E
  sigma2R <- 1 #- (alpha*beta+gamma)^2*sigma2E
  E = rnorm(n, 0, sd=sqrt(sigma2E))
  M = alpha*E + rnorm(n, 0, sd=sqrt(sigma2M))
  R = beta*M + gamma*E + rnorm(n, 0 , sd=sqrt(sigma2R))
  result <- list(E=E, M=M, R=R)
  return(result)
}

#semi-partial distance cov
spdcov.test = function (x, y, z, R){
  if (!(class(x) == "dist")) 
    x <- dist(x)
  if (!(class(y) == "dist")) 
    y <- dist(y)
  if (!(class(z) == "dist")) 
    z <- dist(z)
  Dx <- as.matrix(x)
  Dy <- as.matrix(y)
  Dz <- as.matrix(z)
  n <- nrow(Dx)
  Pxz <- energy:::projection(Dx, Dz)  ##we use ::: since these objects are not exported, otherwise we get this error: Error in { : task 1 failed - "'projection' is not an exported object from 'namespace:energy'"
  Py <- energy:::U_center(Dy)
  
  teststat <- n * energy:::U_product(Pxz, Py)
  den <- sqrt(energy:::U_product(Pxz, Pxz) * energy:::U_product(Py, Py))
  if (den > 0) {
    estimate <- teststat/(n * den)
  }
  else estimate <- 0
  bootfn <- function(Pxz, i, Py) {
    energy:::U_product(Pxz[i, i], Py)
  }
  reps <- replicate(R, expr = {
    i <- sample(1:n)
    bootfn(Pxz, i, Py = Py)
  })
  replicates <- n * reps
  pval <- (1 + sum(replicates >= teststat))/(1 + R)
  dataname <- paste("replicates ", R, sep = "")
  names(estimate) <- "spdcor"
  names(teststat) <- "n V^*"
  e <- list(call = match.call(), method = paste("semi-pdcov test", 
                                                sep = ""), statistic = teststat, estimate = estimate, 
            p.value = pval, n = n, replicates = replicates, data.name = dataname)
  class(e) <- "htest"
  return(e)
}

#MODIMA (Alekseyenko)
##MODIMA stat
mediationStat = function(dE, dM, dR){
  EM <- energy::bcdcor(dE, dM)
  MRE <- energy::pdcor(dM, dR, dE)
  res <- EM*MRE
  return(res)
}
#MODIMA p-vals
energy_p = function(E, M, R, nrep=999){
  dE <- stats::dist(E)
  dM <- stats::dist(M)
  dR <- stats::dist(R)
  p1 = energy::dcor.test(dE, dM, R = nrep)$p.value
  p2 = energy::dcor.test(dM, dR, R = nrep)$p.value
  p3 = spdcov.test(dM, dR, dE, R = nrep)$p.value
  p4 = spdcov.test(dR, dM, dE, R = nrep)$p.value
  p5 = energy::pdcov.test(dM, dR, dE, R = nrep)$p.value
  c(mediationStat(dE, dM, dR), p1, p2, p3, p4, p5, max(p1, p2, p3, p4, p5))
}

#MultiMed::medTest (Sampson)
sampson_p = function(E, M, R){
  medSim <- MultiMed::medTest(E=E, M=M, Y=R, nperm=1000)
  sampson.p <- medSim[,"p"]
}

#MedTest::MedOmniTest (Chen)
chen_p = function(E, M, R){
  chen_med <- list(euc=stats::dist(M)
                   #BC=vegdist(otu.tab, method="bray"),
                   #JAC=as.matrix(vegdist(otu.tab, 'jaccard', binary=TRUE)
                   )
  chen_out <- MedTest::MedOmniTest(x=E, y=R, m.list=chen_med, z = NULL, nperm = 999)#999
  chen.p <- chen_out$permP
}

params = expand.grid(replicate=1:1, #1:1000
                     n = c(20, 50, 100, 150, 200), #c(20, 50, 100, 150)
                     alpha = c(0, 0.25, 0.5, 0.75, 1.0), #c(0, 0.25, 0.5, 0.75, 1.0),
                     beta = c(0, 0.25, 0.5, 0.75, 1.0), #c(0, 0.25, 0.5, 0.75, 1.0), 
                     gamma = c(0, 0.1, 0.25, 0.5) #c(0, 0.2)
                     )  

#system.time(
set.seed(843.762)
res_matrix <- foreach::foreach(i=1:nrow(params), 
                               .combine=rbind, 
                               .packages = c("permute", "lattice", "matrixStats", "vegan", "MedTest", "MultiMed", "energy")
                               ) %dopar%
  {
    df <- simulate_1mediator(n = params[i,]$n,
                                   alpha=params[i,]$alpha,
                                   beta=params[i,]$beta,
                                   gamma=params[i,]$gamma)
    sampson.p <- sampson_p(E=df$E, M=df$M, R=df$R)
    chen.p <- chen_p(E=df$E, M=df$M, R=df$R)
    energy.p <- energy_p(E=df$E, M=df$M, R=df$R)
    result <- c(replicate=params$replicate[i],
                   n=params$n[i],
                   alpha=params$alpha[i],
                   beta=params$beta[i],
                   gamma=params$gamma[i],
                   dataE=df[1],
                   dataM=df[2],
                   dataR=df[3],
                   sampson_multimed_p=sampson.p,
                   chen_medtest_p=chen.p,
                   alekseyenko_modima_p=energy.p
                   )
    return(result)
  }
#)

saveRDS(res_matrix, "simulate_MODIMA_1mediator_1000.rds")

# stopImplicitCluster() #we don't need this as clusters were created automatically by registerdoparallel
# stopCluster()
