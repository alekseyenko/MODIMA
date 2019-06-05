#used to install Jun Chen MedTest::MedOmniTest
# if(!require(MedTest))
#   install.packages("/Users/bashir/Box Sync/multi med/Jun Chen MedTest/MedTest", repos = NULL, type="source")

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
library(dirmult); packageVersion("dirmult")
library(HMP); packageVersion("HMP")
library(phyloseq); packageVersion("phyloseq")
library(ape); packageVersion("ape")
library(GUniFrac); packageVersion("GUniFrac")
library(pkgmaker); packageVersion("pkgmaker")
library(registry); packageVersion("registry")
library(rngtools); packageVersion("rngtools")
library(doRNG); packageVersion("doRNG")
doParallel::registerDoParallel(parallel::makeCluster(12))

data(saliva)
data(tongue)
gammaS=dirmult(saliva)$gamma
gammaT=dirmult(tonsils)$gamma

simulate_mixture_microbiome = function(p = 0.5, nreads = 10000){
  HMP::Dirichlet.multinomial(nreads*p, gammaS) +
  HMP::Dirichlet.multinomial(nreads*(1-p), gammaT)
  }
# #1 mediator case, simulate std norm
# simulate_1mediator = function(n, alpha, beta, gamma){
#   sigma2E <- 1
#   sigma2M <- 1 #- alpha^2*sigma2E
#   sigma2R <- 1 #- (alpha*beta+gamma)^2*sigma2E
#   E = rnorm(n, 0, sd=sqrt(sigma2E))
#   M = alpha*E + rnorm(n, 0, sd=sqrt(sigma2M))
#   R = beta*M + gamma*E + rnorm(n, 0 , sd=sqrt(sigma2R))
#   result <- list(E=E, M=M, R=R)
#   return(result)
# }
#Multi med case, simulate mixture
simulate_microbiome_EMR = function(n, alpha, beta, gamma){
  E <- rnorm(n = n, mean = 0, sd = sqrt(1))
  p <- exp(alpha * E)/(1 + exp(alpha * E))                # choose proportion of community type 1 (saliva)
  M <- simulate_mixture_microbiome(p)
  d <- vegan::diversity(M)                                # compute log of Shannon diversity
  R <- rnorm(n = n, mean = beta*d + gamma*E, sd = 0.01)   # rnorm(n = n, mean = beta*d + gamma*E, sd = 0.01) 
  dE <- stats::dist(E)
  dR <- stats::dist(R)
  dM_euc <- stats::dist(M)
  otu <- phyloseq::otu_table(M, taxa_are_rows = F)
  tree_rooted <- ape::rtree(n = 21, rooted = T, tip.label = phyloseq::taxa_names(otu))
  phy <- phyloseq::phyloseq(otu, phyloseq::phy_tree(tree_rooted))
  dM_jsd <- phyloseq::distance(phy, method="jsd")
  dM_bc <- vegan::vegdist(M, 'bray')
  unifracs <- GUniFrac::GUniFrac(M, tree_rooted)$unifracs
  #  dM_UniFrac <- unifracs[, , c('d_UW')]
  #  dM_GUniFrac <- unifracs[, , c('d_0.5')]
  dM_WUniFrac <- unifracs[, , c('d_1')]
  result <- list(E = E, M = M, R = R, #proportion = p, #div = d, 
                 dE = dE ,                             #6
                 dR = dR,                              #7
                 dM_euc = dM_euc,                      #8     
                 dM_jsd = dM_jsd,                      #9      
                 dM_bc = dM_bc,                        #10   
                 dM_wunifrac = dM_WUniFrac,            #11                 
                 tree = tree_rooted)                   #12
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
#MODIMA p-vals JSD
energy_p_jsd = function(dE, dM, dR, nrep=999){
  # dE <- stats::dist(E)
  # #dM <- stats::dist(M)
  # dM <- phyloseq::distance(otu_table(M, taxa_are_rows = F), method="jsd")
  # dR <- stats::dist(R)
  p1 = energy::dcor.test(dE, dM, R = nrep)$p.value
  p2 = energy::dcor.test(dM, dR, R = nrep)$p.value
  p3 = spdcov.test(dM, dR, dE, R = nrep)$p.value
  p4 = spdcov.test(dR, dM, dE, R = nrep)$p.value
  p5 = energy::pdcov.test(dM, dR, dE, R = nrep)$p.value
  c(p1, p2, p3, p4, p5, maxP=max(p1, p2, p3, p4, p5), mediationStat(dE, dM, dR))
}
#MODIMA p-vals BC
energy_p_bc = function(dE, dM, dR, nrep=999){
  p1 = energy::dcor.test(dE, dM, R = nrep)$p.value
  p2 = energy::dcor.test(dM, dR, R = nrep)$p.value
  p3 = spdcov.test(dM, dR, dE, R = nrep)$p.value
  p4 = spdcov.test(dR, dM, dE, R = nrep)$p.value
  p5 = energy::pdcov.test(dM, dR, dE, R = nrep)$p.value
  c(p1, p2, p3, p4, p5, maxP=max(p1, p2, p3, p4, p5), mediationStat(dE, dM, dR))
}
#MODIMA p-vals Weighted UniFrac
energy_p_wuni = function(dE, dM, dR, nrep=999){
  p1 = energy::dcor.test(dE, dM, R = nrep)$p.value
  p2 = energy::dcor.test(dM, dR, R = nrep)$p.value
  p3 = spdcov.test(dM, dR, dE, R = nrep)$p.value
  p4 = spdcov.test(dR, dM, dE, R = nrep)$p.value
  p5 = energy::pdcov.test(dM, dR, dE, R = nrep)$p.value
  c(p1, p2, p3, p4, p5, maxP=max(p1, p2, p3, p4, p5), mediationStat(dE, dM, dR))
}
# energy_p_wuni(dE=z$dE, dM=z$dM_wunifrac, dR=z$dR, nrep=999)
# energy_p_wuni(dE=z$dE, dM=as.dist(z$dM_wunifrac), dR=z$dR, nrep=999)
# #MultiMed::medTest (Sampson)
# sampson_p = function(E,M,R){
#   medSim <- MultiMed::medTest(E=E, M=M, Y=R, nperm=1000) #1000
#   sampson.p <- medSim[,"p"]
# }
# #MedTest::MedOmniTest (Chen)
# chen_p = function(E, medlist, R){
#   chen_out <- MedTest::MedOmniTest(x=E, y=R, m.list=medlist, z = NULL, nperm = 999)#999
#   #chen.p <- chen_out$permP
#   return(chen_out)
# }
# params = expand.grid(replicate=1:1, #1:1000
#                      n = c(10), #c(20, 50, 100, 150, 200),
#                      alpha = c(0.25), #c(0, 0.25, 0.5, 0.75, 1.0),
#                      beta = c(0.25), 
#                      gamma = c(0.1) #c(0, 0.1, 0.25, 0.5)
#                      )
params = expand.grid(replicate=1:1000, #1:1000
                     n = c(150), #c(20, 50, 100, 150, 200),
                     alpha = c(0, 0.25, 0.5, 0.75, 1.0), #c(0, 0.25, 0.5, 0.75, 1.0),
                     beta = c(0, 0.25, 0.5, 0.75, 1.0), 
                     gamma = c(0, 0.1, 0.25, 0.5) #c(0, 0.1, 0.25, 0.5)
                     )
# #testing it out
# z=simulate_microbiome_EMR(n=10, alpha=0.2, beta=0.2, gamma=0.2)
# sampson.p=sampson_p(E = z$E, M = z$M, R = z$R)
# energy.p.jsd=energy_p_jsd(E = z$E, M = z$M, R = z$R)
# energy.p.jac=energy_p_jsd(E = z$E, M = z$M, R = z$R)
# chen.p=chen_p(E = z$E, M = z$M, R = z$R)

# df<-simulate_microbiome_EMR(#n = 100, alpha = 0, beta = .25, gamma = .5
#                                   n = params$n[1], alpha = params$alpha[1], beta = params$beta[1], gamma = params$gamma[1])
# chen_p <- MedTest::MedOmniTest(x=df$E,
#                                y=df$R,
#                                m.list=list(BC = df$dM_bc,
#                                            WUniFrac = df$dM_wunifrac,
#                                            JSD = df$dM_jsd),
#                                z = NULL, nperm = 999)
# sampson.p=sampson_p(E = test$E, M = test$M, R = test$R)
# energy.p=energy_p(E = test$E, M = test$M, R = test$R)
# chen.p=chen_p(E = test$E, M = test$M, R = test$R)
# dM_binary <- as.matrix(vegan::vegdist(z$M, method='jaccard', binary=T))
# dM_binary
# dM_binary_logical <- as.matrix(vegan::vegdist(0+(z$M<median(z$M)), method='jaccard', binary=T))
# dM_binary_logical
# dM <- as.matrix(vegan::vegdist(z$M, method='jaccard'))
# dM
# dM_logical <- as.matrix(vegan::vegdist(0+(z$M<median(z$M)), method='jaccard'))
# dM_logical
# dM_jsd <- phyloseq::distance(otu_table(z$M, taxa_are_rows = F), method="jsd")
# dM_jsd
# dM_bray <- phyloseq::distance(otu_table(z$M, taxa_are_rows = F), method="bray")
# dM_bray
# dM_bray2 <- as.matrix(vegan::vegdist(z$M, method='bray'))
# dM_bray2
# z=simulate_microbiome_EMR(n=10, alpha=0.25, beta=0.25, gamma=0.1)
# df = as.data.frame(apply(as.data.frame(mixture_1[,-c(6:8)]),2, unlist))
# rownames(df) = rownames(mixture_1)
# apply(df, 2, class)


#system.time(
#set.seed(843.762, kind = "L'Ecuyer-CMRG")
registerDoRNG(843.762)
res_matrix <- foreach::foreach(i=1:nrow(params),
                               .combine=rbind, #rbind #used for 1mediator case
                               .packages = c("permute", "lattice", "matrixStats", "vegan", "MedTest", "MultiMed", "energy", "HMP", "phyloseq")
                               ) %dopar% 
  {
    df <- simulate_microbiome_EMR(n = params[i,]$n,
                                  alpha=params[i,]$alpha,
                                  beta=params[i,]$beta,
                                  gamma=params[i,]$gamma)
    #sampson.p <- sampson_p(E=df$E, M=df$M, R=df$R)
    chen_p <- MedTest::MedOmniTest(x=df$E,
                                   y=df$R,
                                   m.list=list(BC = df$dM_bc,
                                               WUniFrac = df$dM_wunifrac,
                                               JSD = df$dM_jsd),
                                   z = NULL, nperm = 999)
    energy.p.jsd <- energy_p_jsd(dE=df$dE, dM=df$dM_jsd, dR=df$dR)
    energy.p.bc <- energy_p_bc(dE=df$dE, dM=df$dM_bc, dR=df$dR)
    energy.p.wuni <- energy_p_wuni(dE=df$dE, dM=as.dist(df$dM_bc), dR=df$dR)
    result <- c(replicate = params$replicate[i],
                n = params$n[i],
                alpha = params$alpha[i],
                beta = params$beta[i],
                gamma = params$gamma[i],
                df[c(1:3, 10)],
                chen_medtest_margP = chen_p$margPs,
                chen_medtest_permP = chen_p$permP,
                modima_jsd = energy.p.jsd,
                modima_bc = energy.p.bc,
                modima_WUni = energy.p.wuni
                )
    return(result)
  }
saveRDS(res_matrix, "simulate_MODIMA_mixture_1000_rng_cl12_n150.rds")
sessionInfo()
# df = as.data.frame(apply(as.data.frame(res_matrix[,-c(6:17)]), 2, unlist))
# rownames(df) = rownames(res_matrix)
# apply(df, 2, class)
# str(df)
# apply(df[,-c(1:15)], 2, table, useNA="always")
