#' Put main.R & supportFun.R into the identical folder.
#' In the same folder, create two subfolders named figure' & 'Rimage', respectively.

rm(list = ls())
if (!("rstudioapi" %in% rownames(installed.packages()))) 
  install.packages("rstudioapi")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("functions.R")

#### Global settings
options(warn = -1, digits = 4)
set.seed(1)
FVE.x = .99
FVE.y = .99
RR = 200 # number of replication
propTrain = .8 # proportion of real dataset left for training
realCase = 1 # 1 for DTI; 2 for BG
simu = F # T for simulation, F for real data
simuCase = 2 # 1 for beta = P_2(s)P_2(t); 2 for beta = P_4(s)P_4(t)
rho = .9 # control the autocorrelation for epsilon
SNR = 5 # signal-to-noise-ratio for simulation

#### Simulated data ####

if (simu == T) {
  N = 300
  tRange = c(0, 1)
  domain.x = seq(from = tRange[1], to = tRange[2], length.out = 101)
  domain.y = seq(from = tRange[1], to = tRange[2], length.out = 100)
  
  # Beta.true
  Lambda = c(100, 10, 1)
  J = length(Lambda) # number of basis functions
  eigenFun.x = orthoBasis(order = 2:(J+1), denseGrid = domain.x, type = 'shiftedLegendre', normalized = T)
  eigenFun.y = orthoBasis(order = 2:(J+1), denseGrid = domain.y, type = 'shiftedLegendre', normalized = T)
  Beta.true = switch(simuCase,
                     as.matrix(eigenFun.x[1, ]) %*%  t(as.matrix(eigenFun.y[1, ])),
                     as.matrix(eigenFun.x[J, ]) %*%  t(as.matrix(eigenFun.y[J, ]))
  )
  
  # covariance for epsilon

  sigma = switch (simuCase,
    Lambda[1]^.5/SNR,
    Lambda[J]^.5/SNR
  )
  Sigma = matrix(NA, nrow = length(domain.y), ncol = length(domain.y))
  for (i in 1:nrow(Sigma)){
    for (j in 1:ncol(Sigma)){
      Sigma[i, j] = (rho^(abs(domain.y[i] - domain.y[j]))) * sigma^2
    }
  }
  
  # Create x and y
  x.all = list()
  y.all = list()
  for (R in 1:RR) {
    eigenScore = MASS::mvrnorm(n = N, mu = numeric(J), Sigma = diag(Lambda))
    x.all[[R]] = eigenScore %*% eigenFun.x
    epsilon = MASS::mvrnorm(n = N, mu = numeric(length(domain.y)), Sigma = Sigma)
    y.all[[R]] = integral(x.all[[R]], Beta.true, domain = domain.x, type = 222) + epsilon
  }
}


#### Real data ####

if (simu == F){
  if (realCase == 1) {
    # cca vs. rcst (DTI in R package classiFunc)
    
    x = classiFunc::DTI$cca
    y = classiFunc::DTI$rcst
    
    domain.x = (1:ncol(x))/ncol(x)
    domain.y = (1:ncol(y))/ncol(y)
    
    N = nrow(x)
    p.max = pUpper.compu(x, domain.x, FVE.x, basis.name = "bspline")
  }
  
  if (realCase == 2){
    # Boy's gait (gait in R package fda)
    
    x = t(fda::gait[,,"Hip Angle"])
    y = t(fda::gait[,,"Knee Angle"])
    
    domain.x = seq(from=0.025, to=0.975, by=0.05)
    domain.y = seq(from=0.025, to=0.975, by=0.05)
    
    x = smooth.curve(x)
    y = smooth.curve(y)
    N = nrow(x)
    p.max = pUpper.compu(x, domain.x, FVE.x, basis.name = "fourier")
  }
}

# Dangerous! clear existing result
res.fAPLS = list(); time.fAPLS = numeric(); reispe.fAPLS = numeric(); reisee.fAPLS = numeric()
res.SigComp = list(); time.SigComp = numeric(); reispe.SigComp = numeric(); reisee.SigComp = numeric()
res.NIPALS = list(); time.NIPALS = numeric(); reispe.NIPALS = numeric(); reisee.NIPALS = numeric()
res.SIMPLS = list(); time.SIMPLS = numeric(); reispe.SIMPLS = numeric(); reisee.SIMPLS = numeric()

# Check current progress
check.fAPLS = length(time.fAPLS)
check.SigComp = length(time.SigComp)
check.NIPALS = length(time.NIPALS)
check.SIMPLS = length(time.SIMPLS)

################ Replica ##################
for (R in 1:RR) {
  
  if (simu == T) {
    
    sampIdx = sample(1:N, round(N * propTrain))
    x.old = x.all[[R]][sampIdx, ]
    y.old = y.all[[R]][sampIdx, ]
    x.new = x.all[[R]][-sampIdx, ]
    y.new = y.all[[R]][-sampIdx, ]
    p.max = pUpper.compu(x.all[[R]], domain.x, FVE.x, basis.name = "bspline")

    # regression
    if (R > check.NIPALS){
      ptm0 = proc.time()[3]
      res.NIPALS[[R]] = NIPALS(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.NIPALS[R] = ptm1 - ptm0
      
      reispe.NIPALS[R] = res.NIPALS[[R]]$reispe
      integ1 = integral((res.NIPALS[[R]]$Beta - Beta.true)^2, domain = domain.y, type = 201)
      integ2 = integral((res.NIPALS[[R]]$Beta)^2, domain = domain.y, type = 201)
      reisee.NIPALS[R] = integral(integ1, domain = domain.x, type = 100)/
        integral(integ2, domain = domain.x, type = 100)
    }
    
    if (R > check.SIMPLS){
      ptm0 = proc.time()[3]
      res.SIMPLS[[R]] = SIMPLS(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.SIMPLS[R] = ptm1 - ptm0
      
      reispe.SIMPLS[R] = res.SIMPLS[[R]]$reispe
      integ1 = integral((res.SIMPLS[[R]]$Beta - Beta.true)^2, domain = domain.y, type = 201)
      integ2 = integral((res.SIMPLS[[R]]$Beta)^2, domain = domain.y, type = 201)
      reisee.SIMPLS[R] = integral(integ1, domain = domain.x, type = 100)/
        integral(integ2, domain = domain.x, type = 100)
    }
  
    if (R > check.fAPLS){
      ptm0 = proc.time()[3]
      res.fAPLS[[R]] = fAPLS(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max, tune = T)
      ptm1 = proc.time()[3]
      time.fAPLS[R] = ptm1 - ptm0
      
      reispe.fAPLS[R] = res.fAPLS[[R]]$reispe
      integ1 = integral((res.fAPLS[[R]]$Beta - Beta.true)^2, domain = domain.y, type = 201)
      integ2 = integral((res.fAPLS[[R]]$Beta)^2, domain = domain.y, type = 201)
      reisee.fAPLS[R] = integral(integ1, domain = domain.x, type = 100)/
        integral(integ2, domain = domain.x, type = 100)
    }
    
    if (R > check.SigComp){
      ptm0 = proc.time()[3]
      res.SigComp[[R]] = SigComp(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.SigComp[R] = ptm1 - ptm0
      
      reispe.SigComp[R] = res.SigComp[[R]]$reispe
      integ1 = integral((res.SigComp[[R]]$Beta - Beta.true)^2, domain = domain.y, type = 201)
      integ2 = integral((res.SigComp[[R]]$Beta)^2, domain = domain.y, type = 201)
      reisee.SigComp[R] = integral(integ1, domain = domain.x, type = 100)/
        integral(integ2, domain = domain.x, type = 100)
    }
  }
  
  if (simu == F) {
    
    sampIdx = sample(1:N, round(N * propTrain))
    x.old = x[sampIdx, ]
    y.old = y[sampIdx, ]
    x.new = x[-sampIdx, ]
    y.new = y[-sampIdx, ]
    
    # regression
    if (R > check.NIPALS){
      ptm0 = proc.time()[3]
      res.NIPALS[[R]] = NIPALS(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.NIPALS[R] = ptm1 - ptm0

      reispe.NIPALS[R] = res.NIPALS[[R]]$reispe
    }
    
    if (R > check.SIMPLS){
      ptm0 = proc.time()[3]
      res.SIMPLS[[R]] = SIMPLS(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.SIMPLS[R] = ptm1 - ptm0
      
      reispe.SIMPLS[R] = res.SIMPLS[[R]]$reispe
    }
    
    if (R > check.fAPLS){
      ptm0 = proc.time()[3]
      res.fAPLS[[R]] = fAPLS(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max, tune = T)
      ptm1 = proc.time()[3]
      time.fAPLS[R] = ptm1 - ptm0
      
      reispe.fAPLS[R] = res.fAPLS[[R]]$reispe
    }
    
    if (R > check.SigComp){
      ptm0 = proc.time()[3]
      res.SigComp[[R]] = SigComp(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.SigComp[R] = ptm1 - ptm0
      
      reispe.SigComp[R] = res.SigComp[[R]]$reispe
    }
  }
  
  if (R == 1){
    file = ifelse(simu,
                  paste0('Rimage/',
                         'Simu', simuCase, '_',
                         J, 'eigen_',
                         RR, 'repeats_',
                         propTrain * 100, 'train_',
                         rho * 1e1L, 'rho_',
                         SNR, 'SNR_',
                         FVE.x * 1e2L, 'FVEx.RData'),
                  paste0('Rimage/',
                         'Real', realCase, '_', 
                         RR, 'repeats_',
                         propTrain * 100, 'train_',
                         FVE.x * 1e2L, 'FVEx.RData')
    )
  }
  if (R %% 40 == 0){
    cat(R, '\n')
  }else 
    cat(R)
  
  # save the R space
  if (R == RR){
    if (simu == T){
      # shrink the R space for simulation
      rm(x.all, y.all)
    }
    save.image(file = file)
  }
}

################ Plots and summaries of results ##################

# ReISEPE matrix and boxplots
reispeLst = list(
  reispe.fAPLS
  ,reispe.SigComp
  ,reispe.NIPALS
  ,reispe.SIMPLS
)
names(reispeLst) = c(
  'fAPLS'
  ,'SigComp'
  ,'NIPALS'
  ,'SIMPLS'
)
reispeMat = creatErrMat(reispeLst)
round(cbind(
  apply(reispeMat, 2, mean),
  apply(reispeMat, 2, sd)
), digits = 2)
source("functions.R")
boxplotErr(RR, reispeMat, type = 'ReISPE')

# ReISEE matrix and boxplots
if (simu == T){
  reiseeLst = list(
    reisee.fAPLS
    ,reisee.SigComp
    ,reisee.NIPALS
    ,reisee.SIMPLS
  )
  names(reiseeLst) = names(reispeLst)
  reiseeMat = creatErrMat(reiseeLst)
  round(cbind(
    apply(reiseeMat, 2, mean),
    apply(reiseeMat, 2, sd)
  ), digits = 2)
  boxplotErr(RR, reiseeMat, type = 'ReISEE')
}

# Running time
round(c(
  sum(time.fAPLS)
  ,sum(time.SigComp)
  ,sum(time.NIPALS)
  ,sum(time.SIMPLS)
), digits = 1)
