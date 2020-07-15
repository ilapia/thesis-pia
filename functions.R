


###{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library("rstan")
library(matrixStats)
options(mc.cores = parallel::detectCores())

library("raster")
library(Rcpp)
###


### load files

library("rstan")
library(matrixStats)
options(mc.cores = parallel::detectCores())

library("raster")
library(Rcpp)
library(readxl)

sourceCpp("Matern32Covariance.cpp")
sourceCpp("expCovariance.cpp")
  
load("finalTestTrainData")
load("finalDataPred")
###

### functions
###{r raster plot}
load("rastergrids")

# Visualize the study area
raster.plot = function(var, title, points=F, e=ee, r = rr, boxplot=F, m=1,n=1, space=cbind(whitefish.raster$E_etrs89,whitefish.raster$N_etrs89), color=rev(terrain.colors(255))) {
  par(mfrow=c(m,n))
  z <- rasterize(space, r, var, fun=mean)
  plot(z, xlim=cbind(min(whitefish.raster$E_etrs89),max(whitefish.raster$E_etrs89)),
       ylim=cbind(min(whitefish.raster$N_etrs89),max(whitefish.raster$N_etrs89)), main=title, col=color)
  # Plot the locations of sampling sites
  if(points) points(whitefish.dat$E_etrs89,whitefish.dat$N_etrs89, cex=0.7)
  if(boxplot) boxplot(var, main=titles[[i]])
}
###

###{r model fit}
model=function(x,s,y,V, xpr=xpred,
               model.file, dat = whitefish_dat,
               wp=400, it=1000, ch=3, param=c("f", "s2_matern", "l","r") , 
               save=T,model.name=NA, pred.name=NA){
  
  dat <- list(Dx = ncol(x),
              N = nrow(x),
              Ds = ncol(s),
              x = x,
              s = s,
              y = y,
              V = V)
  fit = stan(file =  model.file, data = dat, warmup=wp, iter = it, chains = ch, 
	       #            #    init=list( list(f=as.vector(rep(1,nrow(x))), l=2, s2_matern=1, r=0.2)) ,
		   control = list(adapt_delta = 0.99), pars=param )
  if(save) saveRDS(fit, model.name)	
  m=posterior.check(fit)
  p=prediction(m, expCovariance, x,xpr, s,spred)
  if(save) saveRDS(p, file=pred.name)
  raster.plot(p$Ef, title = "log density" ) 
  
  return(list(fit=fit, pr=p))
}
###

###{r linear cov}
linearCovariance <- function(x1,x2,sigma2) {
  # K = exponentialCovariance(x1,x2,l,sigma2)
  # 
  # A function that takes matrices x1 and x2 of input locations, lengthscale
  # parameter l and magnitude parameter sigma2 and constructs a covariance
  # matrix between f(x1) and f(x2) corresponding to an exponential covariance
  # function. 
  
  K = as.matrix(x1)%*%diag(sigma2)%*%t(as.matrix(x2))
  return(K)
}

###

###{r prediction}
prediction= function(m, Covariance, x, xpred , s, spred, thin=8) {
  m_thinned = as.matrix(m[seq(1,nrow(m),thin),])

  EfMCMC = matrix(nrow=nrow(xpred),ncol=nrow(m_thinned)) # E(f) in each location, for each MCMC iter(thinned)
  sampf_linMCMC = matrix(nrow=ncol(x),ncol=nrow(m_thinned)) # sample f: to get beta estimation for each covariate , in each iter(thinned)
  VarfMCMC = matrix(nrow=nrow(xpred),ncol=nrow(m_thinned)) # Var(f) in each location, for each MCMC iter(thinned)
  
  for (i1 in 1:nrow(m_thinned)){
    l = as.double(m_thinned[i1,"l"])
    sigma2 = as.double(m_thinned[i1,"s2_matern"])
    s2_lin = as.vector(rep(10,ncol(x)))   # Linear covariance function
    
    y = m_thinned[i1,3:(nrow(x)+2)]
  
      Kt = Covariance(s,s,l,sigma2, 1) + linearCovariance(x,x,s2_lin) + diag(rep(1e-1,nrow(x))) #add jiitter 
    Kpt = Covariance(spred, s, l,sigma2,0)
    Kpt_lin = linearCovariance(xpred,x,s2_lin) # xpred * diag(sigma2) * x^T
    Kpt_all = Kpt + Kpt_lin 
    
    # Find the Cholesky Decomposition of the covariance matrix
    L = t( chol(Kt) ) 
    a = solve(t(L),solve(L,y)) #a =(L.t*L)-1 y   (z)
    
    # Plot the different components
    
    #  The full f(t)
    # The posterior predictive mean
    EfMCMC[,i1] = Kpt_all%*%a
    LKtp = solve(L,t(Kpt_all))
    VarfMCMC[,i1] = matrix(sigma2 + xpred^2%*%s2_lin,nrow(Kpt_all),1) - as.matrix( colSums( LKtp*LKtp ) )
    
    # The posterior of the linear weights
    xpred22 = diag(ncol(x))
    Kpt_lin2 = linearCovariance(xpred22,x,s2_lin) # diag(s2)*X^T
    
    Ef_linMCMC = Kpt_lin2%*%a # mean f tilde
    LKtp = solve(L,t(Kpt_lin2))
    Covf =  xpred22%*%diag(s2_lin) - as.matrix(  t(LKtp)%*%LKtp ) # + diag(rep(1e-6,nrow(xpred22)))# cov f tilde
    sampf_linMCMC[,i1] = Ef_linMCMC + chol(Covf)%*%rnorm(nrow(Covf)) # sample beta coeff E(ftilde) +Lz
    #Varf_linMCMC[,i1] = matrix(xpred22^2%*%s2_lin,nrow(Kpt_lin),1) - as.matrix( colSums( LKtp*LKtp ) )
  }
  
  # linear coefficients(weights): beta estimates
  # Posterior mean and standard deviation of fixed effects
    mean=rowMeans(sampf_linMCMC)
   sd=rowSds(sampf_linMCMC)
  # 95% quantiles for fixed effects
  q=apply(sampf_linMCMC, 1, quantile, probs = c(0.025, 0.975))
  
  summary_beta=data.frame(mean,sd,t(q))
  colnames(summary_beta)=c("Mean", "Stand_Dev", "quant2.5%", "quant97.5%")
  if(ncol(x)==15) rownames(summary_beta)=c( "intercept",paste( "BOTTOMCLS", 0:5, sep="_"), colnames(whitefish.dat.cov)[-c(1,2,4,8)])
  else rownames(summary_beta)=c("intercept", paste( "BOTTOMCLS", 0:5, sep="_"), colnames(whitefish.dat.cov)[-c(1,2,4,8)], paste(colnames(whitefish.dat.cov)[-c(1,2,4,8)],"^2", sep=""))
  
  # larvae density f: mean and variance per each location
  Ef = rowMeans(EfMCMC)
  Varf = rowMeans(VarfMCMC) + rowSds(EfMCMC)^2
return( list("summary_beta"=round(summary_beta,3), "Ef"=Ef, "Varf"=Varf, "betas"=sampf_linMCMC, "EfMC"=EfMCMC, "Kpt"=Kpt, "Kpt_lin"=Kpt_lin))
}
###

###{r posterior.check}
posterior.check=function(fit) {
  m = as.matrix(fit)
  #print(fit.w.q)              # Rhat
  print(summary(fit, pars = c("s2_matern", "l", "r"))  )          # summary
  quietgg(stan_trace(fit, pars = c("s2_matern", "l", "r"))  )       # traceplot
  quietgg(stan_ac(fit,inc_warmup = FALSE, lags = 25, pars = c("s2_matern", "l", "r")))   # autocorrelation between Markov chain samples
  quietgg(stan_hist(fit, pars = c("s2_matern", "l", "r")))   # posterior density 
  # scatter plot of parameters of spatial random effect
  quietgg(stan_scat(fit, pars = c("l", "s2_matern"), color = "black", size = 3))
  
  return(m)
}
###

###{r residual}
ee1 <- extent(s2)
rr1 <- raster(ee1, ncol=length(y2), nrow=length(y2))
residual=function(ptest, xpr, ypr, Vpr ){
  p.outlier=which.max(ptest$Ef)
  par(mfrow=c(3,3))
  for (i  in 8:15 ) {
    plot(xpr[,i], ptest$Ef, xlab=covnames[i], main="Predicted target vs cov", pch=18, col="orange")
    points(xpr[p.outlier,i], ptest$Ef[p.outlier], col=2, pch=18, cex=1.3)
  }
  
  par(mfrow=c(3,3))
  #residual plots
  y00=ifelse(ypr==0,0.1,ypr)
  y.logdens= log(y00/Vpr)
#y.logdens[!is.finite(y.logdens)]=0

  resid= y.logdens-ptest$Ef
  resid.std= (resid-mean(resid))/sd(resid)
   for (i  in 8:15 ) {
       plot(xpr[,i], resid.std, xlab=rownames(ptest$summary_beta)[i], main="Std Residual plot", pch=18, col="light green", ylim=c(-5,5))
       abline(h=mean(resid.std), col="dark green")
       #points(x1[which.max(p.tr1$Ef),i], resid[which.max(p.tr1$Ef)], col=2, cex=1.3)
  }
  
  raster.plot(resid, "residuals pred fun", e = ee1, r = rr1, color = rainbow(length(unique(resid))), space = s2 )
  par(mfrow=c(1,2))
  hist(resid)
  boxplot(resid)
  return(resid)
}
###


###{r validate}
validate=function(post,yp,Vp,xp,sp,x,s){
    r = as.matrix(post, pars="r")
  	m = as.matrix(post)
    p = prediction(m, expCovariance, x, xp , s, sp, thin=1)
  lik.yp=matrix(nrow=length(r), ncol=length(yp)) #row=MC iter, col=pred obs

  for (j in 1:length(r)) {
      lik.yp[j,] = dnbinom(yp, mu = Vp*exp(p$EfMC[,j]), size =r[j])
    }
  logLikest = log(colMeans(lik.yp))
  #	return(logLikest)
  KfoldCV = mean(logLikest, na.rm =T)
  return(KfoldCV)
}
###

###{r kfoldcv}
cv1 = function(data, stanmodel, K=5) {
 rows = sample(nrow(data))
  df = data[rows,]
  x = as.matrix(subset(df, select=-c(y,V,s1,s2)))  # target first place
  y = df$y
  V = df$V
  s =cbind(df$s1, df$s2)
  
  B = ceiling(length(y)/K)  # the size of each fold
  n = length(y)
  logLikest = c();
  
  for (i in 1:K) {
  	ind1 = (B*(i-1)+1):min(n, i*B)    # indexes of test data 
  	ind2 = which(y != y[ind1]);  # indexes of training data 
  	# dataset for the ith fold
  	dataset = list(Dx = ncol(x),
  	               Ds = ncol(s),
                   N = nrow(x[ind2,]),
                   x = x[ind2,],
  	               s = s[ind2,],
                   y = y[ind2],
                   V = V[ind2])

  
  	# posterior samples
  	post = stan(model_code =stanmodel, data = dataset, chains = 1, warmup = 150, thin = 1, iter = 500,
  	                control = list(adapt_delta = 0.99), pars=c("f", "s2_matern", "l", "r"))
  	r = as.matrix(post, pars="r")
  	m = as.matrix(post)
  	# log point-wise posterior predictive density estimates
  	Np = length(ind1)
    xp = x[ind1,]
    sp = s[ind1,]
    yp = y[ind1]
    Vp = V[ind1]
    p = prediction(m, expCovariance, x[ind2,], xp , s[ind2,], sp, thin=1)
    
    lik.yp=matrix(nrow=length(r), ncol=length(ind1)) #row=MC iter, col=pred obs
    for (j in 1:length(r)) {
      lik.yp[j,] = dnbinom(yp, mu = Vp*exp(p$EfMC[,j]), size =r[j])
    }
  	logLikest[ind1] = log(colMeans(lik.yp))
  	#print(logLikest)
  }
  print(which(!is.finite(logLikest)))
  KfoldCV = mean(logLikest, na.rm =T)
  return(KfoldCV)
}
###

### r{response curve}
# covariates ranges
# > apply(x[,8:15],2,"min");apply(x[,8:15],2,"max")
# [1] -1.0360879 -0.9819702 -3.2550771 -1.6389181 -1.1424918 -2.9078146 -1.8407923 -1.9312393
# [1] 2.662831 3.075169 1.708342 3.601653 4.047479 1.927289 2.305125 1.982890
# > apply(xpred[,8:15],2,"min");apply(xpred[,8:15],2,"max")
# [1] -1.1164992 -0.9819702 -3.2550771 -1.6572790 -1.1424918 -3.8179518 -1.9496938 -2.2903797
# [1] 3.547355 4.048025 1.708342 4.499599 4.642510 2.042780 2.688552 2.167914

respc= function(pr, XtrainingCont = x1[,8:15], Xminmax = matrix(c(-1.1164992, -0.9819702, -3.2550771, -1.6572790, -1.1424918, -3.8179518, -1.9496938, -2.2903797,
                                        3.547355, 4.048025, 1.708342, 4.499599, 4.642510, 2.042780, 2.688552, 2.167914), ncol=2)){
  
  covariatenames = rownames(pr$summary_beta)[8:15]
  B1 = as.matrix(pr$betas[8:15,])
  par(mfrow=c(3,3))
  boxplot(t(pr$betas[2:7,]), main="Bottom Class")
  for (i in 1:8){
    xtemp = as.vector(seq(Xminmax[i,1],Xminmax[i,2],length.out = 100))
    ftemp = xtemp%*%t(B1[i,])
    Eftemp = apply(ftemp, 1, "mean")
    plot(xtemp,Eftemp,type="l",main=covariatenames[i])
    lines(xtemp,apply(ftemp, 1, "quantile",probs=0.05),lty="dashed")
    lines(xtemp,apply(ftemp, 1, "quantile",probs=0.95),lty="dashed")
    points(XtrainingCont[,i],rep(0,nrow(XtrainingCont)), pch=4)
  }
}