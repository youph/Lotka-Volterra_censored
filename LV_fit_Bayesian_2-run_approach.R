# Bayesian inference of Lotka-Volterra system's parameters from time series data of species populations, using MCMC sampling from 
# unnormalized posterior distribution of parameters.
# The 2-run approach:
# "Preliminary" run: starting from a specified param.init point in the parameter space, generate a few relatively short MCMC chains to explore the 
# parameter space a bit. Find the best point found by those chains, and run an optimizer from that point. This gives an improved "best point"
# "Main" run: starting from the improved best point found in the preliminary run, generate a bunch of long MCMC chains, to explore the parameter space around that point.
# Merge and thin the chains from the main run; the resulting chain is the representatation of joint distribution of the L-V parameters learned from data.

# N.B. Before running this for the 1st time:
#   - set the working directory: Go to Session -> Set Working Directory -> To Source File Location

# Check that all required packages are installed. If not, install the missing packages.
list.of.packages <- c("coda", "deSolve", "Matrix", "mcmc", "adaptMCMC", "runjags", "Rmpi")
missing.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)>0) install.packages(missing.packages)

rm(list=ls())

library(coda)
library(deSolve)
library(Matrix)
#-----------------------------------------------------------------------------------------------------------------------
source('log_uposterior_surface.R')
#-----------------------------------------------------------------------------------------------------------------------
heaviside <- function (x, a = 0) 
{
  result = (sign(x - a) + 1)/2
  result[which(result==0.5)] = 1
  result
}
#-----------------------------------------------------------------------------------------------------------------------
plot.dynamics <- function(pop.data, plot.type="l", line.type="solid", line.width=1, pch=19, xlab="", ylab="", newplot=c(TRUE,FALSE))  #plot population dynamics of species
{
  species.names <- names(pop.data[,-1])  #the 1st column is time, so skip it
  max.pop <- max(pop.data[,-1],na.rm=TRUE)
  tmin=min(pop.data$t); tmax=max(pop.data$t);
  for (i in 1:(ncol(pop.data)-1))
  {
    if (newplot==TRUE && i==1) plot(pop.data$t,pop.data[,species.names[[i]]], type=plot.type, lty=line.type, lwd=line.width, pch=pch, col=i, xlim=c(tmin,tmax), ylim=c(0,max.pop), xlab=xlab, ylab=ylab)
    else points(pop.data$t,pop.data[,species.names[[i]]], type=plot.type, pch=pch, lty=line.type, lwd=line.width, col=i)
  }
}
#-----------------------------------------------------------------------------------------------------------------------
X_density_map <- function(mcmc.chain,N.samples,t.series,conf.level=0.95,units=c("rel.user-hours","abs.user-hours","abs.KB"),
                          max.X=NULL,colors=cm.colors(100,alpha=1),z.scale='log',enhancing.factor=1){
  # Generate X.time.series for N.samples random samples from mcmc.chain, and show them all on one density plot:
  sample.ind <- sample.int(n=nrow(mcmc.chain),size=N.samples,replace=F)
  # plot.dynamics(pop.data=Y.series,plot.type='p',newplot=T) 
  X.series.mcmc <- list()
  for (i in 1:N.samples){
    X.series.mcmc[[i]] <- generate_X_time_series(theta.chain[sample.ind[i],],as.data.frame(t(X0.chain[sample.ind[i],])),t.set=t.series,dt,ode.solver,ode.method)
    #   plot.dynamics(pop.data=X.series.mcmc[[i]],plot.type='l',newplot=F)
  }
  
  Y.series.internal <- Y.series
  Y.series.all.internal <- Y.series.all
  
  if (!is.numeric(max.X)){# Define maximum X value in X.series.mcmc:
    max.X <- 0
    for (specie in 1:N.species){
      max.X <- max(max.X, max(X.series.mcmc[[specie]][,-1],na.rm=TRUE))
    }
  }
  if (data.type=="User-hours") y_label <- "Relative usage"
  if (data.type=="Volumes") y_label <- "Relative traffic volumes"
  
  if (units=="abs.user-hours"){#convert all X's from relative units into absolute units
    coef <- max.Y/sampled.fraction
    for (i in 1:N.samples) X.series.mcmc[[i]][,-1] <- X.series.mcmc[[i]][,-1]*coef    
    max.X <- max.X*coef
    Y.series.internal[,-1] <- Y.series[,-1]*coef
    Y.series.all.internal[,-1] <- Y.series.all[,-1]*coef
    y_label <- "Total user-hours in a week"
  }
  
  if (units=="abs.KB" && data.type=="User-hours"){#convert all X's from relative units into absolute KB in a week
    if (Netflix.introduced) {
      average_bitrates_per_bucket_kBps <- read.csv("./Video_data/average_bitrates_per_bucket+Netflix_kBps.csv")
    } else average_bitrates_per_bucket_kBps <- read.csv("./Video_data/average_bitrates_per_bucket_kBps.csv")
    coef_vec <- max.Y/sampled.fraction*average_bitrates_per_bucket_kBps*3600
    for (i in 1:N.samples) {
      for (row in 1:nrow(X.series.mcmc[[i]])) X.series.mcmc[[i]][row,-1] <- X.series.mcmc[[i]][row,-1]*coef_vec
    }
    #     max.X <- max(max.X*coef_vec)
    for (row in 1:nrow(Y.series)) Y.series.internal[row,-1] <- Y.series[row,-1]*coef_vec
    for (row in 1:nrow(Y.series)) Y.series.all.internal[row,-1] <- Y.series.all[row,-1]*coef_vec
    y_label <- "Total KB in a week"
  }
  
  if (units=="abs.KB" && data.type=="Volumes"){#report absolute KB in a week, as predicted
    coef <- max.Y/sampled.fraction
    for (i in 1:N.samples) X.series.mcmc[[i]][,-1] <- X.series.mcmc[[i]][,-1]*coef    
    Y.series.internal[,-1] <- Y.series[,-1]*coef
    Y.series.all.internal[,-1] <- Y.series.all[,-1]*coef
    y_label <- "Total KB in a week"
  }
  
  # Create a matrix of X densities:
  X.density <- list()   #list of matrices (rows=times,cols=1024)
  X.median <- list()
  X.confint.lb <- list()   #lower and upper boundaries of confidence intervals of species populations
  X.confint.ub <- list()
  X.sd <- list()
  for (specie in 1:N.species){
    X.density[[specie]] <- matrix(0,nrow=length(t.series),ncol=1024)
    X.median[[specie]] <- numeric()
    X.sd[[specie]] <- numeric()
    X.confint.lb[[specie]] <- X.confint.ub[[specie]] <- numeric()
    for (t in t.series){
      ind <- which(X.series.mcmc[[1]]$t == t)   #index corresponding to t, same for all samples
      population.t <- numeric()
      for (i in 1:N.samples){
        population.t[i] <- X.series.mcmc[[i]][ind,specie+1]   #concatenate all species populations at time t, corresponding to all samples from the mcmc chain
      }
      if (sum(complete.cases(population.t))>=2) X.density[[specie]][ind,] <- density(population.t, kernel="gaussian", n=1024, from=0, to=max(max.X,max(Y.series.all.internal[,-1],na.rm=TRUE)), na.rm=TRUE)$y
      X.density[[specie]][ind,] <- X.density[[specie]][ind,]*heaviside(t,t.switch.on)[specie]
      X.median[[specie]][ind] <- median(population.t,na.rm = TRUE)
      X.sd[[specie]][ind] <- sd(population.t, na.rm = TRUE)
      confint.boundaries <- quantile(population.t, probs=c((1-conf.level)/2,1-(1-conf.level)/2), na.rm=TRUE)
      X.confint.lb[[specie]][ind] <- confint.boundaries[1]
      X.confint.ub[[specie]][ind] <- confint.boundaries[2]
    }
  }
  
  X.density.all <- matrix(0,nrow=length(t.series),ncol=1024)
  for (specie in 1:N.species) {
    X.density.all <- X.density.all + X.density[[specie]]
  }
  if (z.scale=='lin') z=X.density.all/max(X.density.all)
  if (z.scale=='log') z=log(X.density.all/max(X.density.all)+1)
  image(x=t.series,y=seq(from=0,to=max(max.X,max(Y.series.all[,-1],na.rm=T)),length.out=1024),z=z,zlim=c(0,max(z)/enhancing.factor),col=colors,oldstyle=F,
        xlab="Week",ylab=y_label)
  plot.dynamics(pop.data=Y.series.all.internal,plot.type='p',pch=1,newplot=F) 
  plot.dynamics(pop.data=Y.series.internal,plot.type='p',pch=19,newplot=F) 
  
  for (specie in 1:N.species){
    #Plot median curve for each specie:
    lines(t.series,X.median[[specie]],lwd=3,col=specie)
    #Plot credible intervals around each specie's median:
    lines(t.series,X.confint.lb[[specie]],lwd=3,lty='dotted',col=specie)
    lines(t.series,X.confint.ub[[specie]],lwd=3,lty='dotted',col=specie)
    #     #Plot standard deviation bands around each specie's median:
    #     lines(t.series,X.median[[specie]]-X.sd[[specie]],lwd=3,lty='dotted',col=specie)
    #     lines(t.series,X.median[[specie]]+X.sd[[specie]],lwd=3,lty='dotted',col=specie)
  }
}
#-----------------------------------------------------------------------------------------------------------------------
LV_params_to_theta <- function(r,K,alpha){  #pack r, K, and alpha (competition matrix) into theta (model's parameter vector)
  theta <- numeric()
  ind <- 2*N.species+1
  for (i in 1:N.species){
    theta[i] <- r[i]
    theta[i+N.species] <- K[i]
    for (j in 1:N.species)  
    {
      if(!j==i){  #do not pack the diagonal element of alpha into theta
        theta[ind] <- alpha[i,j]  #skip the first 2*N.species elements of theta - those are for r and K vectors
        ind <- ind+1
      }
    }
  }
  return(theta)
}
#-----------------------------------------------------------------------------------------------------------------------
theta_to_LV_params <- function(theta){  #extract r, K, and alpha (competition matrix) from theta (model's parameter vector)
  r <- theta[1:N.species]
  K <- theta[(N.species+1):(2*N.species)]
  alpha <- matrix(NA,nrow=N.species,ncol=N.species)
  ind <- 2*N.species+1
  for (i in 1:N.species){
    for (j in 1:N.species)  
    {
      if(j==i) alpha[i,j] <- 1  #add diagonal element to alpha
      else {
        alpha[i,j] <- theta[ind]   #skip the first 2*N.species elements of coef - those are r and K vectors
        ind <- ind+1
      }
    }
  }
  if (alpha.positive) alpha <- abs(alpha)  #enforce positive alphas if alpha.positive==TRUE
  if (alpha.symmetric) alpha <- symmetrize(alpha)   #enforce symmetric alpha if alpha.symmetric==TRUE
  return(list('r'=r,'K'=K,'alpha'=alpha))
}
#-----------------------------------------------------------------------------------------------------------------------
pack_parameters <- function(theta,sigma,X0)  #pack theta, sigma, and X0 into a single parameter vector
{
  params <- theta
  params <- c(params,sigma)
  params <- c(params,unlist(X0))
  return(params)
}
unpack_parameters <- function(params)  #unpack theta, sigma, and X0 from a single parameter vector passed as argument
{
  theta <- params[1:(N.species*(N.species+1))]
  sigma <- params[(N.species*(N.species+1)+1)]
  X0 <- params[(N.species*(N.species+1)+2):length(params)]
  return(list('theta'=theta,'sigma'=sigma,'X0'=X0))
}
#-----------------------------------------------------------------------------------------------------------------------
unpack_mcmc_chain <- function(mcmc.chain)  #unpack theta, sigma, and X0 chains from a single parameter mcmc chain passed as argument
{
  theta.chain <- mcmc.chain[,1:(N.species*(N.species+1))]
  r.chain <- theta.chain[,1:N.species]
  K.chain <- theta.chain[,((N.species+1):(2*N.species))]
  sigma.chain <- mcmc.chain[,(N.species*(N.species+1)+1)]
  X0.chain <- mcmc.chain[,(N.species*(N.species+1)+2):ncol(mcmc.chain)]
  return(list('theta.chain'=theta.chain,'r.chain'=r.chain,'K.chain'=K.chain,'sigma.chain'=sigma.chain,'X0.chain'=X0.chain))
}
#-----------------------------------------------------------------------------------------------------------------------
LV_rhs <- function(r,K,alpha,X,t)  #RHS of L-V equations
  # X: current (at time t) populations of all species (vector)
  #   alpha = theta[2*N.species+1:N.species^2+N.species]: off-diagonal elements of the N.species by N.species competition matrix 
{
  S <- alpha %*% (X*heaviside(t,t.switch.on))  #sum(alpha[i,j]*X[j], over all j) - a vector. Note that we include in the RHS only those species for which t>=t.switch.on
  rhs <- r*(1-S/K)*X
  return(rhs)
}
#-----------------------------------------------------------------------------------------------------------------------
LVeqs <- function (time,X,theta) {
  # Extract r, K, alpha from theta:
  tmp <- theta_to_LV_params(theta)
  r <- tmp$r
  K <- tmp$K
  alpha <- tmp$alpha
  
  dX <- LV_rhs(r,K,alpha,X,time)*heaviside(time,t.switch.on)  #only change those species for which t>=t.switch.on
  return(list(dX))
}
#-----------------------------------------------------------------------------------------------------------------------
LV_discrete_rhs.fast <- function(lambda,lambda_over_K,alpha,X,t,dt)  #RHS of discrete-time approximation of L-V equations
  # X: current (at time t) populations of all species (vector)
  #   alpha = theta[2*N+1:N^2+N]: off-diagonal elements of the N by N competition matrix 
{
  S <- alpha %*% t(X*heaviside(t,t.switch.on))  #sum(alpha[i,j]*X[j], over all j) - a vector. Note that we include in the sum only those species for which t>=t.switch.on
  rhs <- lambda/(1 + lambda_over_K*S)*X   
  return(rhs)
}
#-----------------------------------------------------------------------------------------------------------------------
generate_X_time_series <- function(theta,X0,t.set,dt,ode.solver=c('integrate','iterate'),ode.method){#generate time series X(t) with t from t.set, for given theta,X0
  # Solve the L-V equations over the range of times specified by t.set, and return the values of X(t) at t from t.set, assuming that X0=X(min(t.set))
  # Make sure the first element of t.set is equal to the global t.start!
  if(min(t.set)!=t.start) stop("The 1st element of t.set is not equal to t.start! Stopping.")
  
  if (ode.solver=='iterate'){
    # Extract r, K, alpha from theta:
    tmp <- theta_to_LV_params(theta)
    r <- tmp$r
    K <- tmp$K
    alpha <- tmp$alpha
    
    lambda <- exp(r*dt)              #used in LV_discrete_rhs.fast()
    lambda_over_K <- (lambda-1)/K    #used in LV_discrete_rhs.fast()
    
    t.start <- min(t.set)
    t.end <- max(t.set)
    t.fine <- seq(t.start,t.end,by=dt)
    X.series <- list()         #time series X(t) at t from t.set (output of this function)
    X.series[[1]] <- list('t'=t.start,'X'=X0)
    
    X.prev <- X0  #initial value vector
    ind <- 2
    for (t in t.fine[2:length(t.fine)]){
      #     X <- LV_discrete_rhs(r,K,alpha,X=X.prev,dt)
      X <- X.prev + (LV_discrete_rhs.fast(lambda,lambda_over_K,alpha,X=X.prev,t,dt) - X.prev)*heaviside(t,t.switch.on)  #only change those species for which t>=t.switch.on
      if (t %in% t.set) {
        X.series[[ind]] <- list('t'=t,'X'=X)
        ind <- ind+1
      }
      X.prev <- X
    } 
    X.series <- data.frame(matrix(unlist(X.series),nrow=length(t.set), byrow=TRUE))
  }
  
  if (ode.solver=='integrate'){
    require(deSolve)
    X.series <- as.data.frame(ode(func = LVeqs, y = as.numeric(X0), parms = theta, times = t.set, method=ode.method))
  }
  
  colnames(X.series) <- c('t',paste('Species_',seq(1,N.species),sep=''))
  for (row in 1:nrow(X.series)) X.series[row,-1] <- X.series[row,-1]*ifelse(heaviside(X.series[row,'t'],t.switch.on)==1,1,NA) #if the species hasn't been born yet, set its population to NA
  return(X.series)
}
#-----------------------------------------------------------------------------------------------------------------------
temper.Rmpi <- function(logdensity,params.init,scale,niter,Bmin,swap.interval,...){  #temper function to be executed on slaves
  require(deSolve)
  rank=mpi.comm.rank();
  size=mpi.comm.size();
  params.dim <- length(params.init)
  swap=0;
  swaps.attempted=0;  #counter of attempted inter-chain swaps
  swaps.accepted=0;   #counter of accepted inter-chain swaps
  moves.attempted=0;  #counter of attempted within-chain moves
  moves.accepted=0;   #counter of accepted within-chain moves
  
  #Higher ranks run the higher "temperatures" (~smaller fractional powers)
  B=rep(0,size-1);
  for(r in 1:size-1){
    temp=(r-1)/(size-2);
    B[r]=Bmin^temp;
  }
  
  #Create a list for proposal moves
  prop=rep(0,params.dim);
  params=matrix(NA,niter,params.dim)
  log.p=matrix(NA,niter,1)
  
  params[1,] <- params.init
  log.p[1,] <- logdensity(params[1,])
  for(t in 2:niter){
    
    for(c in 1:length(prop)) prop[c]=params[t-1,c]+rnorm(1,0,scale[c]);
    moves.attempted=moves.attempted+1
    
    #Calculate Log-Density at proposed and current position
    logdensity.current=logdensity(params[t-1,])
    logdensity.prop=logdensity(prop);
    
    #Calculate log acceptance probability
    lalpha=B[rank]*(logdensity.prop-logdensity.current)  #note that the current chain's value of B is defined by its rank
    
    if(log(runif(1))<lalpha){
      #Accept proposed move
      params[t,]=prop;
      logdensity.current=logdensity.prop;
      log.p[t,]=logdensity.current
      moves.accepted=moves.accepted+1
    }else{
      #Otherwise do not move
      params[t,]=params[t-1,];
      log.p[t,]=logdensity.current
    } 
    
    #inter-chain swap:
    if(t%%swap.interval ==0){
      for(evenodd in 0:1){
        swap=0;
        logdensity.partner=0;
        if(rank%%2 == evenodd%%2){
          rank.partner=rank + 1;
          #ranks range from 1:size-1. Cannot have a partner rank == size
          if(0<rank.partner && rank.partner<size){
            #On first iteration, evens receive from above odd
            #On second iteration, odds receive from above evens
            logdensity.partner<-mpi.recv.Robj(rank.partner,rank.partner);
            lalpha = (B[rank]-B[rank.partner])*(logdensity.partner-logdensity.current);
            swaps.attempted=swaps.attempted+1;
            if(log(runif(1))<lalpha){
              swap=1;
              swaps.accepted=swaps.accepted+1;
            }
            mpi.send.Robj(swap,dest=rank.partner,tag=rank)
          }
          if(swap==1){
            paramsswap=params[t,];
            mpi.send.Robj(paramsswap,dest=rank.partner,tag=rank)
            params[t,]=mpi.recv.Robj(rank.partner,rank.partner)
            log.p[t,]=logdensity(params[t,])
          }
        }else{
          rank.partner=rank-1;
          #ranks range from 1:size-1. Cannot have a partner rank ==0
          if(0<rank.partner && rank.partner<size){
            #On first iteration, odds send to evens below
            #On second iteration, evens sent to odds below
            mpi.send.Robj(logdensity.current,dest=rank.partner,tag=rank);
            swap=mpi.recv.Robj(rank.partner,rank.partner);
          }
          if(swap==1){
            paramsswap=params[t,];
            params[t,]=mpi.recv.Robj(rank.partner,rank.partner);
            log.p[t,]=logdensity(params[t,])
            mpi.send.Robj(paramsswap,dest=rank.partner,tag=rank);
          }
        }
      }
    }
  }
  return(list('samples'=params,'log.p'=log.p,'move.accept.rate'=moves.accepted/moves.attempted,
              'swap.accept.rate'=swaps.accepted/swaps.attempted))
}
#-----------------------------------------------------------------------------------------------------------------------
log_likelihood <- function(params)  #log-likelihood of observed data Y.series - a function of theta,sigma,X0 (packed into params vector)
{
  params <- unpack_parameters(params)
  theta <- params$theta
  sigma <- params$sigma
  X0 <- as.data.frame(t(params$X0))
  two_sigma_sq <- 2*sigma^2
  X.series.new <- generate_X_time_series(theta,X0,t.set=t.series,dt,ode.solver,ode.method)   #generate underlying population series for given theta,X0
  loglik <- -sum((Y.series[,-1] - X.series.new[,-1])^2,na.rm=TRUE)/two_sigma_sq -    #remember that the 1st columns in X,Y.series contain time
    N.series*N.species/2*log(pi*two_sigma_sq)   #the 2nd term must be evaluated for every new sigma, so can't be pre-computed
  return(loglik)   
}

log_prior <- function(params,type="uniform",lambda_alpha=ifelse((type=="sparse_alpha")||(type=="sparse_alpha&sigma"),1.0,NA))  #log-likelihood of prior over theta,sigma,X0 (packed into params vector)
{
  if (type=="uniform") {
    logprior <- 0
  } else if (type=="sparse_alpha") {  #sparseness-promoting prior for off-diagonal elements of alpha[i,j]: Laplace distribution with zero mean
    alpha <- theta_to_LV_params(unpack_parameters(params)$theta)$alpha
    logprior <- N.species*(N.species-1)*log(lambda_alpha/2) - lambda_alpha*(sum(abs(alpha))-N.species)   #the 2nd term contains the sum of the off-diagonal elements of alpha, hence we sum all elements of alpha and subtract the sum of the diagonal elements, which is equal to N.species
  } else if (type=="sparse_alpha&sigma") {  #sparseness-promoting prior for off-diagonal elements of alpha[i,j], and for sigma: Laplace distribution with zero mean
    lambda_sigma = lambda_alpha/0.5   #defines the width of prior distribution over sigma, in terms of the width of prior distribution over alphas
    params_list <- unpack_parameters(params)
    alpha <- theta_to_LV_params(params_list$theta)$alpha
    sigma <- params_list$sigma
    logprior <- N.species*(N.species-1)*log(lambda_alpha/2) - lambda_alpha*(sum(abs(alpha))-N.species) + log(lambda_sigma/2) - lambda_sigma*abs(sigma)   #the 2nd term contains the sum of the off-diagonal elements of alpha, hence we sum all elements of alpha and subtract the sum of the diagonal elements, which is equal to N.species
  } else {
    stop("Invalid type of prior specified. Stopping.")
  }
  return(logprior)   
}

log_uposterior <- function(params,lb.par=lb.params,ub.par=ub.params)   #log posterior (unnormalized) = log-likelihood + log-prior
{
  if (all(params >= lb.par & params <= ub.par)) return(log_likelihood(params) + log_prior(params,prior.type))
  else return(-Inf)     #this ensures that the mcmc chain stays within the hyperbox defined by lb.params and ub.params
}

neg_log_uposterior <- function(params)   #negative log posterior (unnormalized)
{
  return(-log_uposterior(params,lb.params,ub.params))
}

ludfun <- function(state,lb.params,ub.params) {#Used only with temper()
  d <- length(lb.params)
  ncomp <- N.temperatures
  stopifnot(is.numeric(state))
  stopifnot(length(state) == d + 1)
  icomp <- state[1]    #number of the chain
  stopifnot(icomp == as.integer(icomp))
  stopifnot(1 <= icomp && icomp <= ncomp)
  params <- state[-1]
  return(log_uposterior(params,lb.params,ub.params))
}

symmetrize <- function(alpha){ #symmetrize the alpha matrix
  require(Matrix)
  return(symmpart(alpha))   #computes (alpha+t(alpha))/2
}

save_mcmc.chain <- function(mcmc.chain,file.name){  #save the generated mcmc.chain to a csv file file.name
  mcmc.chain.dataframe <- as.data.frame(mcmc.chain)
  for (ind in 1:(2*N.species)){
    if (ind<=N.species) colnames(mcmc.chain.dataframe)[ind] <- paste("r[",ind,"]",sep="")
    if (N.species<ind && ind<=2*N.species) colnames(mcmc.chain.dataframe)[ind] <- paste("K[",ind-N.species,"]",sep="")
  }
  ind <- 2*N.species+1
  for (i in 1:N.species){
    for (j in 1:N.species)  
    {
      if(j!=i){
        colnames(mcmc.chain.dataframe)[ind] <- paste("alpha[",i,"][",j,"]",sep="")
        ind <- ind+1
      }
    }
  }
  colnames(mcmc.chain.dataframe)[ind] <- "sigma"
  ind <- ind+1
  for (i in 1:N.species){
    colnames(mcmc.chain.dataframe)[ind] <- paste("X0[",i,"]",sep="")
    ind <- ind+1
  }
  
  mcmc.chain.dataframe <- cbind("Names"=c(names(Y.series)[-1],rep(NA,nrow(mcmc.chain.dataframe)-N.species)),mcmc.chain.dataframe)
  
  write.csv(mcmc.chain.dataframe,file=file.name, row.names=FALSE)
}

################################################################################################################################################
# GENERAL PARAMETERS:
################################################################################################################################################
import.data <- TRUE  # If TRUE, import training data, if FALSE, generate synthetic training data
data.type = "Volumes"   #either "User-hours" or "Volumes" - specify depending on what we are modelling
sampled.fraction = 0.02   # fraction of the total national traffic sampled by DPI
alpha.positive = FALSE   #should all alphas be positive? (positive alphas => competition, no predator-prey relations)
alpha.symmetric = FALSE   #should the symmetry of alpha matrix be enforced? (in general it shouldn't)
Netflix.introduced = FALSE  #introduce Netflix into the system?

dt.min <- 0.1   #minimum time step for numeric approximation of L-V eqs
if (!import.data) {
  N.species <- 3
  # Gaussian noise sigma (used only with import.data==FALSE):
  sigma <- 0.1
}

ode.solver = 'integrate'  #either 'integrate' (strongly recommended) or 'iterate'
ode.method = 'lsoda'  #method of numerical integration to use in ode (see help(ode))

prior.type = "sparse_alpha&sigma"   # "uniform" for uniform prior over all parameters,
# "sparse_alpha" for sparseness-promoting prior over alpha[i,j] (off-diagonal), uniform prior over the rest of parameters
# "sparse_alpha&sigma" for sparseness-promoting prior over alpha[i,j] (off-diagonal) and sigma, and uniform prior over the rest of parameters

################################################################################################################################################

if (import.data==TRUE)
{
  #Import observed population data for 3 species:
  # Y.series <- read.csv("US_car_fleet_by_emission_control_systems.csv", header=T)
  # Y.series <- read.csv("Smartphone OS market share USA.csv", header=T); Y.series[,1] = seq(1,nrow(Y.series));
  
  if (data.type=="User-hours") {
    if (Neflix.introduced) Y.series <- read.csv("./Video_data/User-hours weeks 2-20 with Akamai outliers removed + Netflix estimates.csv", header=T);
    if (!Neflix.introduced) Y.series <- read.csv("./Video_data/User-hours weeks 2-20 with Akamai outliers removed.csv", header=T);
  }
  if (data.type=="Volumes") {
    if (Netflix.introduced) Y.series <- read.csv("./Video_data/Traffic volumes weeks 2-20 with Akamai outliers removed + Netflix estimates.csv", header=T);
#     if (!Netflix.introduced) Y.series <- read.csv("./Video_data/Traffic volumes weeks 2-20 with Akamai outliers removed.csv", header=T);
    if (!Netflix.introduced) Y.series <- read.csv("./Video_data/Traffic volumes weeks 2-20.csv", header=T);
  }
  if ("Start.date" %in% names(Y.series)) Y.series <- Y.series[,-which(names(Y.series)=="Start.date")]   #remove "start date" column from the data
  
  #Normalize Y.series to max(Y.series[,-1]), to avoid large values. This is equivalent to measuring X and K in units of max(Y.series[,-1])
  max.Y <- max(Y.series[,-1],na.rm=TRUE)
  Y.series[,-1] <- Y.series[,-1]/max.Y
  N.species <- ncol(Y.series)-1
  colnames(Y.series) <- c('t',names(Y.series[,-1]))
  Y.series.all <- Y.series  #this will be used to compare the fitted model's "forecast" with the actual data for t>t.end.Y
  t.series <- Y.series[,1]
  t.start <- min(t.series)   #t.start is the starting time of observed time series
  t.end <- max(t.series)
  dt <- min(min(diff(t.series)),dt.min)   #internal time step used to solve L-V eqs
  
  t.start.Y <- min(Y.series.all$t)
  t.end.Y <- max(Y.series.all$t)
  ind <- which(!(Y.series$t>=t.start.Y & Y.series$t<=t.end.Y))
  Y.series[ind,2:ncol(Y.series)] <- NA
  
  # Specify "birth" times for all species (time at which the specie is introduced into the system):
  t.switch.on <- list()
  for (specie in 1:N.species){
#     t.switch.on[[specie]] <- t.start   #switch-on at the earliest time for which the data exists (default)
    t.switch.on[[specie]] <- min(Y.series[complete.cases(Y.series[,(1+specie)]),"t"])  #each specie's "birth" time is the time at which its 1st observation was made
  }
  if (Netflix.introduced) t.switch.on[[N.species]] <- min(Y.series[complete.cases(Y.series[,"Netflix"]),"t"])   # Netflix is introduced on week 41
  t.switch.on <- unlist(t.switch.on)
  
  plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,xlab="Week",ylab="Traffic",newplot=T)      #plot observed (noisy) population dynamics
  plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)
  # Show the range of data used for fitting:
  abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
  abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
  legend('right',names(Y.series[,-1]),lty=1, lwd=2, col=seq(1,N.species))
}

if (import.data==FALSE)  # Generate data for N.species, using discrete-time approximation of L-V equations:
{
  if (N.species==3)
  {
    # 3 species:
    r <- c(0.3,0.2,0.2)
    K <- c(2.5,2,2)
    alpha <- t(matrix(c(1,0.4,0.5,0.14,1,0.94,0.24,0.34,1),N.species, N.species))  
    X0 <- data.frame('Species_1'=0.2,'Species_2'=0.1,'Species_3'=0.01)  #vector of initial populations
  }
  if (N.species==4)
  {
    # 4 species:
    r <- c(1,0.72,1.53,1.27)
    K <- c(1.4,1,1,1)
    alpha <- t(matrix(c(1,1.09,1.52,1e-3,1e-3,1,0.44,1.36,2.33,1e-3,1,0.47,1.21,0.51,0.35,1),N.species, N.species))
    X0 <- data.frame('Species_1'=0.1,'Species_2'=0.01,'Species_3'=0.1,'Species_4'=0.01)
  }
  if (alpha.symmetric) alpha <- symmetrize(alpha)
  theta <- LV_params_to_theta(r,K,alpha)
  t.start <- 0   #t.start is the starting time of observed time series
  t.end <- 50   #t.end is the end time of observed time series
  t.series <- seq(t.start,t.end,by=1)    #set of times to generate time series data at. Doesn't have to be uniform
  N.series <- length(t.series)
  dt <- max(dt.min,min(diff(t.series))/10)   #internal time step used to solve L-V eqs
  
  # Specify "birth" times for all species (time at which the specie is introduced into the system):
  t.switch.on <- list()
  for (specie in 1:N.species){
    t.switch.on[[specie]] <- t.start   #switch-on at the earliest time for which the data exists (default)
  }
  t.switch.on <- unlist(t.switch.on)
  
  X.series <- generate_X_time_series(theta,X0,t.set=t.series,dt,ode.solver,ode.method)   #generate underlying population series (unobserved)
  # Generate observed population series Y.series:
  Y.series <- X.series
  #Add Gaussian observational noise:
  N <- dim(Y.series[,2:ncol(Y.series)])
  Y.series[,2:ncol(Y.series)] <- Y.series[,2:ncol(Y.series)] + rnorm(N[1]*N[2],0,sigma)   #this is the observed data (simulated)
  Y.series.all <- Y.series  #this will be used to compare the fitted model's "forecast" with the actual data for t>t.end.Y
  t.start.Y <- 10
  t.end.Y <- 60
  ind <- which(!(X.series$t>=t.start.Y & X.series$t<=t.end.Y))
  Y.series[ind,2:ncol(Y.series)] <- NA
  
  plot.dynamics(pop.data=X.series,plot.type='l',newplot=T)   #plot the underlying (unobserved) population dynamics
  plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,newplot=F)      #plot observed (noisy) population dynamics
  plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)
  # Show the range of data used for fitting:
  abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
  abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
  # Global maximum of log_likelihood, evaluated at true parameter values, just for reference:
  params.true <- pack_parameters(theta,sigma,X0)
  cat("Global maximum of log_likelihood, evaluated at true parameter values:",log_likelihood(params.true),"\n")
  #Now that we have simulated the data using the specified parameters theta, sigma and X0, forget the parameters (they are to be inferred)
  sigma.true <- sigma
  X0.true <- X0
  rm(list=c('theta','sigma','X0'))
}

t.series <- Y.series$t
N.series <- sum(!is.na(Y.series[,2]))
# t.series.extended <- seq(t.start,t.end+50,by=1)   #time range for which to plot the fitted and forecasted series
t.series.extended <- seq(t.start,t.end+50,by=1)

if (import.data==FALSE){
  # lower and upper bounds on parameters:
  #r:
  max.value <- 3
  eps <- 1e-5
  lb.r <- rep(0,N.species)
  ub.r <- rep(max.value,N.species)
  #K:
  lb.K <- rep(eps,N.species)
  ub.K <- rep(max.value,N.species)
  #alpha:
  if (alpha.positive) lb.alpha <- rep(0,N.species^2-N.species) else lb.alpha <- rep(-max.value,N.species^2-N.species)
  ub.alpha <- rep(max.value,N.species^2-N.species)
  #sigma:
  lb.sigma <- rep(0,1)
  ub.sigma <- rep(0.2,1)
  #X0:
  lb.X0 <- rep(0,N.species)
  ub.X0 <- rep(max.value,N.species)
} else {
  max.value <- 3
  eps <- 1e-5
  lb.r <- rep(0,N.species)
  ub.r <- rep(max.value,N.species)
  #K:
  lb.K <- rep(eps,N.species)
  ub.K <- rep(max.value,N.species)
  #alpha:
  if (alpha.positive) lb.alpha <- rep(0,N.species^2-N.species) else lb.alpha <- rep(-max.value,N.species^2-N.species)
  ub.alpha <- rep(max.value,N.species^2-N.species)
  #sigma:
  lb.sigma <- rep(0,1)
  ub.sigma <- rep(max.value,1)
  #X0:
  lb.X0 <- rep(0,N.species)
  ub.X0 <- rep(max.value,N.species)
}
#overall lb and ub:
lb.params <- c(lb.r,lb.K,lb.alpha,lb.sigma,lb.X0)
ub.params <- c(ub.r,ub.K,ub.alpha,ub.sigma,ub.X0)

################################################################################################################################################
# PRELIMINARY RUN:
cat("\nPRELIMINARY RUN:\n")
################################################################################################################################################
set.seed(357)
# set.seed(579)

ode.solver = 'integrate'  #either 'integrate' (strongly recommended) or 'iterate'
ode.method = 'lsoda'  #method of numerical integration to use in ode (see help(ode))

mcmc.method <- "metropolis_adaptive"   #either "metropolis" for Metropolis_Hasting-Green method (library('mcmc')), 
# or "metropolis_adaptive" for adaptive metropolis sampler (library('adaptMCMC')),
# or "tempering" for parallel tempering
burnin.length <- 1   #length of "burn-in" in mcmc chain
chain.length <- 1e6/10   #length of each mcmc chain
thin <- 100/10

N.temperatures <- 6   #(only used with mcmc.method="tempering") number of temperatures (chains) used in parallel tempering
Bmin=0.005    #(used only with tempering.parallel==TRUE) defines the "temperature schedule" for tempering
tempering.parallel <- TRUE  #toggles between 1-CPU and multi-CPU parallel tempering
scale.factor=1e-2   #used in temper.Rmpi. The smaller the scale, the higher shoul dbe the acceptance rate, but the slower is the exploration of the parameter space
swap.interval=2     #(used only with tempering.parallel==TRUE) interval between inter-chain swap attempts

mcmc.parallel <- TRUE   #(only used with mcmc.method=="metropolis_adaptive"): toggles parallel computation of chains on multiple CPU cores
N.chains <- 8    #number of chains to generate when mcmc.parallel==TRUE
N.cpu <- N.chains       #number of CPU cores to use when mcmc.parallel=TRUE
################################################################################################################################################

#vector of scales applied to proposed steps in mcmc chains:
scale.vec <- c(rep(1,N.species),rep(1,N.species),rep(1,N.species*(N.species-1)),0.5,rep(1,N.species))

# # Random starting point in the parameter space, for mcmc chain to start at:
# r.init <- runif(N.species,min=eps,max=max.value)
# K.init <- runif(N.species,min=eps,max=max(Y.series[,-1]))
# alpha.init <- matrix(runif(N.species*N.species,min=0,max=max.value), nrow=N.species, ncol=N.species)
# theta.init <- LV_params_to_theta(r.init,K.init,alpha.init)
# sigma.init <- runif(1,min=eps,max=0.1)

#Fixed starting point in the parameter space, for mcmc chains to start at:
r.init <- rep(eps,N.species)
K.init <- rep(1,N.species)
alpha.init <- matrix(rep(0.1,N.species*N.species), nrow=N.species, ncol=N.species)
theta.init <- LV_params_to_theta(r.init,K.init,alpha.init)
sigma.init <- rep(1e-2,1)

X0.init <- numeric(length=N.species)
for (col in 1:N.species) X0.init[col] <- abs(Y.series.all[,(col+1)][complete.cases(Y.series.all[,(col+1)])][1])
# X0.init <- rep(eps,N.species)
names(X0.init) <- names(Y.series[,-1])

params.init <- pack_parameters(theta.init,sigma.init,X0.init)
params.varnames <- names(params.init)
# params.init <- params.true   #this is for testing purposes only
#Check that params.init is within the constraints:
if (!all(params.init >= lb.params & params.init <= ub.params)) stop("Initial parameters are outside of the parameter hyperbox! Stopping.")

###------------------------------------------------------------------------------------------------------------------------
if (mcmc.method=="metropolis"){#Metropolis-Hasting-Green method of mcmc sampling
  require(mcmc)
  cat("Burning-in with",burnin.length,"samples...")
  mcmc.burnin <- metrop(log_uposterior, initial=params.init, nbatch=burnin.length, blen=1, nspac=1, scale=scale.vec*scale.factor, debug=FALSE, 
                        lb.par=lb.params, ub.par=ub.params)   #burn-in
  mcmc.chain <- mcmc.burnin
  cat("done, time elapsed:",mcmc.burnin$time['elapsed'],"seconds \n")
  cat("Simulating a chain of length",chain.length,"...")
  mcmc.chain <- metrop(mcmc.chain, nbatch=chain.length, lb.par=lb.params, ub.par=ub.params)    #this chain starts where the previous (burn-in) chain has stopped
  cat("done, time elapsed:",mcmc.chain$time['elapsed'],"\n")
  cat("Chain acceptance rate:",mcmc.chain$accept,"\n")
  params.means <- apply(mcmc.chain$batch, 2, mean)  #the grand means (means of batch means)
  # Extract individual paramater chains:
  r.chain <- as.data.frame(matrix(NA,nrow=mcmc.chain$nbatch,ncol=N.species))
  K.chain <- as.data.frame(matrix(NA,nrow=mcmc.chain$nbatch,ncol=N.species))
  alpha.chain <- as.data.frame(matrix(NA,nrow=mcmc.chain$nbatch,ncol=N.species^2))
  sigma.chain <- as.data.frame(matrix(NA,nrow=mcmc.chain$nbatch,ncol=1))
  X0.chain <- as.data.frame(matrix(NA,nrow=mcmc.chain$nbatch,ncol=N.species))
  for (i in 1:mcmc.chain$nbatch){
    row <- unpack_parameters(mcmc.chain$batch[i,])
    theta_i <- row$theta
    tmp <- theta_to_LV_params(theta_i)
    r.chain[i,] <- tmp$r
    K.chain[i,] <- tmp$K
    alpha.chain[i,] <- tmp$alpha
    sigma.chain[i,] <- row$sigma
    X0.chain[i,] <- row$X0
  }
  
  # Estimate mean values of parameters:
  theta.mean <- unpack_parameters(params.means)$theta
  tmp <- theta_to_LV_params(theta.mean)
  r.mean <- tmp$r
  K.mean <- tmp$K
  alpha.mean <- tmp$alpha
  sigma.mean <- unpack_parameters(params.means)$sigma
  X0.mean <- unpack_parameters(params.means)$X0
  
  # Estimate Markov Chain Standard Errors (MCSEs) of params.means:
  params.mcse <- numeric()
  for (i in 1:length(params.means)){
    params.mcse[i] <- sqrt(initseq(mcmc.chain$batch[,i])$var.con/mcmc.chain$nbatch)
  }
  theta.mcse <- unpack_parameters(params.mcse)$theta
  tmp <- theta_to_LV_params(theta.mcse)
  r.mcse <- tmp$r
  K.mcse <- tmp$K
  alpha.mcse <- tmp$alpha
  sigma.mcse <- unpack_parameters(params.mcse)$sigma
  X0.mcse <- unpack_parameters(params.mcse)$X0
  
  cat("\nMCMC estimates for r:\n")
  for (n in 1:N.species){
    cat("95% confidence interval for r[",n,"] is [",r.mean[n] + (-1)*qnorm(0.975)*r.mcse[n],", ",
        r.mean[n] + qnorm(0.975)*r.mcse[n],"]",sep='') 
    cat(" (true value is ",r[n],")\n",sep='')
  }
  cat("\nMCMC estimates for K:\n")
  for (n in 1:N.species){
    cat("95% confidence interval for K[",n,"] is [",K.mean[n] + (-1)*qnorm(0.975)*K.mcse[n],", ",
        K.mean[n] + qnorm(0.975)*K.mcse[n],"]",sep='')
    cat(" (true value is ",K[n],")\n",sep='')
  }
  cat("\nMCMC estimates for sigma:\n")
  cat("95% confidence interval for sigma is [",sigma.mean + (-1)*qnorm(0.975)*sigma.mcse,", ",
      sigma.mean + qnorm(0.975)*sigma.mcse,"]",sep='')
  cat(" (true value is ",sigma.true,")\n",sep='')
}
###------------------------------------------------------------------------------------------------------------------------

###------------------------------------------------------------------------------------------------------------------------
if (mcmc.method=="metropolis_adaptive"){#Metropolis-Hasting-Green method of mcmc sampling
  require(adaptMCMC)
  require(runjags)
  if (mcmc.parallel) {
    cat("Simulating",N.chains,"chains of length",chain.length,"each...")
    system.time(mcmc.chains <- MCMC.parallel(p=log_uposterior, init=params.init, n=chain.length, n.chain=N.chains, 
                                             n.cpu=N.cpu, scale=scale.vec*scale.factor, adapt=TRUE, acc.rate=0.234, 
                                             lb.par=lb.params, ub.par=ub.params))    
    for (i in 1:length(mcmc.chains)){
      cat("Chain",i,"acceptance rate:",mcmc.chains[[i]]$acceptance.rate,"\n")
      cat("Maximum log-density found by chain",i,":",max(mcmc.chains[[i]]$log.p),"\n")
    }
    #Combine all chains into a single chain, taking into account burn-in and thinning:
    mcmc.chains.truncated <- list()
    for (i in 1:length(mcmc.chains)){
      mcmc.chains.truncated[[i]] <- convert.to.coda(mcmc.chains[[i]])[seq(from=burnin.length,to=chain.length,by=thin),]
    }
    mcmc.chain <- combine.mcmc(mcmc.objects=convert.to.coda(mcmc.chains.truncated))
  }
  else {  #not parallel
    cat("Simulating a chain of length",chain.length,"...")
    mcmc.chains <- MCMC(p=log_uposterior, init=params.init, n=chain.length, scale=scale.vec*scale.factor, adapt=TRUE, acc.rate=0.234, 
                        lb.par=lb.params, ub.par=ub.params)  
    cat("Chain acceptance rate:",mcmc.chains$acceptance.rate,"\n")
    mcmc.chain <- convert.to.coda(mcmc.chains)
    #   Take into account burn-in and thinning:
    mcmc.chain <- mcmc.chain[burnin.length:chain.length,]
  }
}
###------------------------------------------------------------------------------------------------------------------------

###------------------------------------------------------------------------------------------------------------------------
if (mcmc.method=="tempering"){
  if (!tempering.parallel){  #single-CPU implementation of parallel tempering
    require(mcmc)
    # Parallel tempering:
    # k: number of chains
    # p: length of parameter vector
    # For parallel tempering the state of the Markov chain is vector of vectors (x1, . . . , xk), where each x is of length p. 
    # This vector of vectors is represented as a k × p matrix.
    neighbors <- matrix(FALSE, N.temperatures, N.temperatures)
    neighbors[row(neighbors) == col(neighbors) + 1] <- TRUE
    neighbors[row(neighbors) == col(neighbors) - 1] <- TRUE
    params.init <- matrix(params.init,nrow=N.temperatures,ncol=length(params.init),byrow=TRUE)
    cat("Simulating",N.temperatures,"tempered chains of length",chain.length,"each...")
    system.time(mcmc.chains <- temper(obj=ludfun, initial=params.init, neighbors, nbatch=chain.length, blen=1, nspac=1,
                                      scale=scale.vec*scale.factor, parallel=TRUE, 
                                      lb.par=lb.params, ub.par=ub.params))  
    cat("done\nWithin-chain acceptance rates for the chains:",mcmc.chains$acceptx,"\n")
    mcmc.chain <- mcmc(mcmc.chains$batch[,1,])
    #   Take into account burn-in and thinning:
    mcmc.chain <- mcmc.chain[burnin.length:chain.length,]
  } else {#multi-CPU implementation of parallel tempering (1 CPU per tempered chain)
    require(Rmpi)
    #Spawn the slaves:
    mpi.spawn.Rslaves(nslaves=N.temperatures)   #spawn as many slaves as there are tempered chains
    mpi.bcast.Robj2slave(all=TRUE)  #send all objects (data and functions) from master's global environment to slaves
    # Execute temper.Rmpi on the slaves:
    mcmc.chain.list=mpi.remote.exec(temper.Rmpi(logdensity=log_uposterior,params.init,scale.vec*scale.factor,niter=chain.length,Bmin,swap.interval,
                                                lb.params,ub.params))
    mpi.close.Rslaves()
    for (slave in names(mcmc.chain.list)){
      cat("Chain from",slave,": within-chain acceptance rate =",mcmc.chain.list[[slave]]$move.accept.rate,
          ", swap acceptance rate =",mcmc.chain.list[[slave]]$swap.accept.rate,"\n")
    }
    require(adaptMCMC)
    mcmc.chain=convert.to.coda(mcmc.chain.list$slave1)   #this is the chain sampled from the target (actual) density (with temperature=1)
    #   Take into account burn-in and thinning:
    mcmc.chain <- mcmc.chain[burnin.length:chain.length,]
  }
}
###------------------------------------------------------------------------------------------------------------------------

if (mcmc.method!="tempering" && mcmc.parallel) gelman.diag(convert.to.coda(mcmc.chains.truncated))  #Gelman and Rubin diagnostic of inter-chain convergence. Potential scale reduction factors 
# should all be less than 1.1, over the meaningful range of the chains
#Extract stuff from mcmc.chain:
params.means <- apply(mcmc.chain, 2, mean)  #means
# Extract individual paramater chains:
tmp <- unpack_mcmc_chain(mcmc.chain)
theta.chain <- tmp$theta.chain
r.chain <- tmp$r.chain
K.chain <- tmp$K.chain
sigma.chain <- tmp$sigma.chain
X0.chain <- tmp$X0.chain

# Some visual diagnostics of the chains:
plot(mcmc(r.chain), trace=FALSE, main="r")
plot(mcmc(K.chain), trace=FALSE, main="K")
# plot(mcmc(alpha.chain))
plot(mcmc(sigma.chain), main="sigma")
# plot(ts(cbind(r.chain,K.chain), names=c("r1","r2","r3","K1","K2","K3")),xlab="iteration", main="MCMC chains for r and K")
# if (mcmc.method!="tempering") plot(mcmc(mcmc.chains$log.p), main="log(P)")  #plot the values of log-posterior in the combined chain
# if (mcmc.method=="tempering") plot(mcmc(mcmc.chain.list$slave1$log.p), main="log(P)")   #plot the values of log-posterior in the lowest temperature chain

ode.method="lsoda"   #for faster integration
# Estimate mean values of parameters from the chain:
theta.mean <- unpack_parameters(params.means)$theta
tmp <- theta_to_LV_params(theta.mean)
r.mean <- tmp$r
K.mean <- tmp$K
alpha.mean <- tmp$alpha
if (alpha.symmetric) alpha.mean <- symmetrize(alpha.mean)
sigma.mean <- unpack_parameters(params.means)$sigma
X0.mean <- unpack_parameters(params.means)$X0

# Plot the solution corresponding to MCMC estimates of parameter means:
X.series.mean <- generate_X_time_series(theta.mean,as.data.frame(t(X0.mean)),t.set=t.series.extended,dt,ode.solver,ode.method)   #generate population series for the inferred parameters
par(mfrow=c(1,1))
if (import.data==FALSE) {
  plot.dynamics(pop.data=X.series,plot.type='l',newplot=T)   #plot the actual underlying (unobserved) population dynamics
} else plot.dynamics(pop.data=X.series.mean,plot.type='l',line.type="dashed",line.width=3,newplot=T) 
plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,newplot=F)      #plot observed (noisy) population dynamics
plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)
abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
plot.dynamics(pop.data=X.series.mean,plot.type='l',line.type="dashed",line.width=3,newplot=F)  #plot the inferred population dynamics
title(main="PRELIMINARY RUN: L-V curves corresponding to MCMC estimates of parameter means")

# Find the best maximum of log_uposterior found by all chains, and plot the solution corresponding to this maximum:
if (mcmc.method!="tempering"){# for "tempering", the values of the objective function are not being output :(
  obj.max <- -Inf
  if (mcmc.parallel){
    for (i in 1:length(mcmc.chains)){
      chain.max <- max(mcmc.chains[[i]]$log.p)
      ind.max <- min(which(mcmc.chains[[i]]$log.p == chain.max))
      if (chain.max > obj.max){
        params.best <- mcmc.chains[[i]]$samples[ind.max,]
        obj.max <- chain.max
      }
    }
  } else {
    chain.max <- max(mcmc.chains$log.p)
    ind.max <- min(which(mcmc.chains$log.p == chain.max))
    if (chain.max > obj.max){
      params.best <- mcmc.chains$samples[ind.max,]
      obj.max <- chain.max
    }
  }
}

if (mcmc.method=="tempering" & mcmc.parallel==TRUE){
  obj.max <- -Inf
  chain.max <- max(mcmc.chain.list$slave1$log.p)
  ind.max <- min(which(mcmc.chain.list$slave1$log.p == chain.max))
  if (chain.max > obj.max){
    params.best <- mcmc.chain.list$slave1$samples[ind.max,]
    obj.max <- chain.max
  }
}

solution <- params.best
theta.mle <- unpack_parameters(solution)$theta
tmp <- theta_to_LV_params(theta.mle)
r.mle <- tmp$r
K.mle <- tmp$K
alpha.mle <- tmp$alpha
if (alpha.symmetric) alpha.mle <- symmetrize(alpha.mle)
sigma.mle <- unpack_parameters(solution)$sigma
X0.mle <- unpack_parameters(solution)$X0

X.series.mle <- generate_X_time_series(theta.mle,as.data.frame(t(X0.mle)),t.set=t.series.extended,dt,ode.solver,ode.method)   #generate population series for the inferred parameters
if (import.data==FALSE) {
  plot.dynamics(pop.data=X.series,plot.type='l',newplot=T)   #plot the actual underlying (unobserved) population dynamics
  plot.dynamics(pop.data=X.series.mle,plot.type='l',line.type="dashed",line.width=3,newplot=F)  #plot the inferred population dynamics
} else plot.dynamics(pop.data=X.series.mle,plot.type='l',line.type="dashed",line.width=3,newplot=T)  #plot the inferred population dynamics
plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,newplot=F)      #plot observed (noisy) population dynamics
plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)
abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
# legend('left',names(Y.series[,-1]),lty='dashed', lwd=2, col=seq(1,N.species))
title(main="PRELIMINARY RUN: L-V curves corresponding to MLE parameters found by MCMC chains")

# Refine the best solution found by MCMC, using an optimizer:
#optim optimizer:
max.iterations <- 1e4
rel.tol = sqrt(.Machine$double.eps)
# optim.algorithm <- "SANN"
# optim.algorithm <- "BFGS"
# optim.algorithm <- "L-BFGS-B"
optim.algorithm <- "Nelder-Mead"
if (optim.algorithm == "L-BFGS-B"){
  params.mle.optim <- optim(par=params.best,fn=neg_log_uposterior,method=optim.algorithm,lower=lb.params,upper=ub.params,
                            control=list(maxit=max.iterations,reltol=rel.tol,trace=FALSE,REPORT=1))
} else {
  params.mle.optim <- optim(par=params.best,fn=neg_log_uposterior,method=optim.algorithm,
                            control=list(maxit=max.iterations,reltol=rel.tol,trace=FALSE,REPORT=1))
}
solution <- params.mle.optim$par
theta.mle.optim <- unpack_parameters(solution)$theta
tmp <- theta_to_LV_params(theta.mle.optim)
r.mle.optim <- tmp$r
K.mle.optim <- tmp$K
alpha.mle.optim <- tmp$alpha
if (alpha.symmetric) alpha.mle.optim <- symmetrize(alpha.mle.optim)
sigma.mle.optim <- unpack_parameters(solution)$sigma
X0.mle.optim <- unpack_parameters(solution)$X0

X.series.optim <- generate_X_time_series(theta.mle.optim,as.data.frame(t(X0.mle.optim)),t.set=t.series.extended,dt,ode.solver,ode.method)   #generate population series for the inferred parameters
if (import.data==FALSE) {
  plot.dynamics(pop.data=X.series,plot.type='l',newplot=T)   #plot the actual underlying (unobserved) population dynamics
  plot.dynamics(pop.data=X.series.optim,plot.type='l',line.type="dashed",line.width=3,newplot=F)  #plot the inferred population dynamics
} else plot.dynamics(pop.data=X.series.optim,plot.type='l',line.type="dashed",line.width=3,newplot=T)  #plot the inferred population dynamics
plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,newplot=F)      #plot observed (noisy) population dynamics
plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)
abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
# legend('left',names(Y.series[,-1]),lty='dashed', lwd=2, col=seq(1,N.species))
title(main="PRELIMINARY RUN: Best MLE solution found by MCMC, refined by optim")

################################################################################################################################################
# MAIN RUN:
cat("MAIN RUN:\n")
################################################################################################################################################
# set.seed(357)

ode.solver = 'integrate'  #either 'integrate' (strongly recommended) or 'iterate'
ode.method = 'lsoda'  #method of numerical integration to use in ode (see help(ode))

mcmc.method <- "metropolis_adaptive"   #either "metropolis" for Metropolis_Hasting-Green method (library('mcmc')), 
# or "metropolis_adaptive" for adaptive metropolis sampler (library('adaptMCMC')),
# or "tempering" for parallel tempering
burnin.length <- 100   #length of "burn-in" in mcmc chain
chain.length <- 1e6/10   #length of each mcmc chain
thin <- 100/10

N.temperatures <- 6   #(only used with mcmc.method="tempering") number of temperatures (chains) used in parallel tempering
Bmin=0.005    #(used only with tempering.parallel==TRUE) defines the "temperature schedule" for tempering
tempering.parallel <- TRUE  #toggles between 1-CPU and multi-CPU parallel tempering
scale.factor=1e-2   #used in temper.Rmpi. The smaller the scale, the higher should be the acceptance rate, but the slower is the exploration of the parameter space
swap.interval=2     #(used only with tempering.parallel==TRUE) interval between inter-chain swap attempts

mcmc.parallel <- TRUE   #(only used with mcmc.method=="metropolis_adaptive"): toggles parallel computation of chains on multiple CPU cores
N.chains <- 8    #number of chains to generate when mcmc.parallel==TRUE
N.cpu <- N.chains       #number of CPU cores to use when mcmc.parallel=TRUE
################################################################################################################################################

# Start MCMC chains from the MLE parameter estimate, refined by optim in the preliminary run:
params.init <- params.mle.optim$par  
names(params.init) <- params.varnames
#Check that params.init is within the constraints:
if (!all(params.init >= lb.params & params.init <= ub.params)) stop("Initial parameters are outside of the parameter hyperbox! Stopping.")

###------------------------------------------------------------------------------------------------------------------------
if (mcmc.method=="metropolis"){#Metropolis-Hasting-Green method of mcmc sampling
  require(mcmc)
  cat("Burning-in with",burnin.length,"samples...")
  mcmc.burnin <- metrop(log_uposterior, initial=params.init, nbatch=burnin.length, blen=1, nspac=1, scale=scale.vec*scale.factor, debug=FALSE, 
                        lb.par=lb.params, ub.par=ub.params)   #burn-in
  mcmc.chain <- mcmc.burnin
  cat("done, time elapsed:",mcmc.burnin$time['elapsed'],"seconds \n")
  cat("Simulating a chain of length",chain.length,"...")
  mcmc.chain <- metrop(mcmc.chain, nbatch=chain.length, lb.par=lb.params, ub.par=ub.params)    #this chain starts where the previous (burn-in) chain has stopped
  cat("done, time elapsed:",mcmc.chain$time['elapsed'],"\n")
  cat("Chain acceptance rate:",mcmc.chain$accept,"\n")
  params.means <- apply(mcmc.chain$batch, 2, mean)  #the grand means (means of batch means)
  # Extract individual paramater chains:
  r.chain <- as.data.frame(matrix(NA,nrow=mcmc.chain$nbatch,ncol=N.species))
  K.chain <- as.data.frame(matrix(NA,nrow=mcmc.chain$nbatch,ncol=N.species))
  alpha.chain <- as.data.frame(matrix(NA,nrow=mcmc.chain$nbatch,ncol=N.species^2))
  sigma.chain <- as.data.frame(matrix(NA,nrow=mcmc.chain$nbatch,ncol=1))
  X0.chain <- as.data.frame(matrix(NA,nrow=mcmc.chain$nbatch,ncol=N.species))
  for (i in 1:mcmc.chain$nbatch){
    row <- unpack_parameters(mcmc.chain$batch[i,])
    theta_i <- row$theta
    tmp <- theta_to_LV_params(theta_i)
    r.chain[i,] <- tmp$r
    K.chain[i,] <- tmp$K
    alpha.chain[i,] <- tmp$alpha
    sigma.chain[i,] <- row$sigma
    X0.chain[i,] <- row$X0
  }
  
  # Estimate mean values of parameters:
  theta.mean <- unpack_parameters(params.means)$theta
  tmp <- theta_to_LV_params(theta.mean)
  r.mean <- tmp$r
  K.mean <- tmp$K
  alpha.mean <- tmp$alpha
  sigma.mean <- unpack_parameters(params.means)$sigma
  X0.mean <- unpack_parameters(params.means)$X0
  
  # Estimate Markov Chain Standard Errors (MCSEs) of params.means:
  params.mcse <- numeric()
  for (i in 1:length(params.means)){
    params.mcse[i] <- sqrt(initseq(mcmc.chain$batch[,i])$var.con/mcmc.chain$nbatch)
  }
  theta.mcse <- unpack_parameters(params.mcse)$theta
  tmp <- theta_to_LV_params(theta.mcse)
  r.mcse <- tmp$r
  K.mcse <- tmp$K
  alpha.mcse <- tmp$alpha
  sigma.mcse <- unpack_parameters(params.mcse)$sigma
  X0.mcse <- unpack_parameters(params.mcse)$X0
  
  cat("\nMCMC estimates for r:\n")
  for (n in 1:N.species){
    cat("95% confidence interval for r[",n,"] is [",r.mean[n] + (-1)*qnorm(0.975)*r.mcse[n],", ",
        r.mean[n] + qnorm(0.975)*r.mcse[n],"]",sep='') 
    cat(" (true value is ",r[n],")\n",sep='')
  }
  cat("\nMCMC estimates for K:\n")
  for (n in 1:N.species){
    cat("95% confidence interval for K[",n,"] is [",K.mean[n] + (-1)*qnorm(0.975)*K.mcse[n],", ",
        K.mean[n] + qnorm(0.975)*K.mcse[n],"]",sep='')
    cat(" (true value is ",K[n],")\n",sep='')
  }
  cat("\nMCMC estimates for sigma:\n")
  cat("95% confidence interval for sigma is [",sigma.mean + (-1)*qnorm(0.975)*sigma.mcse,", ",
      sigma.mean + qnorm(0.975)*sigma.mcse,"]",sep='')
  cat(" (true value is ",sigma.true,")\n",sep='')
}
###------------------------------------------------------------------------------------------------------------------------

###------------------------------------------------------------------------------------------------------------------------
if (mcmc.method=="metropolis_adaptive"){#Metropolis-Hasting-Green method of mcmc sampling
  require(adaptMCMC)
  require(runjags)
  if (mcmc.parallel) {
    cat("Simulating",N.chains,"chains of length",chain.length,"each...")
    system.time(mcmc.chains <- MCMC.parallel(p=log_uposterior, init=params.init, n=chain.length, n.chain=N.chains, 
                                             n.cpu=N.cpu, scale=scale.vec*scale.factor, adapt=TRUE, acc.rate=0.234, 
                                             lb.par=lb.params, ub.par=ub.params))    
    for (i in 1:length(mcmc.chains)){
      cat("Chain",i,"acceptance rate:",mcmc.chains[[i]]$acceptance.rate,"\n")
      cat("Maximum log-density found by chain",i,":",max(mcmc.chains[[i]]$log.p),"\n")
    }
    #Combine all chains into a single chain, taking into account burn-in and thinning:
    mcmc.chains.truncated <- list()
    for (i in 1:length(mcmc.chains)){
      mcmc.chains.truncated[[i]] <- convert.to.coda(mcmc.chains[[i]])[seq(from=burnin.length,to=chain.length,by=thin),]
    }
    mcmc.chain <- combine.mcmc(mcmc.objects=convert.to.coda(mcmc.chains.truncated))
  }
  else {  #not parallel
    cat("Simulating a chain of length",chain.length,"...")
    mcmc.chains <- MCMC(p=log_uposterior, init=params.init, n=chain.length, scale=scale.vec*scale.factor, adapt=TRUE, acc.rate=0.234, 
                        lb.par=lb.params, ub.par=ub.params)  
    cat("Chain acceptance rate:",mcmc.chains$acceptance.rate,"\n")
    mcmc.chain <- convert.to.coda(mcmc.chains)
    #   Take into account burn-in and thinning:
    mcmc.chain <- mcmc.chain[burnin.length:chain.length,]
  }
}
###------------------------------------------------------------------------------------------------------------------------

###------------------------------------------------------------------------------------------------------------------------
if (mcmc.method=="tempering"){
  if (!tempering.parallel){  #single-CPU implementation of parallel tempering
    require(mcmc)
    # Parallel tempering:
    # k: number of chains
    # p: length of parameter vector
    # For parallel tempering the state of the Markov chain is vector of vectors (x1, . . . , xk), where each x is of length p. 
    # This vector of vectors is represented as a k × p matrix.
    neighbors <- matrix(FALSE, N.temperatures, N.temperatures)
    neighbors[row(neighbors) == col(neighbors) + 1] <- TRUE
    neighbors[row(neighbors) == col(neighbors) - 1] <- TRUE
    params.init <- matrix(params.init,nrow=N.temperatures,ncol=length(params.init),byrow=TRUE)
    cat("Simulating",N.temperatures,"tempered chains of length",chain.length,"each...")
    system.time(mcmc.chains <- temper(obj=ludfun, initial=params.init, neighbors, nbatch=chain.length, blen=1, nspac=1,
                                      scale=scale.vec*scale.factor, parallel=TRUE, 
                                      lb.par=lb.params, ub.par=ub.params))  
    cat("done\nWithin-chain acceptance rates for the chains:",mcmc.chains$acceptx,"\n")
    mcmc.chain <- mcmc(mcmc.chains$batch[,1,])
    #   Take into account burn-in and thinning:
    mcmc.chain <- mcmc.chain[burnin.length:chain.length,]
  } else {#multi-CPU implementation of parallel tempering (1 CPU per tempered chain)
    require(Rmpi)
    #Spawn the slaves:
    mpi.spawn.Rslaves(nslaves=N.temperatures)   #spawn as many slaves as there are tempered chains
    mpi.bcast.Robj2slave(all=TRUE)  #send all objects (data and functions) from master's global environment to slaves
    # Execute temper.Rmpi on the slaves:
    mcmc.chain.list=mpi.remote.exec(temper.Rmpi(logdensity=log_uposterior,params.init,scale.vec*scale.factor,niter=chain.length,Bmin,swap.interval,
                                                lb.params,ub.params))
    mpi.close.Rslaves()
    for (slave in names(mcmc.chain.list)){
      cat("Chain from",slave,": within-chain acceptance rate =",mcmc.chain.list[[slave]]$move.accept.rate,
          ", swap acceptance rate =",mcmc.chain.list[[slave]]$swap.accept.rate,"\n")
    }
    require(adaptMCMC)
    mcmc.chain=convert.to.coda(mcmc.chain.list$slave1)   #this is the chain sampled from the target (actual) density (with temperature=1)
    #   Take into account burn-in and thinning:
    mcmc.chain <- mcmc.chain[burnin.length:chain.length,]
  }
}
###------------------------------------------------------------------------------------------------------------------------

if (mcmc.method!="tempering" && mcmc.parallel) gelman.diag(convert.to.coda(mcmc.chains.truncated))  #Gelman and Rubin diagnostic of inter-chain convergence. Potential scale reduction factors 
# should all be less than 1.1, over the meaningful range of the chains
#Extract stuff from mcmc.chain:
params.means <- apply(mcmc.chain, 2, mean)  #means
# Extract individual paramater chains:
tmp <- unpack_mcmc_chain(mcmc.chain)
theta.chain <- tmp$theta.chain
r.chain <- tmp$r.chain
K.chain <- tmp$K.chain
sigma.chain <- tmp$sigma.chain
X0.chain <- tmp$X0.chain

# Some visual diagnostics of the chains:
plot(mcmc(r.chain), trace=FALSE, main="r")
plot(mcmc(K.chain), trace=FALSE, main="K")
# plot(mcmc(alpha.chain))
plot(mcmc(sigma.chain), main="sigma")
# plot(ts(cbind(r.chain,K.chain), names=c("r1","r2","r3","K1","K2","K3")),xlab="iteration", main="MCMC chains for r and K")
# if (mcmc.method!="tempering") plot(mcmc(mcmc.chains$log.p), main="log(P)")  #plot the values of log-uposterior in the combined chain
# if (mcmc.method=="tempering") plot(mcmc(mcmc.chain.list$slave1$log.p), main="log(P)")   #plot the values of log-posterior in the lowest temperature chain

ode.method="lsoda"   #for faster integration in stiff parameter regions
# Estimate mean values of parameters from the chain:
theta.mean <- unpack_parameters(params.means)$theta
tmp <- theta_to_LV_params(theta.mean)
r.mean <- tmp$r
K.mean <- tmp$K
alpha.mean <- tmp$alpha
if (alpha.symmetric) alpha.mean <- symmetrize(alpha.mean)
sigma.mean <- unpack_parameters(params.means)$sigma
X0.mean <- unpack_parameters(params.means)$X0

# Plot the solution corresponding to MCMC estimates of parameter means:
X.series.mean <- generate_X_time_series(theta.mean,as.data.frame(t(X0.mean)),t.set=t.series.extended,dt,ode.solver,ode.method)   #generate population series for the inferred parameters
par(mfrow=c(1,1))
if (import.data==FALSE) {
  plot.dynamics(pop.data=X.series,plot.type='l',newplot=T)   #plot the actual underlying (unobserved) population dynamics
} else plot.dynamics(pop.data=X.series.mean,plot.type='l',line.type="dashed",line.width=3,newplot=T) 
plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,newplot=F)      #plot observed (noisy) population dynamics
plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)
abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
plot.dynamics(pop.data=X.series.mean,plot.type='l',line.type="dashed",line.width=3,newplot=F)  #plot the inferred population dynamics
title(main="MAIN RUN: L-V curves corresponding to MCMC estimates of parameter means")

# Find the best maximum of log_uposterior found by all chains, and plot the solution corresponding to this maximum:
if (mcmc.method!="tempering"){# for "tempering", the values of the objective function are not being output :(
  obj.max <- -Inf
  if (mcmc.parallel){
    for (i in 1:length(mcmc.chains)){
      chain.max <- max(mcmc.chains[[i]]$log.p)
      ind.max <- min(which(mcmc.chains[[i]]$log.p == chain.max))
      if (chain.max > obj.max){
        params.best <- mcmc.chains[[i]]$samples[ind.max,]
        obj.max <- chain.max
      }
    }
  } else {
    chain.max <- max(mcmc.chains$log.p)
    ind.max <- min(which(mcmc.chains$log.p == chain.max))
    if (chain.max > obj.max){
      params.best <- mcmc.chains$samples[ind.max,]
      obj.max <- chain.max
    }
  }
}

if (mcmc.method=="tempering" & mcmc.parallel==TRUE){
  obj.max <- -Inf
  chain.max <- max(mcmc.chain.list$slave1$log.p)
  ind.max <- min(which(mcmc.chain.list$slave1$log.p == chain.max))
  if (chain.max > obj.max){
    params.best <- mcmc.chain.list$slave1$samples[ind.max,]
    obj.max <- chain.max
  }
}

solution <- params.best
theta.mle <- unpack_parameters(solution)$theta
tmp <- theta_to_LV_params(theta.mle)
r.mle <- tmp$r
K.mle <- tmp$K
alpha.mle <- tmp$alpha
if (alpha.symmetric) alpha.mle <- symmetrize(alpha.mle)
sigma.mle <- unpack_parameters(solution)$sigma
X0.mle <- unpack_parameters(solution)$X0

X.series.mle <- generate_X_time_series(theta.mle,as.data.frame(t(X0.mle)),t.set=t.series.extended,dt,ode.solver,ode.method)   #generate population series for the inferred parameters
if (import.data==FALSE) {
  plot.dynamics(pop.data=X.series,plot.type='l',newplot=T)   #plot the actual underlying (unobserved) population dynamics
  plot.dynamics(pop.data=X.series.mle,plot.type='l',line.type="dashed",line.width=3,newplot=F)  #plot the inferred population dynamics
} else plot.dynamics(pop.data=X.series.mle,plot.type='l',line.type="dashed",line.width=3,newplot=T)  #plot the inferred population dynamics
plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,newplot=F)      #plot observed (noisy) population dynamics
plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)
abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
legend('left',names(Y.series[,-1]),lty='dashed', lwd=2, col=seq(1,N.species))
title(main="MAIN RUN: L-V curves corresponding to MLE parameters found by MCMC chains")

# Refine the best solution found by MCMC, using an optimizer:
#optim optimizer:
max.iterations <- 1e4
rel.tol = sqrt(.Machine$double.eps)
# optim.algorithm <- "SANN"
# optim.algorithm <- "BFGS"
# optim.algorithm <- "L-BFGS-B"
optim.algorithm <- "Nelder-Mead"
if (optim.algorithm == "L-BFGS-B"){
  params.mle.optim <- optim(par=params.best,fn=neg_log_uposterior,method=optim.algorithm,lower=lb.params,upper=ub.params,
                            control=list(maxit=max.iterations,reltol=rel.tol,trace=FALSE,REPORT=1))
} else {
  params.mle.optim <- optim(par=params.best,fn=neg_log_uposterior,method=optim.algorithm,
                            control=list(maxit=max.iterations,reltol=rel.tol,trace=FALSE,REPORT=1))
}
solution <- params.mle.optim$par
theta.mle.optim <- unpack_parameters(solution)$theta
tmp <- theta_to_LV_params(theta.mle.optim)
r.mle.optim <- tmp$r
K.mle.optim <- tmp$K
alpha.mle.optim <- tmp$alpha
if (alpha.symmetric) alpha.mle.optim <- symmetrize(alpha.mle.optim)
sigma.mle.optim <- unpack_parameters(solution)$sigma
X0.mle.optim <- unpack_parameters(solution)$X0

X.series.optim <- generate_X_time_series(theta.mle.optim,as.data.frame(t(X0.mle.optim)),t.set=t.series.extended,dt,ode.solver,ode.method)   #generate population series for the inferred parameters
if (import.data==FALSE) {
  plot.dynamics(pop.data=X.series,plot.type='l',newplot=T)   #plot the actual underlying (unobserved) population dynamics
  plot.dynamics(pop.data=X.series.optim,plot.type='l',line.type="dashed",line.width=3,newplot=F)  #plot the inferred population dynamics
} else plot.dynamics(pop.data=X.series.optim,plot.type='l',line.type="dashed",line.width=3,newplot=T)  #plot the inferred population dynamics
plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,newplot=F)      #plot observed (noisy) population dynamics
plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)
abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
# legend('left',names(Y.series[,-1]),lty='dashed', lwd=2, col=seq(1,N.species))
title(main="MAIN RUN: Best MLE solution found by MCMC, refined by optim")

# Generate X.time.series for a bunch of random samples from mcmc.chain, and show them all on one density plot:
if (import.data){
  if (data.type=="User-hours"){
    # Relative user-hours:
    N.samples <- 100
    #     set.seed(1234567)
    X_density_map(mcmc.chain,N.samples,t.series=t.series.extended,conf.level=0.8,units="rel.user-hours",max.X=1.2,colors=cm.colors(100,alpha=1),z.scale='log',enhancing.factor=3)
    abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
    abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
    # legend('left',names(Y.series[,-1]), lty=1, lwd=2, col=seq(1,N.species))
    title("MAIN RUN: Probabilistic forecast of user-hours in a week")
    
    # Absolute KB in a week:
    #     set.seed(1234567)
    X_density_map(mcmc.chain,N.samples,t.series=t.series.extended,conf.level=0.8,units="abs.KB",max.X=1.25e13,colors=cm.colors(100,alpha=1),z.scale='log',enhancing.factor=5)
    abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
    abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
    legend('topleft',names(Y.series[,-1]), lty=1, lwd=2, col=seq(1,N.species))
    title("MAIN RUN: Probabilistic forecast of KB in a week") 
  }
  
  if (data.type=="Volumes"){
    # Relative volumes:
    N.samples <- 1000
    #     set.seed(1234567)
    X_density_map(mcmc.chain,N.samples,t.series=t.series.extended,conf.level=0.8,units="rel.volumes",max.X=1.2,colors=cm.colors(100,alpha=1),z.scale='log',enhancing.factor=3)
    abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
    abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
    # legend('left',names(Y.series[,-1]), lty=1, lwd=2, col=seq(1,N.species))
    title("MAIN RUN: Probabilistic forecast of relative volumes in a week")
    
    # Absolute KB in a week:
    #     set.seed(1234567)
    X_density_map(mcmc.chain,N.samples,t.series=t.series.extended,conf.level=0.8,units="abs.KB",max.X=3.25e13,colors=cm.colors(100,alpha=1),z.scale='log',enhancing.factor=5)
    abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
    abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
    legend('topleft',names(Y.series[,-1]), lty=1, lwd=2, col=seq(1,N.species))
    title("MAIN RUN: Probabilistic forecast of KB in a week") 
  } 
}

if (!import.data){
  N.samples <- 100
  set.seed(1234567)
  X_density_map(mcmc.chain,N.samples,t.series=t.series.extended,conf.level=0.8,units="rel.user-hours",colors=cm.colors(100,alpha=1),z.scale='log',enhancing.factor=3)
  abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
  abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
  # legend('left',names(Y.series[,-1]), lty=1, lwd=2, col=seq(1,N.species))
  title("MAIN RUN: Probabilistic forecast of relative usage in a week")
}

# specie=4  #choose the specie for which to make the r vs K surface slice
# log_uposterior_surface(params=params.mle.optim$par,slice.dimension=c(specie,N.species+specie),
#                        n.gridpoints=100,zlim=c(0,220), x.name=paste('r[',specie,']',sep=''),y.name=paste('K[',specie,']',sep=''))