# Parallel implementation of MLE inference of Lotka-Volterra system's parameters from time series data of species populations
rm(list=ls())
library('nloptr')
library(deSolve)
if (!is.loaded("mpi_initialize")) {
  library("Rmpi")
}

.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}

plot.dynamics <- function(pop.data, plot.type="l", line.type="solid", line.width=1, pch=19, newplot=c(TRUE,FALSE))  #plot population dynamics of species
{
  species.names <- names(pop.data[,-1])  #the 1st column is time, so skip it
  max.pop <- max(pop.data[,-1],na.rm=TRUE)
  tmin=min(pop.data$t); tmax=max(pop.data$t);
  for (i in 1:(ncol(pop.data)-1))
  {
    if (newplot==TRUE && i==1) plot(pop.data$t,pop.data[,species.names[[i]]], type=plot.type, lty=line.type, lwd=line.width, pch=pch, col=i, xlim=c(tmin,tmax), ylim=c(0,max.pop), xlab="Time", ylab="Population")
    else points(pop.data$t,pop.data[,species.names[[i]]], type=plot.type, pch=pch, lty=line.type, lwd=line.width, col=i)
  }
}

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
  return(list('r'=r,'K'=K,'alpha'=alpha))
}

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

LV_rhs <- function(r,K,alpha,X)  #RHS of discrete-time approximation of L-V equations
  # X: current (at time t) populations of all species (vector)
  #   alpha = theta[2*N+1:N^2+N]: off-diagonal elements of the N by N competition matrix 
{
  S <- alpha %*% X  #sum(alpha[i,j]*X[j], over all j) - a vector
  rhs <- r*(1-S/K)*X
  return(rhs)
}

LVeqs <- function (time,X,theta) {
  # Extract r, K, alpha from theta:
  tmp <- theta_to_LV_params(theta)
  r <- tmp$r
  K <- tmp$K
  alpha <- tmp$alpha
  
  dX <- LV_rhs(r,K,alpha,X)
  return(list(dX))
}

LV_discrete_rhs <- function(r,K,alpha,X,dt)  #RHS of discrete-time approximation of L-V equations
  # X: current (at time t) populations of all species (vector)
  # theta: vector of L-V model's coefficients, with the following components:
  #   r = theta[1:N]: "birth rates" of all N species 
  #   K = theta[N+1:2*N]: carrying capacities of all species 
  #   alpha = theta[2*N+1:N^2+N]: off-diagonal elements of the N by N competition matrix 
{
  lambda <- exp(r*dt)  #vector
  S <- alpha %*% t(X)  #sum(alpha[i,j]*X[j], over all j) - a vector
  rhs <- lambda*X/(1 + (lambda-1)/K*S)
  return(rhs)
}

LV_discrete_rhs.fast <- function(lambda,lambda_over_K,alpha,X,dt)  #RHS of discrete-time approximation of L-V equations
  # X: current (at time t) populations of all species (vector)
  #   alpha = theta[2*N+1:N^2+N]: off-diagonal elements of the N by N competition matrix 
{
  S <- alpha %*% t(X)  #sum(alpha[i,j]*X[j], over all j) - a vector
  rhs <- lambda/(1 + lambda_over_K*S)*X
  return(rhs)
}

generate_X_time_series <- function(theta,X0,t.set,dt,method=c('integrate','iterate')){#generate time series X(t) with t from t.set, for given theta,X0
  # Solve the L-V equations over the range of times specified by t.set, and return the values of X(t) at t from t.set, assuming that X0=X(min(t.set))
  # Make sure the first element of t.set is equal to the global t.start!
  if(min(t.set)!=t.start) stop("The 1st element of t.set is not equal to t.start! Stopping.")
  
  # Extract r, K, alpha from theta:
  tmp <- theta_to_LV_params(theta)
  r <- tmp$r
  K <- tmp$K
  alpha <- tmp$alpha
  if (method=='iterate'){
    lambda <- exp(r*dt)              #used in LV_discrete_rhs.fast()
    lambda_over_K <- (lambda-1)/K    #used in LV_discrete_rhs.fast()
    
    t.start <- min(t.set)
    t.end <- max(t.set)
    t.fine <- seq(t.start,t.end,by=dt)
    X.series <- list()         #time series X(t) at t from t.set (output of this function)
    X.series[[1]] <- list('t'=t.start,'X'=X0)
    
    X.prev <- X0  #initial value
    ind <- 2
    for (t in t.fine[2:length(t.fine)]){
      #     X <- LV_discrete_rhs(r,K,alpha,X=X.prev,dt)
      X <- LV_discrete_rhs.fast(lambda,lambda_over_K,alpha,X=X.prev,dt)
      if (t %in% t.set) {
        X.series[[ind]] <- list('t'=t,'X'=X)
        ind <- ind+1
      }
      X.prev <- X
    } 
    X.series <- data.frame(matrix(unlist(X.series),nrow=length(t.set), byrow=TRUE))
  }
  
  if (method=='integrate'){
    X.series <- as.data.frame(ode(func = LVeqs, y = as.numeric(X0), parms = theta, times = t.set, method="ode45"))
  }
  
  colnames(X.series) <- c('t',paste('Species_',seq(1,N.species),sep=''))
  return(X.series)
}

neg_log_likelihood <- function(params)  #log-likelihood of observed data Y.series - a function of theta,sigma,X0 (packed into params vector)
{
  params <- unpack_parameters(params)
  theta <- params$theta
  sigma <- params$sigma
  X0 <- as.data.frame(t(params$X0))
  two_sigma_sq <- 2*sigma^2
  X.series.new <- generate_X_time_series(theta,X0,t.set=t.series,dt,ode.solver)   #generate underlying population series for given theta,X0
  loglik <- -sum((Y.series[,-1] - X.series.new[,-1])^2,na.rm=TRUE)/two_sigma_sq -     #remember that the 1st columns in X,Y.series contain time
    N.series*N.species/2*log(pi*two_sigma_sq)
  return(-loglik)   #return -loglik, since optimizer (nlopt or optim) looks for minimum by default
}

# Function the slaves will call to perform optimization of the objective function, starting from a random initial point in the parameter
# space, picked by the master.
# Assumes: params.init (initial point in the parameter space), algorithm (optimization algoritms)
# Also assumes all the required functions have been broadcast to the slave.
slave_optimize <- function() {
  require(nloptr)
  require(deSolve)
  # Note the use of the tag for sent messages: 
  #     1=ready_for_task, 2=done_task, 3=exiting 
  # Note the use of the tag for received messages: 
  #     1=task, 2=done_tasks 
  junk <- 0
  done <- FALSE
  while (!done){
    # Signal the master of being ready to receive a new task 
    mpi.send.Robj(junk,0,1)
    # Receive a task from the master: 
    task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
    task_info <- mpi.get.sourcetag() 
    tag <- task_info[2]
    
    if (tag == 1){  #the received task is an actual task
      params.init <- task$params.init   #starting point for the optimizer
      lb.params <- task$lb.params       #lower boundaries for all parameters
      ub.params <- task$ub.params       #upper boundaries for all parameters
      algorithm <- task$algorithm       #optimization algorithm
      maxeval <- task$maxeval           #maximum number of iterations (-1 = no limit)
      #Perform optimization:
      params.mle.slave <- nloptr(x0=params.init,eval_f=neg_log_likelihood,lb=lb.params,ub=ub.params,
                                 opts=list("algorithm"=algorithm, maxeval=maxeval,"print_level" = 0))   #minimize -loglik <=> maximize loglik
      #Send result back to the master:
      result <- list(params.mle=params.mle.slave,params.init=params.init,algorithm=algorithm)
      mpi.send.Robj(result,0,2)   # tag=2 for 'task done'
    }
    else if (tag == 2){  #master says all tasks are done
      done <- TRUE
    } 
    #ignore messages with unknown tags
  }
  # When done, tell the master the slave is done
  mpi.send.Robj(junk,0,3)
}

################################################################################################################################################
import.data <- FALSE
N.species <- 3
dt.min <- 0.1   #minumum time step for numeric approximation of L-V eqs (irrelevant if ode.solver=='integrate')
# Gaussian noise sigma:
sigma <- 0.05

ode.solver = 'integrate'  #either 'integrate' (strongly recommended) or 'iterate'

N.starts <- 100    #number of random initial points in the parameter space to try. 
N.slaves <- 10     #number of slaves to use for solution search (each slave runs an independent optimization task with a random starting point)
################################################################################################################################################

if (import.data==TRUE)
{
  
  #Import observed population data for 3 species:
  setwd("~/Projects/Lotka-Volterra")
  # Y.series <- read.csv("LV variables with bitrate.csv", header=T)
  # Y.series <- read.csv("US_car_fleet_by_emission_control_systems.csv", header=T)
  Y.series <- read.csv("Smartphone OS market share USA.csv", header=T); Y.series[,1] = seq(1,nrow(Y.series));
  #   #Normalize Y.series to max(Y.series[,-1]), to avoid large values. This is equivalent to measuring X and K in units of max(Y.series[,-1])
  #   max.Y <- max(Y.series[,-1],na.rm=TRUE)
  #   Y.series[,-1] <- Y.series[,-1]/max.Y
  N.species <- ncol(Y.series)-1
  colnames(Y.series) <- c('t',names(Y.series[,-1]))
  Y.series.all <- Y.series  #this will be used to compare the fitted model's "forecast" with the actual data for t>t.end.Y
  t.series <- Y.series[,1]
  t.start <- min(t.series)   #t.start is the starting time of observed time series
  t.end <- max(t.series)
  dt <- min(min(diff(t.series)),dt.min)   #internal time step used to solve L-V eqs when ode.solver=='iterate' (irrelevant if ode.solver=='integrate')
  
  t.start.Y <- min(Y.series.all[complete.cases(Y.series.all),]$t)
  t.end.Y <- max(Y.series.all[complete.cases(Y.series.all),]$t)
  ind <- which(!(Y.series$t>=t.start.Y & Y.series$t<=t.end.Y))
  Y.series[ind,2:ncol(Y.series)] <- NA
  
  # Specify "switching-on" times for all species (time at which the specie is introduced into the system):
  t.switch.on <- list()
  for (specie in 1:N.species){
    t.switch.on[[specie]] <- t.start   #switch-on at the starting time (default)
  }
  t.switch.on[[3]] <- 1972   #Specie 3 is introduced into the system at t=1972
  
  plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,newplot=T)      #plot observed (noisy) population dynamics
  plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)
  # Show the range of data used for fitting:
  abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
  abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
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
  theta <- LV_params_to_theta(r,K,alpha)
  t.start <- 0   #t.start is the starting time of observed time series
  t.end <- 100   #t.end is the end time of observed time series
  t.series <- seq(t.start,t.end,by=1)    #set of times to generate time series data at. Doesn't have to be uniform
  N.series <- length(t.series)
  dt <- min(dt.min,min(diff(t.series))/10)   #internal time step used to solve L-V eqs
  X.series <- generate_X_time_series(theta,X0,t.set=t.series,dt,ode.solver)   #generate underlying population series (unobserved)
  # Generate observed population series Y.series:
  Y.series <- X.series
  #Add Gaussian observational noise:
  N <- dim(Y.series[,2:ncol(Y.series)])
  Y.series[,2:ncol(Y.series)] <- Y.series[,2:ncol(Y.series)] + rnorm(N[1]*N[2],0,sigma)   #this is the observed data (simulated)
  Y.series.all <- Y.series  #this will be used to compare the fitted model's "forecast" with the actual data for t>t.end.Y
  t.start.Y <- 0
  t.end.Y <- t.end
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
  cat("Global maximum of log_likelihood, evaluated at true parameter values:",-neg_log_likelihood(params.true),"\n")
  #Now that we have simulated the data using the specified parameters theta, sigma and X0, forget the parameters (they are to be inferred)
  sigma.true <- sigma
  X0.true <- X0
  rm(list=c('theta','sigma','X0'))
}

t.series <- Y.series$t
N.series <- sum(!is.na(Y.series[,2]))

# NLOPT optimizer:
# nlopt.algorithm <- "NLOPT_LN_NELDERMEAD"
nlopt.algorithm <- "NLOPT_LN_COBYLA"
# nlopt.algorithm <- "NLOPT_LN_BOBYQA"
# nlopt.algorithm <- "NLOPT_LN_SBPLX"
# nlopt.algorithm <- "NLOPT_GN_DIRECT_L"   #global optimization algorithm, does not require gradients

maxeval=1e6   #maximum number of iterations for nloptr. If -1, the number of iterations is unlimited

# lower and upper bounds on parameters:
#r:
max.value <- 5
eps <- 1e-5
lb.r <- rep(eps,N.species)
ub.r <- rep(max.value,N.species)
#K:
lb.K <- rep(eps,N.species)
ub.K <- rep(max.value,N.species)
#alpha:
lb.alpha <- rep(-max.value,N.species^2-N.species)
ub.alpha <- rep(max.value,N.species^2-N.species)
#sigma:
lb.sigma <- rep(eps,1)
ub.sigma <- rep(max.value,1)
#X0:
lb.X0 <- rep(eps,N.species)
ub.X0 <- rep(max.value,N.species)
#overall lb and ub:
lb.params <- c(lb.r,lb.K,lb.alpha,lb.sigma,lb.X0)
ub.params <- c(ub.r,ub.K,ub.alpha,ub.sigma,ub.X0)

# Global optimizer is too slow, so use a local optimizer with a number of random starting points, and then choose the smallest of the 
# obtained minima, hoping it will be the global minimum:
mpi.spawn.Rslaves(nslaves=N.slaves)   #spawn as many slaves as possible
# Push all the required data and functions to the slaves:
mpi.bcast.Robj2slave(all=TRUE)  #send all objects (data and functions) from master's global environment to slaves
# Send the task handling function slave_optimize() to the slaves:
mpi.bcast.Robj2slave(slave_optimize)
# Call the function slave_optimize() in all the slaves to get them ready to undertake tasks:
mpi.bcast.cmd(slave_optimize())

# Create task list:
tasks <- vector('list')
for (attempt in 1:N.starts){
  # Random starting point in the parameter space, for optimizer to start at:
  r.init <- runif(N.species,min=eps,max=max.value)
  K.init <- runif(N.species,min=eps,max=2*max(Y.series.all[,-1],na.rm=TRUE))
  alpha.init <- matrix(runif(N.species*N.species,min=0,max=max.value), nrow=N.species, ncol=N.species)
  theta.init <- LV_params_to_theta(r.init,K.init,alpha.init)
  sigma.init <- runif(1,min=eps,max=max.value)
  # X0.init <- abs(Y.series.all[1,-1])
  X0.init <- rep(eps,N.species)
  names(X0.init) <- names(Y.series[,-1])
  params.init <- pack_parameters(theta.init,sigma.init,X0.init)
  #Check that params.init is within the constraints:
  if (!all(params.init >= lb.params & params.init <= ub.params)) stop("Initial parameters are outside of the hyperbox! Stopping.")
  tasks[[attempt]] <- list(params.init=params.init,lb.params=lb.params,ub.params=ub.params,algorithm=nlopt.algorithm,maxeval=maxeval)
}
junk <- 0
finished_slaves <- 0
n_slaves <- mpi.comm.size()-1
# list to store the results returned by slaves:
params.mle <- list()
min.objective <- Inf
attempt <- 1   #counts the optimization attempts, there should be N.starts attempts in total.
while (finished_slaves < n_slaves){
  # Receive a message from a slave:
  message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
  message_info <- mpi.get.sourcetag() 
  slave_id <- message_info[1] 
  tag <- message_info[2]
  if (tag == 1){
    # slave is ready for a task. Give it the next task, or tell it tasks 
    cat("Slave",slave_id,"is ready for a new task...")
    # are done if there are none.
    if (length(tasks)>0){
      # Send a task, and then remove it from the task list:
      mpi.send.Robj(tasks[[1]], slave_id, 1); 
      tasks[[1]] <- NULL 
      cat("assigned a new task to slave",slave_id,".",length(tasks),"tasks remaining\n")
    }
    else {
      mpi.send.Robj(junk, slave_id, 2)   #let the slave know there are no more tasks (tag=2)
      cat("but there are no more tasks for slave",slave_id,"\n")
    }
  }
  else if (tag==2){
    # The message received from the slave contains results.
    params.mle[[attempt]] <- message$params.mle
    if (params.mle[[attempt]]$objective < min.objective){
      params.mle.best <- params.mle[[attempt]]
      min.objective <- params.mle.best$objective
    }
    cat("Slave",slave_id,"returned a newly found minimum:",message$params.mle$objective,". Best minimum found so far:",min.objective,"\n")
    attempt <- attempt+1
  }
  else if (tag==3) {
    #The slave has finished.
    cat("Slave",slave_id,"finished work.\n")
    finished_slaves <- finished_slaves + 1
  }
}
mpi.close.Rslaves()  #close all slaves

#Extract MLEs of parameters from the best solution found:
print(params.mle.best)   #print a summary of the best found solution
solution <- params.mle.best$solution
theta.mle <- unpack_parameters(solution)$theta
tmp <- theta_to_LV_params(theta.mle)
r.mle <- tmp$r
K.mle <- tmp$K
alpha.mle <- tmp$alpha
sigma.mle <- unpack_parameters(solution)$sigma
X0.mle <- unpack_parameters(solution)$X0

X.series.nlopt <- generate_X_time_series(theta.mle,as.data.frame(t(X0.mle)),t.set=t.series,dt,ode.solver)   #generate population series for the inferred parameters
if (import.data==FALSE) {
  plot.dynamics(pop.data=X.series,plot.type='l',newplot=T)   #plot the actual underlying (unobserved) population dynamics
} else plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,newplot=T)
plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)      #plot observed (noisy) population dynamics
abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
plot.dynamics(pop.data=X.series.nlopt,plot.type='l',line.type="dashed",line.width=3,newplot=F)  #plot the inferred population dynamics
title(main="Best solution found by nlopt")

# Refine the best found minimum:
params.init <- solution
#optim optimizer:
max.iterations <- 1e4
rel.tol = sqrt(.Machine$double.eps)
# optim.algorithm <- "SANN"
optim.algorithm <- "BFGS"
# optim.algorithm <- "L-BFGS-B"
# optim.algorithm <- "Nelder-Mead"
if (optim.algorithm == "L-BFGS-B"){
  params.mle.optim <- optim(par=params.init,fn=neg_log_likelihood,method=optim.algorithm,lower=lb.params,upper=ub.params,
                            control=list(maxit=max.iterations,reltol=rel.tol,trace=TRUE,REPORT=1))
} else {
  params.mle.optim <- optim(par=params.init,fn=neg_log_likelihood,method=optim.algorithm,
                            control=list(maxit=max.iterations,reltol=rel.tol,trace=TRUE,REPORT=1))
}
solution <- params.mle.optim$par
theta.mle.optim <- unpack_parameters(solution)$theta
tmp <- theta_to_LV_params(theta.mle.optim)
r.mle.optim <- tmp$r
K.mle.optim <- tmp$K
alpha.mle.optim <- tmp$alpha
sigma.mle.optim <- unpack_parameters(solution)$sigma
X0.mle.optim <- unpack_parameters(solution)$X0

X.series.optim <- generate_X_time_series(theta.mle.optim,as.data.frame(t(X0.mle.optim)),t.set=t.series,dt,ode.solver)   #generate population series for the inferred parameters
if (import.data==FALSE) {
  plot.dynamics(pop.data=X.series,plot.type='l',newplot=T)   #plot the actual underlying (unobserved) population dynamics
} else plot.dynamics(pop.data=Y.series.all,plot.type='p',pch=1,newplot=T)
plot.dynamics(pop.data=Y.series,plot.type='p',pch=19,newplot=F)      #plot observed (noisy) population dynamics
abline(v=t.start.Y,lty='dotted',lwd=2,col='red')
abline(v=t.end.Y,lty='dotted',lwd=2,col='red')
plot.dynamics(pop.data=X.series.optim,plot.type='l',line.type="dashed",line.width=3,newplot=F)  #plot the inferred population dynamics
legend('right',names(Y.series[,-1]),lty=1, lwd=2, col=seq(1,N.species))
title(main="Best solution found by nlopt, refined by optim")

# mpi.quit(save="no")

