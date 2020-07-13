log_uposterior_slice <- function(params,slice.dimension,slice.range.lb=NA,slice.range.ub=NA,n.points=100L)   
# Plot a 1-d slice of log_uposterior along slice.dimension, 
# taking the rest of parameters from the 'params' argument
{
  stopifnot(is.numeric(slice.dimension))
  stopifnot(all(slice.dimension <= length(params)))
  n.plots <- length(slice.dimension)
  
  if (n.plots==1) par(mfrow=c(1,1)) else par(mfrow=c(ceiling(n.plots/3),3))
  count <- 1
  for (d in slice.dimension){
    if ((is.na(slice.range.lb[count]) && (is.na(slice.range.ub[count])))){
      slice.range.lb[count] <- lb.params[d]
      slice.range.ub[count] <- ub.params[d]
    }
    slice.var.step <- (slice.range.ub[count]-slice.range.lb[count])/n.points
    if (slice.var.step<=0) stop("slice.var.step should be >0! Stopping.")
    slice.frame <- data.frame('slice.var'=rep(NA,(n.points+1)),'log_uposterior'=rep(NA,(n.points+1)))
    for (i in 1:(n.points+1)){
      slice.frame[i,'slice.var'] <- slice.range.lb[count] + (i-1)*slice.var.step
      slice.params <- params
      slice.params[d] <- slice.frame[i,'slice.var']
      slice.frame[i,'log_uposterior'] <- log_uposterior(slice.params,lb.params,ub.params)
    }
    plot(slice.frame,type='l',col='blue',xlab=paste("parameter",d),ylab="log-uposterior")
    abline(v=params[d],lty='dotted',col='red')
    count <- count + 1
  }
}