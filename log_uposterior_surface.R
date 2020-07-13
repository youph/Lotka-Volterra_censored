log_uposterior_surface <- function(params,slice.dimension,slice.range.lb=NULL,slice.range.ub=NULL,n.gridpoints=10L,n.cuts=20,
                                   zlim=NULL,add.contours=TRUE,x.name=paste("parameter",slice.dimension[1]),y.name=paste("parameter",slice.dimension[2]))   
  # Plot a 2-d surface of log_uposterior vs x=slice.dimension[1], y=slice.dimension[2],
  # taking the rest of parameters from the 'params' argument
{
  require(lattice)
  stopifnot(is.numeric(slice.dimension))
  stopifnot(all(slice.dimension <= length(params)))
  stopifnot(length(slice.dimension)==2)
  
  if (!(is.numeric(slice.range.lb) && (is.numeric(slice.range.ub)))){
    slice.range.lb <- lb.params[slice.dimension]
    slice.range.ub <- ub.params[slice.dimension]
  }
  
  par(mfrow=c(1,1))
  x <- seq(slice.range.lb[1], slice.range.ub[1], length.out=n.gridpoints+1)
  y <- seq(slice.range.lb[2], slice.range.ub[2], length.out=n.gridpoints+1)
  grid <- expand.grid(x=x,y=y)
  z <- matrix(NA,nrow=length(x),ncol=length(y))
  slice.params <- params
  ind <- 1
  for (i in 1:length(x)){
    slice.params[slice.dimension[1]] <- x[i]
    for (j in 1:length(y)){
      slice.params[slice.dimension[2]] <- y[j]
      z[i,j] <- log_uposterior(slice.params,lb.params,ub.params)
      grid$z[ind] <- z[j,i]
      ind <- ind + 1
    }
  }
  
  #   contourplot(z~x*y,grid,xlab=x.name,ylab=y.name,zlab="log-uposterior", cuts=n.cuts, colorkey = TRUE, region = TRUE)
  # Draw the point with x=params[slice.dimension[1]], y=params[slice.dimension[2]]:
  #   trellis.focus("panel", 1, 1, highlight=T) 
  #   lpoints(x=params[slice.dimension[1]],y=params[slice.dimension[2]],pch=19,col='red') 
  #   trellis.unfocus() 
  
  if (!is.numeric(zlim)) zlim=range(z, finite=TRUE)
  if (add.contours) {
    filled.contour(x=x,y=y,z=z, nlevels=n.cuts,xlab=x.name,ylab=y.name,main="log-uposterior",zlim=zlim,color.palette = topo.colors,
                   plot.axes = {contour(x,y,z,zlim=zlim,add=TRUE); axis(1); axis(2); points(params[slice.dimension[1]], params[slice.dimension[2]],pch=19,col='red') })
  } else {
    filled.contour(x=x,y=y,z=z, nlevels=n.cuts,xlab=x.name,ylab=y.name,main="log-uposterior",zlim=zlim,color.palette = topo.colors,
                   plot.axes = {axis(1); axis(2); points(params[slice.dimension[1]], params[slice.dimension[2]],pch=19,col='red') })
  }
}