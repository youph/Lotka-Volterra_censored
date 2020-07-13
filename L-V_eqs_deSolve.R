library(deSolve)

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*(alpha - beta*y)
    dy = -y*(gamma - delta*x)
    return(list(c(dx, dy)))
  })
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

out <- as.data.frame(ode(func = LVeqs, y = as.numeric(X0), parms = theta, times = t.series, method="ode45"))

matplot(out[,-1], type = "l", xlab = "time", ylab = "population")
legend("topright", c("Cute bunnies", "Rabid foxes"), lty = c(1,2), col = c(1,2), box.lwd = 0)



Pars <- c(alpha = 2, beta = .5, gamma = .2, delta = .6)
State <- c(x = 10, y = 10)
