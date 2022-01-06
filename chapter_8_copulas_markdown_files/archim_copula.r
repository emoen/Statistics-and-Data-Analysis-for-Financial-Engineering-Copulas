# Plot of generator function for frank cupola
library(copula)
u= seq(0.000001, 1, length=500)
frank = iPsi(copula=archmCopula(family="frank", param=1), u)
plot(u, frank, type="l", lwd=3, ylab=expression(phi(u)))
abline(h=0)
abline(v=0)

## Scatter plot of 9 bivariate frank copulas
set.seed(5640)
theta = c(-100, -50, -10, -1, 0, 5, 20, 50, 500)
par(mfrow=c(3,3), cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
for(i in 1:9){
  U= rCopula(n=200, copula=archmCopula(family="frank", param=theta[i]))
  plot(U, xlab=expression(u[1]), ylab=expression(u[2]), 
                          main=eval(substitute(expression(paste(theta, " = ", j)),
                         list(j = as.character(theta[i])))))
}