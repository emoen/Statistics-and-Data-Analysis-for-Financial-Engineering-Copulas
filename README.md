Statatistics and Data Analysis for Financial Engineering - Chapter 8
Copulas
================

## Archimedean Copulas

Archimedean copula with generator function:

![
C(u_1,..,u_d) = \\phi^{-1}(\\phi(u_1)+..+\\phi(u_d))
](https://latex.codecogs.com/png.latex?%0AC%28u_1%2C..%2Cu_d%29%20%3D%20%5Cphi%5E%7B-1%7D%28%5Cphi%28u_1%29%2B..%2B%5Cphi%28u_d%29%29%0A "
C(u_1,..,u_d) = \phi^{-1}(\phi(u_1)+..+\phi(u_d))
")

#### Frank copula

Generator function:

![
C(u\|\\theta) = -ln(\\frac{e^{-\\theta u} -1)}{e^{-\\theta} -1)}), -\\infty \< \\theta \< \\infty
](https://latex.codecogs.com/png.latex?%0AC%28u%7C%5Ctheta%29%20%3D%20-ln%28%5Cfrac%7Be%5E%7B-%5Ctheta%20u%7D%20-1%29%7D%7Be%5E%7B-%5Ctheta%7D%20-1%29%7D%29%2C%20-%5Cinfty%20%3C%20%5Ctheta%20%3C%20%5Cinfty%0A "
C(u|\theta) = -ln(\frac{e^{-\theta u} -1)}{e^{-\theta} -1)}), -\infty < \theta < \infty
")

#### Crayton copula

Generator function:

![
C(u\|\\theta) = \\frac{1}{\\theta}(u^{-\\theta}-1), \\theta > 0
](https://latex.codecogs.com/png.latex?%0AC%28u%7C%5Ctheta%29%20%3D%20%5Cfrac%7B1%7D%7B%5Ctheta%7D%28u%5E%7B-%5Ctheta%7D-1%29%2C%20%5Ctheta%20%3E%200%0A "
C(u|\theta) = \frac{1}{\theta}(u^{-\theta}-1), \theta > 0
")

#### Gumbel Copula adsfads

Generator function:

![
C(u\|\\theta) = -ln(u)^{\\theta}, \\theta \\geq 1
](https://latex.codecogs.com/png.latex?%0AC%28u%7C%5Ctheta%29%20%3D%20-ln%28u%29%5E%7B%5Ctheta%7D%2C%20%5Ctheta%20%5Cgeq%201%0A "
C(u|\theta) = -ln(u)^{\theta}, \theta \geq 1
")

#### Joe copula

Generator function:

![
C(u\|\\theta) = -ln(1-(1-\\theta)^{\\theta}), \\theta \\geq 1
](https://latex.codecogs.com/png.latex?%0AC%28u%7C%5Ctheta%29%20%3D%20-ln%281-%281-%5Ctheta%29%5E%7B%5Ctheta%7D%29%2C%20%5Ctheta%20%5Cgeq%201%0A "
C(u|\theta) = -ln(1-(1-\theta)^{\theta}), \theta \geq 1
")

### Plot of generator function for Frank cupola

``` r
library(copula)
u= seq(0.000001, 1, length=500)
frank = iPsi(copula=archmCopula(family="frank", param=1), u)
plot(u, frank, type="l", lwd=3, ylab=expression(phi(u)))
abline(h=0)
abline(v=0)
```

![](readme_files/figure-gfm/Frank%20Copula%20generator%20function-1.png)<!-- -->

### Scatter plot of 9 bivariate Frank copulas

``` r
set.seed(5640)
theta = c(-100, -50, -10, -1, 0, 5, 20, 50, 500)
par(mfrow=c(3,3), cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
for(i in 1:9){
  U= rCopula(n=200, copula=archmCopula(family="frank", param=theta[i]))
  plot(U, xlab=expression(u[1]), ylab=expression(u[2]), 
                          main=eval(substitute(expression(paste(theta, " = ", j)),
                         list(j = as.character(theta[i])))))
}
```

    ## parameter at boundary ==> returning indepCopula()

![](readme_files/figure-gfm/Scatterplot%20of%209%20frank%20copulas-1.png)<!-- -->
Figure 2: Bivariate random samples of size 200 from various Frank
copulas.

### Scatter plot of 9 bivariate Clayton copulas

``` r
set.seed(5640)
theta = c(-0.98, -0.7, -0.3, -0.1, 0.1, 1, 5, 15, 100)
par(mfrow=c(3,3), cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
for(i in 1:9){
  U= rCopula(n=200, copula=archmCopula(family="clayton", param=theta[i]))
  plot(U, xlab=expression(u[1]), ylab=expression(u[2]), 
                          main=eval(substitute(expression(paste(theta, " = ", j)),
                         list(j = as.character(theta[i])))))
}
```

![](readme_files/figure-gfm/Scatterplot%20of%209%20clayton%20copulas-1.png)<!-- -->
Figure 2: Bivariate random samples of size 200 from various clayton
copulas.
