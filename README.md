Statatistics and Data Analysis for Financial Engineering - Chapter 8
Copulas
================

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

### Plot of generator function for frank cupola

``` r
library(copula)
u= seq(0.000001, 1, length=500)
frank = iPsi(copula=archmCopula(family="frank", param=1), u)
plot(u, frank, type="l", lwd=3, ylab=expression(phi(u)))
abline(h=0)
abline(v=0)
```

![](readme_files/figure-gfm/Frank%20Copula%20generator%20function-1.png)<!-- -->

## Including Plots

You can also embed plots, for example:

![](readme_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
