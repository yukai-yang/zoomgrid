<!-- README.md is generated from README.Rmd. Please edit that file -->
zoomgrid version 1.0.0 (Red Grid)
=================================

The package implements provides the grid search algorithm with a zoom. The grid search algorithm with a zoom aims to help solving difficult optimization problem where there are many local optimizers inside the domain of the target function. It offers suitable initial or starting value for the following optimization procedure, provided that the global optimum exists in the neighbourhood of the initial or starting value. The grid search algorithm with a zoom saves time tremendously in cases with high-dimenstional arguments.

and the corresponding paper

[Modelling Nonlinear Vector Economic Time Series](https://pure.au.dk/ws/files/45638557/Yukai_Yang_PhD_Thesis.pdf)

See section 1.5.4.

How to install
--------------

You can install the development version from GitHub

``` r
devtools::install_github("yukai-yang/zoomgrid")
```

provided that the package "devtools" has been installed beforehand.

Example
-------

After installing the package, you need to load (attach better say) it by running the code

``` r
library(zoomgrid)
```

You can take a look at all the available functions and data in the package

``` r
ls( grep("zoomgrid", search()) ) 
#> [1] "build_grid"  "grid_search"
```

### Motivation

Consider the two-dimensional **Rastrigin function** is a non-convex function which is widely used for testing the performances of some optimization algorithms.

$$f(x\_1, x\_2) = A \\cdot n + \\sum\_{i=1}^2 \\left( x\_i^2 - A \\cos(2 \\pi x\_i) \\right)$$

where *x*<sub>*i*</sub> ∈ \[ − 5.12, 5.12\]. It has many local minimum and its global minimum is at (0, 0) with the minimum value 0. ![](Rastrigin_function.png)

Graph source: [Rastrigin function @ WIKIPEDIA](https://en.wikipedia.org/wiki/Rastrigin_function). We give the function in R:

``` r
# Rastrigin function
ndim = 2 # number of dimension
nA = 10 # parameter A
# vx in [-5.12, 5.12]

# minimizer = rep(0, ndim)
# minimum = 0
Rastrigin <- function(vx) return(nA * ndim + sum(vx*vx - nA * cos(2*pi*vx)))
```

Then let us try the optimization algorithms available in the `optim` function.

``` r
# set seed and initialize the initial or starting value
set.seed(1)
par = runif(ndim, -5.12, 5.12)
cat("start from", par)
#> start from -2.401191 -1.309451

# results from different optimization algorithms
tmp1 = optim(par = par, Rastrigin, method='Nelder-Mead')
tmp2 = optim(par = par, Rastrigin, method='BFGS')
tmp3 = optim(par = par, Rastrigin, method='L-BFGS-B')
tmp4 = optim(par = par, Rastrigin, method='SANN')

tmp1$par; tmp1$value
#> [1] -1.9899136 -0.9949483
#> [1] 4.97479
tmp2$par; tmp2$value
#> [1] -0.9949586  0.9949586
#> [1] 1.989918
tmp3$par; tmp3$value
#> [1] -1.989912e+00  2.913338e-09
#> [1] 3.979831
tmp4$par; tmp4$value
#> [1] 0.97915333 0.01486102
#> [1] 1.088185
```

None of them are satisfactory...

### Build the grid

We need to build grid first for the grid search. For details, see

``` r
?build_grid
```

We build the grid by running

``` r
# build the grid
bin = c(from=-5.12, to=5.12, by=.1)
grid = build_grid(bin,bin)
```

### Grid search

We can first try the sequential (no parallel) grid search

``` r
# serial computation
ret1 = grid_search(Rastrigin, grid, silent=FALSE)
#> ###########################################################################
#> zoomgrid version 1.0.0 (Red Grid)
#> ---------------------------------------------------------------------------
#> The Grid Search of 0 zoom-in layers with 1 points each gives 1 results.
#> The minimizer is believed to be in the neighbourhood of -0.02 -0.02.
#> ---------------------------------------------------------------------------
#>    user  system elapsed 
#>   8.701   0.027   8.764 
#> ###########################################################################
ret1$par
#> [1] -0.02 -0.02
```

Then we run the parallel one

``` r
# parallel computation
ret2 = grid_search(Rastrigin, grid, num=2, parallel=TRUE, silent=FALSE)
#> ###########################################################################
#> zoomgrid version 1.0.0 (Red Grid)
#> ---------------------------------------------------------------------------
#> Parallel computation runs with 4 cores.
#> The Grid Search of 0 zoom-in layers with 2 points each gives 2 results.
#> The minimizer is believed to be in the neighbourhood of -0.02 -0.02.
#> ---------------------------------------------------------------------------
#>    user  system elapsed 
#>   15.38    0.29    4.72 
#> ###########################################################################
ret2$par
#> [1] -0.02 -0.02
```

Try the grid search with a zoom!

``` r
# grid search with a zoom!
ret3 = grid_search(Rastrigin, grid, zoom=2, num=2, parallel=TRUE, silent=FALSE)
#> ###########################################################################
#> zoomgrid version 1.0.0 (Red Grid)
#> ---------------------------------------------------------------------------
#> Parallel computation runs with 4 cores.
#> The Grid Search of 2 zoom-in layers with 2 points each gives 14 results.
#> The minimizer is believed to be in the neighbourhood of 5.590496e-05 5.590496e-05.
#> ---------------------------------------------------------------------------
#>    user  system elapsed 
#>  26.195   1.603   7.687 
#> ###########################################################################
ret3$par
#> [1] 5.590496e-05 5.590496e-05
```
