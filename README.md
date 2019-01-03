<!-- README.md is generated from README.Rmd. Please edit that file -->
zoomgrid version 1.0.0 (Red Grid)
=================================

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/zoomgrid?color=green)](https://cran.r-project.org/package=zoomgrid) ![](http://cranlogs.r-pkg.org/badges/grand-total/zoomgrid?color=green) ![](http://cranlogs.r-pkg.org/badges/zoomgrid?color=green) ![](http://cranlogs.r-pkg.org/badges/last-week/zoomgrid?color=green)

The package implements provides the grid search algorithm with a zoom. The grid search algorithm with a zoom aims to help solving difficult optimization problem where there are many local optimizers inside the domain of the target function. It offers suitable initial or starting value for the following optimization procedure, provided that the global optimum exists in the neighbourhood of the initial or starting value. The grid search algorithm with a zoom saves time tremendously in cases with high-dimenstional arguments.

You can also find the package on CRAN, see

[zoomgrid@CRAN](https://CRAN.R-project.org/package=zoomgrid)

and the corresponding paper

[Modelling Nonlinear Vector Economic Time Series](https://pure.au.dk/ws/files/45638557/Yukai_Yang_PhD_Thesis.pdf)

See section 1.5.4.

How to install
--------------

You can either install the stable version from CRAN

``` r
install.packages("zoomgrid")
```

or install the development version from GitHub

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
#> [1] "build_grid"        "grid_search"       "grid_search_check"
```

### Motivation

Consider the two-dimensional **Rastrigin function** is a non-convex function which is widely used for testing the performances of some optimization algorithms.

![](https://latex.codecogs.com/gif.latex?f%28x_1%2C%20x_2%29%20%3D%202%20A%20+%20%5Csum_%7Bi%3D1%7D%5E2%20%5Cleft%28%20x_i%5E2%20-%20A%20%5Ccos%282%20%5Cpi%20x_i%29%20%5Cright%29)

where *x*<sub>*i*</sub> ∈ \[ − 5.12, 5.12\] and *A* = 10. It has many local minimum and its global minimum is at (0, 0) with the minimum value 0. <a title="Diegotorquemada [Public domain], from Wikimedia Commons" href="https://commons.wikimedia.org/wiki/File:Rastrigin_function.png"><img alt="Rastrigin function" src="https://upload.wikimedia.org/wikipedia/commons/thumb/8/8b/Rastrigin_function.png/512px-Rastrigin_function.png"></a>

Graph source: [Rastrigin function @ WIKIPEDIA](https://en.wikipedia.org/wiki/Rastrigin_function).

We give the function in R:

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
#>   8.621   0.039   8.713 
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
#>  14.504   0.437   4.159 
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
#>  27.881   1.945   9.325 
#> ###########################################################################
ret3$par
#> [1] 5.590496e-05 5.590496e-05
```

Sometimes it is strongly recommended to check the time consumed by running the grid search first. This is extremely useful when the user is going to run on some super-computing server and need to know approximately how long time it will take in order to specify the corresponding settings according to some batch system like SLURM for example. So you can do as follows

``` r
ret3 = grid_search_check(Rastrigin, grid, zoom=2, num=2, parallel=TRUE, silent=FALSE)
#> ###########################################################################
#> zoomgrid version 1.0.0 (Red Grid)
#> ---------------------------------------------------------------------------
#> Parallel computation runs with 4 cores.
#> The expected time consumed by running the grid search is around 9.557705 seconds.
#> ###########################################################################
ret3 = grid_search(Rastrigin, grid, zoom=2, num=2, parallel=TRUE, silent=FALSE)
#> ###########################################################################
#> zoomgrid version 1.0.0 (Red Grid)
#> ---------------------------------------------------------------------------
#> Parallel computation runs with 4 cores.
#> The Grid Search of 2 zoom-in layers with 2 points each gives 14 results.
#> The minimizer is believed to be in the neighbourhood of 5.590496e-05 5.590496e-05.
#> ---------------------------------------------------------------------------
#>    user  system elapsed 
#>  26.975   1.591  10.203 
#> ###########################################################################
```
