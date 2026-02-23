<!-- README.md is generated from README.Rmd. Please edit that file -->

# zoomgrid version 1.1.0 (Red Grid)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/zoomgrid?color=green)](https://cran.r-project.org/package=zoomgrid)
![](http://cranlogs.r-pkg.org/badges/grand-total/zoomgrid?color=green)
![](http://cranlogs.r-pkg.org/badges/zoomgrid?color=green)
![](http://cranlogs.r-pkg.org/badges/last-week/zoomgrid?color=green)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-green.svg)](https://www.gnu.org/licenses/gpl-3.0)

The package implements the grid search algorithm with a zoom. The grid
search algorithm with a zoom aims to help solving difficult optimization
problem where there are many local optimisers inside the domain of the
target function. It offers suitable initial or starting value for the
following optimization procedure, provided that the global optimum
exists in the neighbourhood of the initial or starting value. The grid
search algorithm with a zoom saves time tremendously in cases with
high-dimensional arguments.

You can also find the package on CRAN, see

[zoomgrid@CRAN](https://CRAN.R-project.org/package=zoomgrid)

and the corresponding paper

[Modelling Nonlinear Vector Economic Time
Series](https://pure.au.dk/ws/files/45638557/Yukai_Yang_PhD_Thesis.pdf)

See section 1.5.4.

## How to install

You can either install the stable version from CRAN

``` r
install.packages("zoomgrid")
```

or install the development version from GitHub

``` r
devtools::install_github("yukai-yang/zoomgrid")
```

provided that the package “devtools” has been installed beforehand.

## Example

After installing the package, you need to load (attach better say) it by
running the code

``` r
library(zoomgrid)
```

You can take a look at all the available functions and data in the
package

``` r
ls("package:zoomgrid")
#> [1] "build_grid"        "grid_search"       "grid_search_check"
```

### Motivation

Consider the two-dimensional **Rastrigin function**, which is a
non-convex function widely used for testing optimisation algorithms.

![](man/figures/gif.latex.gif)

where *x*<sub>*i*</sub> ∈ \[−5.12, 5.12\] and *A* = 10. It has many
local minima and its global minimum is at (0, 0) with the minimum value
0.

<figure>
<img src="man/figures/Rastrigin_function.png"
alt="Diegotorquemada [Public domain], from Wikimedia Commons" />
<figcaption aria-hidden="true">Diegotorquemada [Public domain], from
Wikimedia Commons</figcaption>
</figure>

Graph source: [Rastrigin function @
WIKIPEDIA](https://en.wikipedia.org/wiki/Rastrigin_function).

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

Then let us try the optimization algorithms available in the `optim`
function.

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
#> [1] -1.989912e+00  2.913342e-09
#> [1] 3.979831
tmp4$par; tmp4$value
#> [1] 0.97915333 0.01486102
#> [1] 1.088185
```

None of them are satisfactory…

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
#> 
#> ── zoomgrid version 1.1.0 (Red Grid) ───────────────────────────────────────────
#> ✔ Grid search with 0 zoom-in layers and 1 point each produced 1 result.
#> The minimiser is believed to be in the neighbourhood of -0.0199999999999996 and
#> -0.0199999999999996.
#> ℹ Elapsed: 3.837s (user: 3.798s, system: 0.021s).
ret1$par
#> [1] -0.02 -0.02
```

Then we run the parallel one. Parallel execution uses the future
framework and works on all major platforms including Windows.

``` r
# parallel computation
ret2 = grid_search(Rastrigin, grid, num=2, parallel=TRUE, silent=FALSE)
#> 
#> ── zoomgrid version 1.1.0 (Red Grid) ───────────────────────────────────────────
#> ℹ Parallel computation runs with 2 workers.
#> ✔ Grid search with 0 zoom-in layers and 2 points each produced 2 results.
#> The minimiser is believed to be in the neighbourhood of -0.0199999999999996 and
#> -0.0199999999999996.
#> ℹ Elapsed: 2.605s (user: 0.606s, system: 0.014s).
ret2$par
#> [1] -0.02 -0.02
```

Try the grid search with a zoom!

``` r
# grid search with a zoom!
ret3 = grid_search(Rastrigin, grid, zoom=2, num=2, parallel=TRUE, silent=FALSE)
#> 
#> ── zoomgrid version 1.1.0 (Red Grid) ───────────────────────────────────────────
#> ℹ Parallel computation runs with 2 workers.
#> ✔ Grid search with 2 zoom-in layers and 2 points each produced 14 results.
#> The minimiser is believed to be in the neighbourhood of 5.59049615653446e-05
#> and 5.59049615653446e-05.
#> ℹ Elapsed: 4.855s (user: 1.277s, system: 0.041s).
ret3$par
#> [1] 5.590496e-05 5.590496e-05
```

Sometimes it is strongly recommended to check the time consumed by
running the grid search first. This is extremely useful when the user is
going to run on some super-computing server and need to know
approximately how long it will take in order to specify the
corresponding settings according to some batch system like SLURM for
example. So you can do as follows

``` r
ret3 = grid_search_check(Rastrigin, grid, zoom=2, num=2, parallel=TRUE, silent=FALSE)
#> 
#> ── zoomgrid version 1.1.0 (Red Grid) ───────────────────────────────────────────
#> ℹ Parallel computation runs with 2 workers.
#> ℹ The expected time consumed by running the grid search is around 4.61959 seconds.
ret3 = grid_search(Rastrigin, grid, zoom=2, num=2, parallel=TRUE, silent=FALSE)
#> 
#> ── zoomgrid version 1.1.0 (Red Grid) ───────────────────────────────────────────
#> ℹ Parallel computation runs with 2 workers.
#> ✔ Grid search with 2 zoom-in layers and 2 points each produced 14 results.
#> The minimiser is believed to be in the neighbourhood of 5.59049615653446e-05
#> and 5.59049615653446e-05.
#> ℹ Elapsed: 4.935s (user: 1.317s, system: 0.042s).
```

## Citation

If you use the zoomgrid package in your research, please cite the
underlying methodology.

### Methodology

Yang, Y. (2012). Modelling Nonlinear Vector Economic Time Series. PhD
thesis, Aarhus University, Department of Economics and Business,
CREATES. Available at:
<https://pure.au.dk/ws/files/428886102/Yukai_Yang_PhD_Thesis.pdf>

You can obtain the citation information directly from R by running:

``` r
citation("zoomgrid")
#> To cite package 'zoomgrid' in publications use:
#> 
#>   Yang Y (2012). _Modelling Nonlinear Vector Economic Time Series_.
#>   Ph.D. thesis, Aarhus University, Aarhus, Denmark. PhD thesis,
#>   Department of Economics and Business, CREATES,
#>   <https://pure.au.dk/ws/files/428886102/Yukai_Yang_PhD_Thesis.pdf>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @PhdThesis{,
#>     title = {Modelling Nonlinear Vector Economic Time Series},
#>     author = {Yukai Yang},
#>     year = {2012},
#>     school = {Aarhus University},
#>     address = {Aarhus, Denmark},
#>     note = {PhD thesis, Department of Economics and Business, CREATES},
#>     url = {https://pure.au.dk/ws/files/428886102/Yukai_Yang_PhD_Thesis.pdf},
#>   }
```
