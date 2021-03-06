% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithms.R
\name{grid_search}
\alias{grid_search}
\title{Carry out the grid search algorithm with a zoom.}
\usage{
grid_search(FUN, grid, MoreArgs = NULL, zoom = 0, decay = 0.5,
  num = 1, parallel = FALSE, cores = NULL, silent = TRUE)
}
\arguments{
\item{FUN}{the target function to be minimized.}

\item{grid}{an object of the class GRID from \code{\link{build_grid}}.}

\item{MoreArgs}{a list of other arguments to \code{FUN}, see \code{\link{mapply}}.}

\item{zoom}{number of (additional) rounds or layers of the zoom-in, 0 by default.}

\item{decay}{a number in between 0 and 1 representing the decay rate of the grid sizes of the zoom.}

\item{num}{number of points to return, i.e. the smallest \code{num} points, 1 by default the minimum.}

\item{parallel}{a boolean indicating if the parallel computation is carried out, by default \code{FALSE}.}

\item{cores}{The number of cores to use, i.e. at most how many child processes will be run simultaneously. For details, see \code{mcmapply} in \code{parallel} package.}

\item{silent}{a boolean indicating if the information regarding the omputation is printed.}
}
\value{
a list containing the results from the grid search with a zoom.

The list contains the following components:
\item{par}{the approximate global minimizer}
\item{points}{all the local minimizer points found by the grid search with a zoom}
}
\description{
This function carries out the grid search algorithm with a zoom.
}
\details{
The target function \code{FUN} to be minimized is a scalar real valued function with multiple arguments.
The maximization can be achieved by multiplying -1 to the original function,
and then input the new function to \code{FUN}.

The \code{grid} must be created by using the function \code{\link{build_grid}} function in the package.

Any other invariant arugments to the function \code{FUN} can be specified in \code{MoreArgs} by using a list with variable names.

The common grid search first build a grid within some bounded area.
And then the target function will be evaluated at each point in the grid.
The points that produce the smallest \code{num} target function values shall be returned.

\code{zoom = 0} implies that no zoom-in will be applied and the corresponding grid search is the common one introduced above,
while any integer \code{zoom > 0} implies that the zoom-in will be applied in the grid search.
When a grid search with a zoom is applied, \code{zoom > 0} is actually the number of rounds or layers of the algorithm,
and therefore the grid search algorithm with a zoom consists of
\deqn{ n^{0} + n^{1} + n^{2} + ... + n^{z}   }
grid searches, where \eqn{n} = \code{num} and \eqn{z} = \code{zoom}.

As mentioned above, in each grid search, \code{num} is the number of points that will be returned.
And therefore, in the end, there will be
\deqn{ n^{1} + n^{2} + n^{3} + ... + n^{z+1}   }
points returned, where \eqn{n} = \code{num} and \eqn{z} = \code{zoom}.

Each time when the algorithm zooms in, it will automatically build subgrids
based on the points that have been found in the super gird search.
Due to the exhaustive property of the grid search algorithm,
it is desirable to make fewer points in the subgrid.
The decay rate \code{decay} provides the opportunity to control the number of points in the subgrids.
The number of points for each argument of the target function in the subgrid will be
\code{max} [ \code{Int} (\code{decay} \eqn{*} \eqn{N}), 3 ].

Parallel computation is implemented in the function, which can be activated by setting \code{parallel = TRUE}.

\code{cores}, which represents the number of cores, works only when \code{parallel = TRUE}.
By default \code{cores=NULL} implies that the function will detect the number of cores and use it.

The boolean \code{silent} controls if there will be output in the console.
}
\examples{
# Rastrigin function
ndim = 2 # number of dimension
nA = 10 # parameter A
# vx in [-5.12, 5.12]

# minimizer = rep(0, ndim)
# minimum = 0
Rastrigin <- function(vx) return(nA * ndim + sum(vx*vx - nA * cos(2*pi*vx)))

\donttest{
# set seed and initialize the initial or starting value
set.seed(1)
par = runif(ndim, -5.12, 5.12)
cat("start from", par)

# results from different optimization algorithms
optim(par = par, Rastrigin, method='Nelder-Mead')
optim(par = par, Rastrigin, method='BFGS')
optim(par = par, Rastrigin, method='L-BFGS-B')
optim(par = par, Rastrigin, method='SANN')
}

# a toy example
# build the grid first
bin = c(from=-5.12, to=5.12, by=.5)
grid = build_grid(bin,bin)
# so this is a relatively sparse grid

# serial computation
ret0 = grid_search(Rastrigin, grid, silent=FALSE)
ret0$par

\donttest{

# We can build a finer grid
bin = c(from=-5.12, to=5.12, by=.1)
grid = build_grid(bin,bin)

# serial computation
ret1 = grid_search(Rastrigin, grid, silent=FALSE)
ret1$par

# parallel computation
ret2 = grid_search(Rastrigin, grid, num=2, parallel=TRUE, silent=FALSE)
ret2$par

# grid search with a zoom!
ret3 = grid_search(Rastrigin, grid, zoom=2, num=2, parallel=TRUE, silent=FALSE)
ret3$par
}

}
\seealso{
\code{\link{build_grid}}, \code{\link{grid_search_check}}
}
\author{
Yukai Yang, \email{yukai.yang@statistik.uu.se}
}
\keyword{algorithms}
