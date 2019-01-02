#################################################################################
## functions implementing the algorithms
#################################################################################

#' Build the grid for the grid search algorithm with a zoom.
#'
#' This function builds the grid for the grid search algorithm with a zoom.
#'
#' The argument \code{...} is a sequence of vectors or lists contatining the information about the grid to be built.
#' Each element in the sequence is either a vector or a list taking one of the following forms
#'
#' - \code{x}, if \code{x} is already an sequence of the grid points for the corresponding argument.
#'
#' - \code{c(from=, to=, by=)}
#'
#' - \code{c(from=, to=, length=)}
#'
#' - \code{list(from=, to=, by=)}
#'
#' - \code{list(from=, to=, length=)}
#'
#' where
#'
#' - \code{from}: the min of the argument of the target function
#'
#' - \code{to}: the max of the argument of the target function
#'
#' - \code{by}: the increment of the sequence
#'
#' - \code{length}: desired length.
#'
#'
#' There are many different ways to organize the points on the grid for certain argument of the target function,
#' the user can make them freely and input directly by \code{build_grid(x, ...)}.
#' Notice that \code{x} does not need to be increasing, as the function will sort it.
#' The design that \code{x} does not need to be increasing makes it convenient for the user
#' to interpolate more points at some region without considering to sort it all the time.
#'
#' When \code{by} is provided, the \code{length} will be ignored.
#' So if the user wanna specify the \code{length}, please do not use \code{by}.
#'
#' The order of the sequence \code{...} matters as it represents the order of the corresponding arguments of the target function to be optimized.
#'
#' @param ... a sequence of vectors or lists containing the information about the grid to be built, see Usage and Details.
#'
#' @return a new object of the class GRID with the grid ready for the grid search with a zoom.
#'
#' The object contains the following components:
#' \item{grid}{the grid}
#' \item{size}{number of points in the grid}
#' \item{npar}{number of arguments or parameters}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso \code{\link{grid_search_check}}, \code{\link{grid_search}}
#' @keywords algorithms
#'
#' @examples
#' vx = 1:5
#' build_grid(vx, c(from=1, to=2, by=.2), list(from=3, to=4, length=5))
#'
#' @export
build_grid <- function(...){
  args = list(...)
  ret = list()

  grid_base = list(); iter = 1
  for(arg in args){
    if(is.null(arg)) next

    if(is.list(arg)){
      nfrom = arg$from; if(is.null(nfrom)) next
      nto = arg$to; if(is.null(nto)) next
      if(nfrom > nto)
        stop(simpleError("'from' is bigger than 'to'! Please check what you have input..."))

      nby = arg$by
      if(!is.null(nby)){
        if(nby <= 0) stop(simpleError("'by' is not positive..."))

        # with from, to, by
        tmp = seq(from=nfrom, to=nto, by=nby)
      }else{
        nlength = arg$length; if(is.null(nlength)) next
        if(nlength <= 0) stop(simpleError("'length' is not positive..."))

        # with from, to, length
        tmp = seq(from=nfrom, to=nto, length.out=nlength)
      }
    }else if(is.vector(arg)){
      if(!is.numeric(arg)) next

      nfrom = arg["from"]
      if(is.na(nfrom)){
        # x case
        tmp = sort(arg)
      }else{
        nto = arg["to"]; if(is.na(nto)) next
        if(nfrom > nto)
          stop(simpleError("'from' is bigger than 'to'! Please check what you have input..."))

        nby = arg["by"]
        if(!is.na(nby)){
          if(nby <= 0) stop(simpleError("'by' is not positive..."))

          # with from, to, by
          tmp = seq(from=nfrom, to=nto, by=nby)
        }else{
          nlength = arg["length"]; if(is.na(nlength)) next
          if(nlength <= 0) stop(simpleError("'length' is not positive..."))

          # with from, to, length
          tmp = seq(from=nfrom, to=nto, length.out=nlength)
        }
      }
    }else next

    grid_base[[iter]] = tmp; iter = iter+1
  }

  ret$grid = expand.grid(grid_base, KEEP.OUT.ATTRS=FALSE)
  ret$grid_base = grid_base
  ret$size = nrow(ret$grid)
  ret$npar = length(grid_base)
  class(ret) = "GRID"
  return(ret)
}


#' Carry out the grid search algorithm with a zoom.
#'
#' This function carries out the grid search algorithm with a zoom.
#'
#' The target function \code{FUN} to be minimized is a scalar real valued function with multiple arguments.
#' The maximization can be achieved by multiplying -1 to the original function,
#' and then input the new function to \code{FUN}.
#'
#' The \code{grid} must be created by using the function \code{\link{build_grid}} function in the package.
#'
#' Any other invariant arugments to the function \code{FUN} can be specified in \code{MoreArgs} by using a list with variable names.
#'
#' The common grid search first build a grid within some bounded area.
#' And then the target function will be evaluated at each point in the grid.
#' The points that produce the smallest \code{num} target function values shall be returned.
#'
#' \code{zoom = 0} implies that no zoom-in will be applied and the corresponding grid search is the common one introduced above,
#' while any integer \code{zoom > 0} implies that the zoom-in will be applied in the grid search.
#' When a grid search with a zoom is applied, \code{zoom > 0} is actually the number of rounds or layers of the algorithm,
#' and therefore the grid search algorithm with a zoom consists of
#' \deqn{ n^{0} + n^{1} + n^{2} + ... + n^{z}   }
#' grid searches, where \eqn{n} = \code{num} and \eqn{z} = \code{zoom}.
#'
#' As mentioned above, in each grid search, \code{num} is the number of points that will be returned.
#' And therefore, in the end, there will be
#' \deqn{ n^{1} + n^{2} + n^{3} + ... + n^{z+1}   }
#' points returned, where \eqn{n} = \code{num} and \eqn{z} = \code{zoom}.
#'
#' Each time when the algorithm zooms in, it will automatically build subgrids
#' based on the points that have been found in the super gird search.
#' Due to the exhaustive property of the grid search algorithm,
#' it is desirable to make fewer points in the subgrid.
#' The decay rate \code{decay} provides the opportunity to control the number of points in the subgrids.
#' The number of points for each argument of the target function in the subgrid will be
#' \code{max} [ \code{Int} (\code{decay} \eqn{*} \eqn{N}), 3 ].
#'
#' Parallel computation is implemented in the function, which can be activated by setting \code{parallel = TRUE}.
#'
#' \code{cores}, which represents the number of cores, works only when \code{parallel = TRUE}.
#' By default \code{cores=NULL} implies that the function will detect the number of cores and use it.
#'
#' The boolean \code{silent} controls if there will be output in the console.
#'
#'
#' @param FUN the target function to be minimized.
#' @param grid an object of the class GRID from \code{\link{build_grid}}.
#' @param MoreArgs a list of other arguments to \code{FUN}, see \code{\link{mapply}}.
#' @param zoom number of (additional) rounds or layers of the zoom-in, 0 by default.
#' @param decay a number in between 0 and 1 representing the decay rate of the grid sizes of the zoom.
#' @param num number of points to return, i.e. the smallest \code{num} points, 1 by default the minimum.
#' @param parallel a boolean indicating if the parallel computation is carried out, by default \code{FALSE}.
#' @param cores The number of cores to use, i.e. at most how many child processes will be run simultaneously. For details, see \code{mcmapply} in \code{parallel} package.
#' @param silent a boolean indicating if the information regarding the omputation is printed.
#'
#' @return a list containing the results from the grid search with a zoom.
#'
#' The list contains the following components:
#' \item{par}{the approximate global minimizer}
#' \item{points}{all the local minimizer points found by the grid search with a zoom}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso \code{\link{build_grid}}, \code{\link{grid_search_check}}
#' @keywords algorithms
#'
#' @examples
#' # Rastrigin function
#' ndim = 2 # number of dimension
#' nA = 10 # parameter A
#' # vx in [-5.12, 5.12]
#'
#' # minimizer = rep(0, ndim)
#' # minimum = 0
#' Rastrigin <- function(vx) return(nA * ndim + sum(vx*vx - nA * cos(2*pi*vx)))
#'
#' \donttest{
#' # set seed and initialize the initial or starting value
#' set.seed(1)
#' par = runif(ndim, -5.12, 5.12)
#' cat("start from", par)
#'
#' # results from different optimization algorithms
#' optim(par = par, Rastrigin, method='Nelder-Mead')
#' optim(par = par, Rastrigin, method='BFGS')
#' optim(par = par, Rastrigin, method='L-BFGS-B')
#' optim(par = par, Rastrigin, method='SANN')
#' }
#'
#' # a toy example
#' # build the grid first
#' bin = c(from=-5.12, to=5.12, by=.5)
#' grid = build_grid(bin,bin)
#' # so this is a relatively sparse grid
#'
#' # serial computation
#' ret0 = grid_search(Rastrigin, grid, silent=FALSE)
#' ret0$par
#'
#' \donttest{
#'
#' # We can build a finer grid
#' bin = c(from=-5.12, to=5.12, by=.1)
#' grid = build_grid(bin,bin)
#'
#' # serial computation
#' ret1 = grid_search(Rastrigin, grid, silent=FALSE)
#' ret1$par
#'
#' # parallel computation
#' ret2 = grid_search(Rastrigin, grid, num=2, parallel=TRUE, silent=FALSE)
#' ret2$par
#'
#' # grid search with a zoom!
#' ret3 = grid_search(Rastrigin, grid, zoom=2, num=2, parallel=TRUE, silent=FALSE)
#' ret3$par
#' }
#'
#' @export
grid_search <- function(FUN, grid, MoreArgs=NULL, zoom=0, decay=0.5, num=1, parallel=FALSE, cores=NULL, silent=TRUE){
  if(class(grid)!="GRID")
    stop(simpleError("The argument 'grid' is not an object of class 'GRID'."))

  if(!silent){
    cat0(paste0(rep("#",getOption("width")),collapse=''))
    cat0("zoomgrid version ",vnum," ",packname)
    cat0(paste0(rep("-",getOption("width")),collapse=''))
  }

  # tmp = split(grid$grid, seq_len(grid$size))

  if(is.null(cores)){
    cores = parallel::detectCores()
    if(is.na(cores)){
      if(!silent) cat0("No cores are detected! Anyway, let's use one core...")
      parallel=FALSE
    }
  }

  if(parallel){
    ftmp = grid_peval
    if(!silent) cat0("Parallel computation runs with ",cores," cores.")
  }else ftmp = grid_seval

  if(!silent) ptm = proc.time()

  ret = recursive_search(ftmp=ftmp,FUN=FUN,grid=grid,MoreArgs=MoreArgs,zoom=zoom,decay=decay,num=num,cores=cores)

  if(!silent) cat0("The Grid Search of ",zoom," zoom-in layers with ",num," points each gives ",length(ret), " results.")
  par = ftmp(FUN=FUN,grid=ret,MoreArgs=MoreArgs,num=1,cores=cores)[[1]]
  if(!silent) cat0("The minimizer is believed to be in the neighbourhood of ",par,".")

  if(!silent){
    cat0(paste0(rep("-",getOption("width")),collapse=''))
    print(proc.time() - ptm)
    cat0(paste0(rep("#",getOption("width")),collapse=''))
  }

  return(list(par=par, points=ret))
}


#' Check the time consumed by running the grid search algorithm with a zoom.
#'
#' This function checks the time consumed by running the grid search algorithm with a zoom as well as some other conditions.
#'
#' The running of this function takes only several seconds.
#' So it is recommended to run this function before \code{\link{grid_search}} to check the approximate time consumed
#' by \code{\link{grid_search}} by using exactly the same arguments.
#'
#' This function is extremely useful when the user is going to run \code{\link{grid_search}} on some super-computing server
#' and need to know approximately how long time it will take in order to specify the corresponding settings
#' according to some batch system like SLURM for example.
#'
#' The boolean \code{silent} controls if there will be output in the console.
#'
#' For details, see \code{\link{grid_search}}.
#'
#' @param FUN the target function to be minimized.
#' @param grid an object of the class GRID from \code{\link{build_grid}}.
#' @param MoreArgs a list of other arguments to \code{FUN}, see \code{\link{mapply}}.
#' @param zoom number of (additional) rounds or layers of the zoom-in, 0 by default.
#' @param decay a number in between 0 and 1 representing the decay rate of the grid sizes of the zoom.
#' @param num number of points to return, i.e. the smallest \code{num} points, 1 by default the minimum.
#' @param parallel a boolean indicating if the parallel computation is carried out, by default \code{FALSE}.
#' @param cores The number of cores to use, i.e. at most how many child processes will be run simultaneously. For details, see \code{mcmapply} in \code{parallel} package.
#' @param silent a boolean indicating if the information regarding the computation is printed.
#'
#' @return a number of the time in seconds.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso \code{\link{build_grid}}, \code{\link{grid_search}}
#' @keywords algorithms
#'
#' @examples
#' # Rastrigin function
#' ndim = 2 # number of dimension
#' nA = 10 # parameter A
#' # vx in [-5.12, 5.12]
#'
#' # minimizer = rep(0, ndim)
#' # minimum = 0
#' Rastrigin <- function(vx) return(nA * ndim + sum(vx*vx - nA * cos(2*pi*vx)))
#'
#' # a toy example
#' # build the grid first
#' bin = c(from=-5.12, to=5.12, by=.5)
#' grid = build_grid(bin,bin)
#' # so this is a relatively sparse grid
#'
#' # serial computation
#' ret0 = grid_search(Rastrigin, grid, silent=FALSE)
#' ret0$par
#'
#' \donttest{
#' # If we expand the grid to allow for more points
#' bin = c(from=-5.12, to=5.12, by=.1)
#' grid = build_grid(bin,bin)
#'
#' # run the check before the grid search
#' ret1 = grid_search_check(Rastrigin, grid, silent=FALSE)
#' ret1 = grid_search(Rastrigin, grid, silent=FALSE)
#' }
#'
#' @export
grid_search_check <- function(FUN, grid, MoreArgs=NULL, zoom=0, decay=0.5, num=1, parallel=FALSE, cores=NULL, silent=TRUE){
  if(class(grid)!="GRID")
    stop(simpleError("The argument 'grid' is not an object of class 'GRID'."))

  if(!silent){
    cat0(paste0(rep("#",getOption("width")),collapse=''))
    cat0("zoomgrid version ",vnum," ",packname)
    cat0(paste0(rep("-",getOption("width")),collapse=''))
  }

  if(is.null(cores)){
    cores = parallel::detectCores()
    if(is.na(cores)){
      if(!silent) cat0("No cores are detected! Anyway, let's use one core...")
      parallel=FALSE
    }
  }

  if(parallel){
    ftmp = grid_peval
    if(!silent) cat0("Parallel computation runs with ",cores," cores.")
  }else ftmp = grid_seval

  set.seed(1)

  # test 1
  np = 2000**(1/grid$npar)
  tmp = NULL
  for(base in grid$grid_base){
    tmp = paste(tmp, paste0("c(",paste(sort(sample(base,min(np,length(base)))),collapse=','),")"), sep=", ")
  }
  tmp = substr(tmp, 3, nchar(tmp))
  new_grid = eval(parse(text=paste0("build_grid(",tmp,")")))

  ptm = proc.time()
  ret = recursive_search(ftmp=ftmp,FUN=FUN,grid=new_grid,MoreArgs=MoreArgs,zoom=0,num=num,cores=cores)
  ptm = proc.time() - ptm

  tx1 = ptm[3]; nn1 = new_grid$size/cores

  # test 2
  np = 3000**(1/grid$npar)
  tmp = NULL
  for(base in grid$grid_base){
    tmp = paste(tmp, paste0("c(",paste(sort(sample(base,min(np,length(base)))),collapse=','),")"), sep=", ")
  }
  tmp = substr(tmp, 3, nchar(tmp))
  new_grid = eval(parse(text=paste0("build_grid(",tmp,")")))

  ptm = proc.time()
  ret = recursive_search(ftmp=ftmp,FUN=FUN,grid=new_grid,MoreArgs=MoreArgs,zoom=0,num=num,cores=cores)
  ptm = proc.time() - ptm

  tx2 = ptm[3]; nn2 = new_grid$size/cores

  # forecast the time

  tmp = solve(matrix(c(nn1,nn2,1,1),2,2), c(tx1,tx2))
  rr = tmp[1]; aa = tmp[2]
  nn = grid$size/cores

  tmp = 0:zoom
  ret = sum( (nn*(decay**(tmp*grid$npar)) * rr + aa) * (num**tmp) )
  if(!silent){
    cat0("The expected time consumed by running the grid search is around ",ret," seconds.")
    cat0(paste0(rep("#",getOption("width")),collapse=''))
  }

  return(ret)
}
