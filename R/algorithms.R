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
#' The target function \code{FUN} to be minimized is a scalar real-valued function.
#' Maximization can be achieved by multiplying \code{-1} to the original function
#' and then passing the new function to \code{FUN}.
#'
#' The \code{grid} must be created by \code{\link{build_grid}}.
#'
#' Any other invariant arguments to \code{FUN} can be specified in \code{MoreArgs}
#' using a named list, see \code{\link{mapply}}.
#'
#' The common grid search first builds a grid within a bounded region, evaluates \code{FUN}
#' at each grid point, and returns the \code{num} points that yield the smallest values.
#'
#' \code{zoom = 0} implies no zoom-in (a single grid search). Any integer \code{zoom > 0}
#' applies additional zoom-in layers. With \code{zoom > 0}, the algorithm performs
#' \deqn{ n^{0} + n^{1} + n^{2} + \cdots + n^{z} }
#' grid searches, where \eqn{n} is \code{num} and \eqn{z} is \code{zoom}.
#' Consequently, the total number of returned points is
#' \deqn{ n^{1} + n^{2} + n^{3} + \cdots + n^{z+1} }.
#'
#' At each zoom-in layer, the algorithm builds subgrids around the best points found
#' in the previous layer. To limit the computational burden, the subgrid size is reduced
#' by the decay rate \code{decay}. For each parameter, the number of points in the subgrid is
#' \code{max(Int(decay * N), 3)}, where \code{N} is the number of points in the original grid
#' for that parameter.
#'
#' Parallel computation can be enabled by setting \code{parallel = TRUE}. In that case,
#' the function uses the \pkg{future} framework with \code{future::multisession}
#' (cross-platform). The number of workers is determined as follows:
#' \enumerate{
#'   \item Let \eqn{n} be the value returned by \code{future::availableCores()}.
#'   \item Let \eqn{m} be the user input \code{cores}. If \code{cores = NULL}, set \eqn{m = 2}.
#'   \item The number of workers is \eqn{\min(m, n)}.
#' }
#' If \code{parallel = TRUE}, the packages \pkg{future} and \pkg{future.apply} must be installed.
#'
#' The boolean \code{silent} controls whether progress information is printed to the console.
#'
#' @param FUN the target function to be minimized.
#' @param grid an object of class \code{GRID} created by \code{\link{build_grid}}.
#' @param MoreArgs a named list of additional arguments to \code{FUN}, see \code{\link{mapply}}.
#' @param zoom number of (additional) zoom-in layers, \code{0} by default.
#' @param decay a number in \eqn{(0,1)} controlling the decay of subgrid sizes.
#' @param num number of points to return at each grid search, \code{1} by default.
#' @param parallel a logical; if \code{TRUE}, parallel computation is used.
#' @param cores an integer specifying the requested number of workers when \code{parallel = TRUE}.
#'   If \code{NULL}, the function uses \code{2} workers by default (subject to \code{future::availableCores()}).
#' @param silent a logical indicating whether progress information is printed.
#'
#' @return a list with components:
#' \item{par}{the approximate global minimizer}
#' \item{points}{all candidate points found by the grid search with zoom-in layers}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso \code{\link{build_grid}}, \code{\link{grid_search_check}}
#' @keywords algorithms
#'
#' @examples
#' # Rastrigin function
#' ndim = 2
#' nA = 10
#' Rastrigin <- function(vx) nA * ndim + sum(vx * vx - nA * cos(2 * pi * vx))
#'
#' # Build a grid
#' bin = c(from = -5.12, to = 5.12, by = .5)
#' grid = build_grid(bin, bin)
#'
#' # Serial computation
#' ret0 = grid_search(Rastrigin, grid, silent = FALSE)
#' ret0$par
#'
#' \donttest{
#' # Finer grid
#' bin = c(from = -5.12, to = 5.12, by = .1)
#' grid = build_grid(bin, bin)
#'
#' # Serial computation
#' ret1 = grid_search(Rastrigin, grid, silent = FALSE)
#' ret1$par
#'
#' # Parallel computation (requires future and future.apply)
#' ret2 = grid_search(Rastrigin, grid, num = 2, parallel = TRUE, cores = 2, silent = FALSE)
#' ret2$par
#'
#' # Grid search with zoom-in layers
#' ret3 = grid_search(Rastrigin, grid, zoom = 2, num = 2, parallel = TRUE, cores = 2, silent = FALSE)
#' ret3$par
#' }
#'
#' @export
grid_search <- function(FUN, grid, MoreArgs=NULL, zoom=0, decay=0.5, num=1,
                        parallel=FALSE, cores=NULL, silent=TRUE){
  
  if (!inherits(grid, "GRID"))
    stop("The argument 'grid' is not an object of class 'GRID'.")
  
  if (!silent) {
    cat0(paste0(rep("#", getOption("width")), collapse = ""))
    cat0("zoomgrid version ", vnum, " ", packname)
    cat0(paste0(rep("-", getOption("width")), collapse = ""))
  }
  
  # Decide evaluation function and (optionally) configure futures once at the top level
  if (parallel) {
    
    # Ensure optional parallel backend is available (future is in Suggests)
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Parallel execution requires packages 'future' and 'future.apply'.")
    }
    
    # Detect available cores once via future
    n <- future::availableCores()
    
    # User-requested cores: if NULL, default to 2
    m <- cores
    if (is.null(m)) m <- 2L
    
    # Validate m
    if (!is.numeric(m) || length(m) != 1 || is.na(m) || m < 1) {
      m <- 1L
    } else {
      m <- as.integer(m)
    }
    
    # Choose workers as min(m, n)
    cores <- min(m, as.integer(n))
    if (is.na(cores) || cores < 1L) cores <- 1L
    
    # Set multisession workers once and restore user's plan on exit
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)
    
    ftmp <- grid_peval
    if (!silent) cat0("Parallel computation runs with ", cores, " workers.")
  } else {
    
    ftmp <- grid_seval
  }
  
  if (!silent) ptm <- proc.time()
  
  ret <- recursive_search(
    ftmp = ftmp, FUN = FUN, grid = grid, MoreArgs = MoreArgs,
    zoom = zoom, decay = decay, num = num, cores = cores
  )
  
  if (!silent) cat0("The Grid Search of ", zoom, " zoom-in layers with ", num,
                    " points each gives ", length(ret), " results.")
  
  par <- ftmp(FUN = FUN, grid = ret, MoreArgs = MoreArgs, num = 1, cores = cores)[[1]]
  
  if (!silent) cat0("The minimizer is believed to be in the neighbourhood of ", par, ".")
  
  if (!silent) {
    cat0(paste0(rep("-", getOption("width")), collapse = ""))
    print(proc.time() - ptm)
    cat0(paste0(rep("#", getOption("width")), collapse = ""))
  }
  
  return(list(par = par, points = ret))
}


#' Check the time consumed by running the grid search algorithm with a zoom.
#'
#' This function provides a quick runtime estimate for \code{\link{grid_search}} under the same settings.
#' It performs two short pilot runs on smaller grids (with \code{zoom = 0}) and extrapolates the expected
#' time for the full grid and the requested number of zoom-in layers.
#'
#' This is useful before launching a large run, for example on a compute server or under a batch system
#' such as SLURM, where an approximate runtime is needed to request resources.
#'
#' The boolean \code{silent} controls whether progress information is printed to the console.
#' For details on the algorithm and the meaning of the arguments, see \code{\link{grid_search}}.
#'
#' @param FUN the target function to be minimized.
#' @param grid an object of class \code{GRID} created by \code{\link{build_grid}}.
#' @param MoreArgs a named list of additional arguments to \code{FUN}, see \code{\link{mapply}}.
#' @param zoom number of (additional) zoom-in layers, \code{0} by default.
#' @param decay a number in \eqn{(0,1)} controlling the decay of subgrid sizes.
#' @param num number of points to return at each grid search, \code{1} by default.
#' @param parallel a logical; if \code{TRUE}, parallel computation is used.
#' @param cores an integer specifying the requested number of workers when \code{parallel = TRUE}.
#'   If \code{NULL}, the function uses \code{2} workers by default (subject to \code{future::availableCores()}).
#'   The number of workers used is \code{min(cores, future::availableCores())}.
#' @param silent a logical indicating whether progress information is printed.
#'
#' @return a numeric value giving the estimated runtime in seconds.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso \code{\link{build_grid}}, \code{\link{grid_search}}
#' @keywords algorithms
#'
#' @examples
#' # Rastrigin function
#' ndim <- 2
#' nA <- 10
#' Rastrigin <- function(vx) nA * ndim + sum(vx * vx - nA * cos(2 * pi * vx))
#'
#' # Build a grid
#' bin <- c(from = -5.12, to = 5.12, by = .5)
#' grid <- build_grid(bin, bin)
#'
#' # Estimate runtime (serial)
#' t_est <- grid_search_check(Rastrigin, grid, silent = FALSE)
#' t_est
#'
#' \donttest{
#' # Finer grid
#' bin <- c(from = -5.12, to = 5.12, by = .1)
#' grid <- build_grid(bin, bin)
#'
#' # Estimate runtime, then run the search
#' t_est <- grid_search_check(Rastrigin, grid, parallel = TRUE, cores = 2, silent = FALSE)
#' ret   <- grid_search(Rastrigin, grid, parallel = TRUE, cores = 2, silent = FALSE)
#' }
#'
#' @export
grid_search_check <- function(FUN, grid, MoreArgs=NULL, zoom=0, decay=0.5, num=1,
                              parallel=FALSE, cores=NULL, silent=TRUE){
  
  if (!inherits(grid, "GRID"))
    stop("The argument 'grid' is not an object of class 'GRID'.")
  
  if (!silent) {
    cat0(paste0(rep("#", getOption("width")), collapse = ""))
    cat0("zoomgrid version ", vnum, " ", packname)
    cat0(paste0(rep("-", getOption("width")), collapse = ""))
  }
  
  # Decide evaluation function and (optionally) configure futures once at the top level
  if (parallel) {
    
    # Ensure optional parallel backend is available (future is in Suggests)
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Parallel execution requires packages 'future' and 'future.apply'.")
    }
    
    # Detect available cores once via future
    n <- future::availableCores()
    
    # User-requested cores: if NULL, default to 2
    m <- cores
    if (is.null(m)) m <- 2L
    
    # Validate m
    if (!is.numeric(m) || length(m) != 1 || is.na(m) || m < 1) {
      m <- 1L
    } else {
      m <- as.integer(m)
    }
    
    # Choose workers as min(m, n)
    cores <- min(m, as.integer(n))
    if (is.na(cores) || cores < 1L) cores <- 1L
    
    # Set multisession workers once and restore user's plan on exit
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)
    
    ftmp <- grid_peval
    if (!silent) cat0("Parallel computation runs with ", cores, " workers.")
  } else {
    
    ftmp <- grid_seval
    
    # In serial mode, treat "cores" as 1 for the timing model
    cores <- 1L
  }
  
  set.seed(1)
  
  # test 1
  np <- 2000^(1 / grid$npar)
  tmp <- NULL
  for (base in grid$grid_base) {
    tmp <- paste(tmp,
                 paste0("c(", paste(sort(sample(base, min(np, length(base)))), collapse = ","), ")"),
                 sep = ", ")
  }
  tmp <- substr(tmp, 3, nchar(tmp))
  new_grid <- eval(parse(text = paste0("build_grid(", tmp, ")")))
  
  ptm <- proc.time()
  ret <- recursive_search(ftmp = ftmp, FUN = FUN, grid = new_grid, MoreArgs = MoreArgs,
                          zoom = 0, num = num, cores = cores)
  ptm <- proc.time() - ptm
  
  tx1 <- ptm[3]
  nn1 <- new_grid$size / cores
  
  # test 2
  np <- 3000^(1 / grid$npar)
  tmp <- NULL
  for (base in grid$grid_base) {
    tmp <- paste(tmp,
                 paste0("c(", paste(sort(sample(base, min(np, length(base)))), collapse = ","), ")"),
                 sep = ", ")
  }
  tmp <- substr(tmp, 3, nchar(tmp))
  new_grid <- eval(parse(text = paste0("build_grid(", tmp, ")")))
  
  ptm <- proc.time()
  ret <- recursive_search(ftmp = ftmp, FUN = FUN, grid = new_grid, MoreArgs = MoreArgs,
                          zoom = 0, num = num, cores = cores)
  ptm <- proc.time() - ptm
  
  tx2 <- ptm[3]
  nn2 <- new_grid$size / cores
  
  # forecast the time
  tmp <- solve(matrix(c(nn1, nn2, 1, 1), 2, 2), c(tx1, tx2))
  rr <- tmp[1]
  aa <- tmp[2]
  nn <- grid$size / cores
  
  tmp <- 0:zoom
  ret <- sum((nn * (decay^(tmp * grid$npar)) * rr + aa) * (num^tmp))
  
  if (!silent) {
    cat0("The expected time consumed by running the grid search is around ", ret, " seconds.")
    cat0(paste0(rep("#", getOption("width")), collapse = ""))
  }
  
  return(ret)
}
