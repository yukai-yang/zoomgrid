#################################################################################
## utility functions
#################################################################################

vnum = "1.0.0"
packname = "(Red Grid)"

# simple cat
cat0 <- function(...)
{
  words = list(...)
  for(tmp in words) cat(tmp)
  cat("\n")
}


grid_seval <- function(FUN, grid, MoreArgs, num, ...){
  tmp = mapply(FUN=FUN, grid, MoreArgs=MoreArgs)

  mins = list(); length(mins) = num
  tmax = max(tmp)
  for(nter in seq_len(num)){
    index = which.min(tmp)
    mins[[nter]] = as.numeric(grid[[index]])
    tmp[index] = tmax
  }

  return(mins)
}


grid_peval <- function(FUN, grid, MoreArgs, num, cores){
  
  # Parallel evaluation over a list of grid points using the current future plan
  tmp <- future.apply::future_lapply(
    grid,
    function(par_vec) {
      if (is.null(MoreArgs)) {
        FUN(par_vec)
      } else {
        do.call(FUN, c(list(par_vec), MoreArgs))
      }
    },
    future.seed = TRUE
  )
  
  tmp <- unlist(tmp, use.names = FALSE)
  
  # Select the smallest 'num' values
  mins <- vector("list", num)
  tmax <- max(tmp)
  
  for (nter in seq_len(num)) {
    index <- which.min(tmp)
    mins[[nter]] <- as.numeric(grid[[index]])
    tmp[index] <- tmax
  }
  
  return(mins)
}


build_subgrids <- function(mins, grid_base, decay = 0.5) {
  ret <- vector("list", length(mins))
  
  # keep your behaviour: nbin >= 3
  nbin <- pmax(as.integer(ceiling(vapply(grid_base, length, integer(1)) * decay)), 3L)
  
  for (jter in seq_along(mins)) {
    tmp_base <- vector("list", length(grid_base))
    
    for (iter in seq_along(grid_base)) {
      gb <- grid_base[[iter]]
      
      # Contract assumptions (internal assertions, not user-facing validation)
      # 1) gb must have at least 3 points if we require an interior index
      if (length(gb) < 3L) {
        stop(sprintf("build_subgrids(): grid_base[[%d]] has length < 3; cannot form interior neighbourhood.", iter))
      }
      
      # 2) mins point must match exactly one grid point
      idx <- which(mins[[jter]][iter] == gb)
      if (length(idx) != 1L) {
        stop(sprintf(
          "build_subgrids(): mins[[%d]][%d] does not match exactly one point in grid_base[[%d]] (matches=%d).",
          jter, iter, iter, length(idx)
        ))
      }
      
      # 3) you said boundary case does not occur: assert strict interior
      if (idx == 1L || idx == length(gb)) {
        stop(sprintf(
          "build_subgrids(): mins[[%d]][%d] lies on boundary of grid_base[[%d]] (index=%d).",
          jter, iter, iter, idx
        ))
      }
      
      tmp_base[[iter]] <-
        seq(from = gb[idx - 1L], to = gb[idx + 1L], length.out = nbin[iter] + 2L)[2:(nbin[iter] + 1L)]
    }
    
    tmp <- list()
    tmp$grid <- expand.grid(tmp_base, KEEP.OUT.ATTRS = FALSE)
    tmp$grid_base <- tmp_base
    tmp$size <- nrow(tmp$grid)
    tmp$npar <- length(tmp_base)
    class(tmp) <- "GRID"
    
    ret[[jter]] <- tmp
  }
  
  ret
}


recursive_search <- function(ftmp, FUN, grid, MoreArgs, zoom, decay, num, cores){
  split_grid = split(grid$grid, seq_len(grid$size))
  ret = ftmp(FUN=FUN, grid=split_grid, MoreArgs=MoreArgs, num=num, cores=cores)

  if(zoom > 0){
    # build subgrids
    subgrids = build_subgrids(mins=ret,grid_base=grid$grid_base,decay=decay)

    zoom = zoom -1
    for(subgrid in subgrids){
      ret = c(ret,recursive_search(ftmp=ftmp,FUN=FUN,grid=subgrid,MoreArgs=MoreArgs,zoom=zoom,decay=decay,num=num,cores=cores))
    }
  }

  return(ret)
}
