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
  tmp = parallel::mcmapply(FUN=FUN, grid, MoreArgs=MoreArgs, mc.cores=cores)

  mins = list(); length(mins) = num
  tmax = max(tmp)
  for(nter in seq_len(num)){
    index = which.min(tmp)
    mins[[nter]] = as.numeric(grid[[index]])
    tmp[index] = tmax
  }

  return(mins)
}


build_subgrids <- function(mins, grid_base, decay=.5){
  ret = list()
  nbin = pmax(sapply(grid_base, length)*decay, 3)

  for(jter in seq_along(mins)){
    tmp = list()
    tmp_base = list(); length(tmp_base) = length(grid_base)

    for(iter in seq_along(grid_base)){
      index = which(mins[[jter]][iter] == grid_base[[iter]])
      tmp_base[[iter]] = seq(from=grid_base[[iter]][index-1],to=grid_base[[iter]][index+1],
                             length.out=nbin[iter]+2)[2:(nbin[iter]+1)]
    }

    tmp$grid = expand.grid(tmp_base,KEEP.OUT.ATTRS=FALSE)
    tmp$grid_base = tmp_base
    tmp$size = nrow(tmp$grid)
    tmp$npar = length(tmp_base)
    class(tmp) = "GRID"
    ret[[jter]] = tmp
  }

  return(ret)
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
