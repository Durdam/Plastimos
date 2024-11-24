# function splitmultimask ##############################################################################################

# multimask = cmask.multi[[1]]; mode = "points";

splitMultiMask = function(multimask, mode = "points") {

  # input check
  if(!is.numeric(multimask)) stop("splitMultiMask implemented for numeric types only")
  if(!is.array(multimask)) stop("splitMultiMask implemented for matrices/arrays only")
  ndim = length(dim(multimask))
  if(is.na(ndim) | (length(ndim) != 1)) stop("ndim = NA or not of length 1")
  if(!(mode %in% c("points","mask"))) stop("mode not in points, mask")
  dims = dim(multimask)

  # split multimask
  unim = sort(setdiff(unique(c(multimask)),0))
  if(length(unim) < 1) return(list())
  splitmask = setNames(lapply(unim, function(u) {
    ixu = which(multimask == u)
    ixuInd = arrayInd(ixu, dims)
    bbox = sapply(1:ndim,function(k) setNames(range(ixuInd[,k]),c("min","max")))
    colnames(bbox) = paste0("D",1:ndim)
    if(mode == "points") {
      mask = ixu
    } else if(mode == "mask") {
      mask = multimask
      mask[] = 0
      mask[ixu] = 1
    }
    return(list("bbox" = bbox, "mask" = mask, "dim" = dims, "mode" = mode))
  }),unim)

  # return
  return(splitmask)
}
