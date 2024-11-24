# function interMask ###################################################################################################
# computes intersection between image masks (distinction by =0 and >0 only)
# mode="mask" -> intersecion mask, mode="size" -> size of intersection only (maybe faster)

# mask1 = nmask; mask2 = cmask; mode="size"
# mask1 = cmask.multi[[1]]; mask2 = cmask.multi[[4]]; mode="size"; multimask = TRUE
# mask1 = cmask.multi[[1]]; mask2 = cmask.multi[[4]]; mode="mask"; multimask = TRUE

interMask <- function(mask1,mask2,mode="size",multimask = FALSE) {

  # input check
  if(!(mode %in% c("size","mask"))) stop("mode not in size, mask")

  if(multimask) {
    if(mode == "size") {
      mask1 = splitMultiMask(mask1, mode = "points")
      mask2 = splitMultiMask(mask2, mode = "points")
    } else if(mode == "mask") {
      mask1 = splitMultiMask(mask1, mode = "mask")
      mask2 = splitMultiMask(mask2, mode = "mask")
    }
  }
  n1 = length(mask1)
  n2 = length(mask2)

  # check dimension
  dim = mask1[[1]]$dim
  ndim = length(dim)
  m1m2 = append(mask1,mask2)
  if(!all(sapply(m1m2,function(m) (length(m$dim) == ndim)))) stop("masks differ in dimension")
  if(!all(sapply(m1m2,function(m) all(m$dim == dim)))) stop("masks differ in dimension")

  # check mode
  if(mode == "size") {
    if(!all(sapply(m1m2,function(m) (m$mode == "points"))))  stop("not all mode = points")
  } else if(mode == "mask") {
    if(!all(sapply(m1m2,function(m) (m$mode == "mask"))))  stop("not all mode = mask")
  }

  # size of intersecion between two masks ------------------------------------------------------------------------------

  if(mode == "size") {

    m1m2 = matrix(0,nrow=n1,ncol=n2,dimnames=list(names(mask1),names(mask2)))
    for(i1 in 1:n1) {                                             # mask1 loop
      b1 = mask1[[i1]]$bbox
      for(i2 in 1:n2) {                                           # mask2 loop
        b2 = mask2[[i2]]$bbox
        # quick bounding box overlap check
        bmin = sapply(1:ndim,function(k) max(b1["min",k],b2["min",k]))
        bmax = sapply(1:ndim,function(k) min(b1["max",k],b2["max",k]))
        if(any(bmax < bmin)) next                                 # no overlap

        # explicit point intersection
        m1m2[i1,i2] = length(intersect(mask1[[i1]]$mask,mask2[[i2]]$mask))
      }
    }

  }

  # intersection mask between two masks --------------------------------------------------------------------------------

  if(mode == "mask") {

    m0 = array(0,dim=dim(mask1[[1]]$mask))
    m1m2 = lapply(1:n1,function(i1) {                             # mask1 loop
      b1 = mask1[[i1]]$bbox
      m1 = mask1[[i1]]$mask
      lapply(1:n2,function(i2) {                                  # mask2 loop
        b2 = mask2[[i2]]$bbox
        # quick bounding box overlap check
        bmin = sapply(1:ndim,function(k) max(b1["min",k],b2["min",k]))
        bmax = sapply(1:ndim,function(k) min(b1["max",k],b2["max",k]))
        if(any(bmax < bmin)) return(m0)                           # no overlap
        m2 = mask2[[i2]]$mask

        # get index of bounding box overlap and intersect
        ixbb = lapply(1:ndim,function(k) bmin[k]:bmax[k])
        ixbb = as.matrix(do.call(expand.grid,ixbb))
        m1m2 = m0
        m1m2[ixbb] = ifelse((m1[ixbb]>0)&(m2[ixbb]>0),1,0)
        return(m1m2)
      })
    })

  }

  # return
  return(m1m2)
}
