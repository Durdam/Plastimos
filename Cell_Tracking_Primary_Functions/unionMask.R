# function unionMask ###################################################################################################

# mask1 = nmask; mask2 = cmask; mode="size"
# mask1 = cmask.mult[[1]]; mask2 = cmask.mult[[4]]; mode="size"; multimask = TRUE
# mask1 = cmask.mult[[1]]; mask2 = cmask.mult[[4]]; mode="mask"; multimask = TRUE

unionMask = function(mask1,mask2,mode="size",multimask = FALSE) {

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


  # size of union of two masks -----------------------------------------------------------------------------------------

  if(mode == "size") {

    m1m2 = matrix(0,nrow=n1,ncol=n2,dimnames=list(names(mask1),names(mask2)))
    for(i1 in 1:n1) {                                             # mask1 loop
      b1 = mask1[[i1]]$bbox
      for(i2 in 1:n2) {                                           # mask2 loop
        b2 = mask2[[i2]]$bbox
        # quick bounding box overlap check
        bmin = sapply(1:ndim,function(k) max(b1["min",k],b2["min",k]))
        bmax = sapply(1:ndim,function(k) min(b1["max",k],b2["max",k]))
        if(any(bmax < bmin)) {                                    # no overlap
          m1m2[i1,i2] = length(mask1[[i1]]$mask) + length(mask2[[i2]]$mask)
          next
        }

        # explicit point union
        m1m2[i1,i2] = length(union(mask1[[i1]]$mask,mask2[[i2]]$mask))
      }
    }

  }

  # union mask of two masks --------------------------------------------------------------------------------------------

  if(mode == "mask") {

    m0 = matrix(0,nrow=nrow(mask1[[1]]$mask),ncol=ncol(mask1[[1]]$mask))
    m1m2 = lapply(1:n1,function(i1) {                             # mask1 loop
      b1 = mask1[[i1]]$bbox
      m1 = mask1[[i1]]$mask
      ixbb1 = lapply(1:ndim,function(k) b1["min",k]:b1["max",k])
      ixbb1 =  as.matrix(do.call(expand.grid,ixbb1))
      lapply(1:n2,function(i2) {                                  # mask2 loop
        b2 = mask2[[i2]]$bbox
        m2 = mask2[[i2]]$mask

        # get index of bounding box union and compute union
        ixbb2 = lapply(1:ndim,function(k) b2["min",k]:b2["max",k])
        ixbb2 = as.matrix(do.call(expand.grid,ixbb2))
        ixbb = rbind(ixbb1,ixbb2)

        m1m2 = m0
        m1m2[ixbb] = ifelse((m1[ixbb]>0)|(m2[ixbb]>0),1,0)
        return(m1m2)
      })
    })

  }

  # return
  return(m1m2)
}
