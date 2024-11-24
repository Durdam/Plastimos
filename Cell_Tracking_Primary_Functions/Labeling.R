
# mask = cmask.multi; multimask = TRUE; ncoressys = 1;

Labeling <- function(mask, multimask = TRUE, ncoressys = 1) {

  # input check
  nframes = length(mask)
  if(nframes < 2) stop("no tracking with less than two frames")

  # set up parallel computing
  parallelOK = (requireNamespace("doParallel",quietly=TRUE) & requireNamespace("foreach",quietly=TRUE))
  if(parallelOK) {
    ncores = parallel::detectCores()
    ncorespar = max(1,min(ncores,ncores-ncoressys))
    if(ncorespar > 1) {
      registerDoParallel(ncorespar)
    } else {
      parallelOK = FALSE
    }
  }

  # edges of bipartite graph

  if(parallelOK) {   # parallel computing if available -----------------------------------------------------------------

    Edges_all = foreach(ifra = 1:(nframes-1)) %dopar% {

      # splitting of consecutive frames required if processing is not known to be sequential
      if(multimask) {
        split1 = splitMultiMask(mask[[ifra]])
        split2 = splitMultiMask(mask[[ifra+1]])
      } else {
        split1 = mask[[ifra]]
        split2 = mask[[ifra+1]]
      }

      # calculate edges
      return(.get_Edges(split1,split2))
    }

  } else {    # serial computing ---------------------------------------------------------------------------------------

    # save splitting time if processing is known to be sequential
    if(multimask) {
      split1 = splitMultiMask(mask[[1]])
    } else {
      split1 = mask[[1]]
    }

    # calculate edges
    Edges_all = vector('list',nframes-1)
    for(ifra in 1:(nframes-1)) {
      if(multimask) {
        split2 = splitMultiMask(mask[[ifra+1]])
      } else {
        split2 = mask[[ifra+1]]
      }

      # calculate edges
      Edges_all[[ifra]] = .get_Edges(split1,split2)

      # shuffle splits
      split1 = split2
    }

  } # parallel ---------------------------------------------------------------------------------------------------------

  # generate globally unique cell labels ###############################################################################

  Labels_all = setNames(vector('list',nframes),names(mask))

  # initialize labels
  Labels_all[[1]] = Edges_all[[1]][,c("index t","index t"),drop=FALSE]
  colnames(Labels_all[[1]]) = c("index t","label")
  ix = which(!duplicated(Labels_all[[1]][,"index t"]) & (Labels_all[[1]][,"index t"]>0))
  Labels_all[[1]] = Labels_all[[1]][ix,,drop=FALSE]
  labmax = max(Labels_all[[1]][,"label"])
  for(ifra in 2:nframes) {

    lab = Labels_all[[ifra-1]]
    trans = Edges_all[[ifra-1]][,c("index t","index t+1"),drop=FALSE]
    trans0 = (apply(trans,1,min) == 0)
    nocc = sapply(trans[,"index t"],function(t) length(which(trans[,"index t"]==t)))
    ix11 = which((nocc==1) & (!trans0))
    ix12 = which((nocc==2) & (!trans0))
    ix01 = which((trans[,"index t"] == 0) & (trans[,"index t+1"] > 0))
    ix10 = which((trans[,"index t+1"] == 0) & (trans[,"index t"] > 0))
    ix = c(ix11,ix12,ix01,ix10)
    if((length(ix) != nrow(trans)) | !setequal(ix,1:nrow(trans))) stop("set not equal")

    # initialize new labels
    Labels = matrix(0, nrow=0, ncol=2,dimnames=list(NULL,c("index t","label")))

    # transition 1-1
    nx11 = length(ix11)
    if(nx11 > 0) {
      ix11m = match(trans[ix11,"index t"],lab[,"index t"])
      if(any(is.na(ix11m))) stop("ix11 not matching")
      lab = cbind(trans[ix11,"index t+1"],lab[ix11m,"label"])
      Labels = rbind(Labels,lab)
    }

    # transition 1-2
    nx12 = length(ix12)
    if(nx12 > 0) {
      lab = cbind(trans[ix12,"index t+1"],labmax+(1:nx12))
      Labels = rbind(Labels,lab)
      labmax = labmax+nx12
    }

    # transition 0-1
    nx01 = length(ix01)
    if(nx01 > 0) {
      lab = cbind(trans[ix01,"index t+1"],labmax+(1:nx01))
      Labels = rbind(Labels,lab)
      labmax = labmax+nx01
    }

    rownames(Labels) = NULL
    Labels_all[[ifra]] = Labels[order(Labels[,"index t"]),,drop=FALSE]
  }

  # t in index t not needed anymore
  for(ifra in 1:nframes) colnames(Labels_all[[ifra]]) = gsub("index t","index",colnames(Labels_all[[ifra]]))

  return(Labels_all)
}

.get_Edges <- function(split1,split2) {
  Meps = .Machine$double.eps

  JaccInd = interMask(split1,split2)/pmax(unionMask(split1,split2),Meps)
  n1 = nrow(JaccInd)
  n2 = ncol(JaccInd)

  # cases ###########################################################

  # 1-1 (tracking isolated cell) and 1-2 (cell division) assignment according to max JaccInd

  Edges = which(JaccInd > Meps, arr.ind = TRUE)
  Edges = cbind(Edges,JaccInd[Edges])
  colnames(Edges) = c("index t","index t+1","JaccInd")
  Edges = Edges[order(Edges[,"JaccInd"],decreasing=TRUE),,drop=FALSE]
  nEdges = nrow(Edges)
  if(nEdges > 1) {
    eval = sapply(2:nEdges, function(iE) {
      nocc1 = sum(Edges[1:iE,1] == Edges[iE,1])
      nocc2 = sum(Edges[1:iE,2] == Edges[iE,2])
      ((nocc1 <= 2) & (nocc2 <= 1)) })                    # nocc1 = 2: cell division
    Edges = Edges[c(TRUE,eval),,drop=FALSE]
  }

  # 1-0 out of frame
  ixout = setdiff(1:n1,Edges[,"index t"])
  nxout = length(ixout)
  if(nxout > 0) {
    Edges = rbind(Edges,cbind(ixout,rep(0,nxout),rep(0,nxout)))
  }

  # 0-1 new in frame
  ixin = setdiff(1:n2,Edges[,"index t+1"])
  nxin = length(ixin)
  if(nxin > 0) {
    Edges = rbind(Edges,cbind(rep(0,nxin),ixin,rep(0,nxin)))
  }

  rownames(Edges) = NULL
  Edges = Edges[order(Edges[,"index t"]),,drop=FALSE]

  return(Edges)
}

