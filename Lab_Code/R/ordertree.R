#!/usr/bin/env Rgetopt.py
library(Rgetopt)

argspec <- c("ordertree.R - order a matrix's rows by hierarchical clustering then
switiching branch order to minimize distance between adjacent rows.
Assumes a datafile with a header line and a single column of rownames.

WARNING: this can be very slow for matrices with large number of rows (>400)!

Usage:
    ordertree.R [options] matrix_file

Options:",
             "dist=s      distance measure from either cor() or dist()",
             "abs         calculate distances using absolute values of data",
             "fillMissing replace missing distances with 2X the max. distance",
             "columns     order the columns rather than rows",
             "quick|q     just hierarchically cluster, don't order")
main <- function() {
  o <- Rgetopt(argspec, defaults=list(dist="euclidean"))

  if(length(o$argv) != 1) usage("Need an input matrix", argspec=argspec)
  df <- readTabFile(o$argv)

  d <- as.matrix(df[,2:ncol(df)])
  if (!is.null(o$columns)) d <- t(d)
  if (!is.null(o$abs)) d <- abs(d)

  if (o$dist %in% c("pearson", "kendall", "spearman")) {
    dist.d <- as.dist(1 - cor(t(d), use="pairwise.complete", method=o$dist))
  } else {
    dist.d <- dist(d, method=o$dist)
  }

  if (any(is.na(dist.d))) {
    if (is.null(o$fillMissing)) {
      stop("Some distances can not be calculated, consider --fillMissing")
    } else {
      m <- max(dist.d, na.rm=TRUE)
      dist.d[is.na(dist.d)] <- m
    }
  }

  hc <- hclust(dist.d)
  if(is.null(o$quick)) {
    reorder <- ordertree(dist.d, hc)
  } else {
    reorder <- hc$order
  }
  df <- if (is.null(o$columns)) df[reorder,] else df[,c(1,reorder+1)]
  
  write.table(df, file=stdout(), quote=F, sep="\t", row.names=FALSE)
}

###
### Outputs an optimal ordering for the leaves of an heirarchical clustering
###
ordertree <- function(distobj, hc) {
  d <- as.matrix(distobj)
  dmin <- min(d)
  if (dmin < 0) { d <- d - dmin + 1 }
  n <- length(hc$order)
  
  m <- matrix(nrow = n, ncol = n)  # the min matrix
  height <- matrix(nrow = n, ncol = n) # which v at m[i,j]
  bl <- matrix(nrow = n, ncol = n)  # the backtrace matrix
  br <- matrix(nrow = n, ncol = n)  # the backtrace matrix

  child <- list()
  children <- function(x) { if (x < 0) -x else child[[x]] }

  # initialization
#  for (i in 1:n) { m[i,i] = 0; height[i,i] = -i}
  diag(m) <- 0
  diag(height) <- -(1:nrow(height))

  # iteration
  for (v in 1:nrow(hc$merge)) {
#    if (v == 4) browser()

    left <- hc$merge[v,1]
    right <- hc$merge[v,2]

    cl <- children(left)
    cr <- children(right)
    temp <- matrix(nrow=length(cl), ncol=length(cr))
    t.min <- matrix(nrow=length(cl), ncol=length(cr))
    for (i in 1:length(cl)) {
      for (l in 1:length(cr)) {
        scores <- m[cl[i], cl] + d[cl, cr[l]]
        scores[height[cl[i],cl] != hc$merge[v,1]] <- Inf
        temp[i,l] <- min(scores)
        t.min[i,l] <- cl[which.min(scores)]
      }
    }
    if (length(cl) == 1) m[cl,cl] <- Inf
    for (i in 1:length(cl)) {
      for (j in 1:length(cr)) {
        scores <- temp[i, ] + m[cr,cr[j]]
        scores[height[cr[j],cr] != hc$merge[v,2]] <- Inf
        I <- cl[i]; J <- cr[j]
        m[I,J] <- m[J,I] <- min(scores)
        bl[I,J] <- br[J,I] <- t.min[i,which.min(scores)]
        br[I,J] <- bl[J,I] <- cr[which.min(scores)]
        height[I,J] <- height[J,I] <- v
      }
    }
    if (length(cr) == 1) m[cr,cr] <- Inf
    child[[v]] <- c(children(left), children(right))
  }
#browser()

  # traceback
  result <- vector(length=n)
  overall.min <- min(m[height==(n-1)])
  which.omin <- which(m==overall.min)[1]
  leftnode <- col(m)[which.omin]
  rightnode <- row(m)[which.omin]
  if (leftnode == 0) leftnode <- n

  traceback <- function(start, leftnode, rightnode) {
    p <- height[leftnode,rightnode]
    if(p < 0) return();
    if (leftnode %in% children(hc$merge[p,1])) {
      n <- length(children(hc$merge[p,1]))
    } else { n <- length(children(hc$merge[p,2])) }
    result[start+n-1] <<- bl[leftnode, rightnode]
    traceback(start, leftnode, bl[leftnode,rightnode])
    result[start+n] <<- br[leftnode,rightnode]
    traceback(start+n, br[leftnode,rightnode], rightnode)
  }

  traceback(1, leftnode, rightnode)

  result
}

if (Sys.getenv("RGETOPT_DEBUG") != "") {debug(main)}

main()
