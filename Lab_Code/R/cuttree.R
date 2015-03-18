#!/usr/bin/env Rgetopt.py
library(Rgetopt)

argspec <- c("cuttree.R - cut a hierachical clustering from Eisen Cluster
Be sure to specify .cdt on the end of the clustering name, and make sure
that both the .cdt and .gtr files are there.

WARNING: this program is pretty hackish for a .R script, and that's saying a
lot.  Be sure to sanity check your results.

Usage:
    cuttree.R [-h height|-k number_clusters] clustering.gtr",
             "k=i   number of clusters",
             "h=f   tree height to cut")

main <- function(argv) {
  if (missing(argv)) argv <- RgetArgvEnvironment()[-1]
  o <- Rgetopt(argv=argv, argspec=argspec)
  
  if (length(o$argv) == 0) {
    usage("must specify a file", argspec=argspec)
  }
  if (length(o$argv) > 1) {
    stop("Too many input files specified, must only have one")
  }
  infile <- o$argv[1]

  if (is.null(o$k) && is.null(o$h)) {
    stop("neither -h or -k set")
  }

  dump <- function(x,f=stdout()) {
    write.table(x,file=f,sep="\t",col.names=FALSE,quote=FALSE)
  }

  treefile <- paste(substr(infile, 1, nchar(infile)-4), ".gtr", sep='')

  hier <- read.cdt.tree(cdt.file=infile, tree.file=treefile)
  
  if (is.null(o$h)) {
    dump(cutree(hier, k=o$k))
  } else {
    dump(cutree(hier, h=o$h))
  }
}

#
# Need three attributes: merge, height, and labels
#  This should hopefully be easy....
#
# What still needs to be done when I have more time:  read in cdt data,
#  and also array tree

read.cdt.tree <- function(tree.file, cdt.file, cdt.headerlines=3){
  tree.data <- read.delim(tree.file, header=FALSE)
  
  merge.labels <- tree.data[,1]
  
  gene.labels <- system(paste("cut -f1", cdt.file), intern=TRUE)
  gene.labels <- gene.labels[-(1:cdt.headerlines)]
  pretty.labels <- system(paste("cut -f2", cdt.file), intern=TRUE)
  pretty.labels <- pretty.labels[-(1:cdt.headerlines)]
  
  get.tree.index <- function(x, neg, pos) {
    nm <- -match(x, neg)
    pm <- match(x, pos)
    
    if(any(is.na(nm) == is.na(pm))) stop("overlapping gene and node names")
    pm[is.na(pm)] <- nm[!is.na(nm)]
    return(pm)
  }
  merge.col1 <- get.tree.index(tree.data[,2], gene.labels, merge.labels)
  merge.col2 <- get.tree.index(tree.data[,3], gene.labels, merge.labels)

  list(merge=cbind(merge.col1, merge.col2),
       labels=pretty.labels,
       height=(1-tree.data[,4])) # this is for Pearson correlation...
}



read.cluster.tree <- function(tree.file, gene.labels, pretty.labels) {
  tree.data <- read.delim(tree.file, header=FALSE)
  merge.labels <- tree.data[,1]
  
  get.tree.index <- function(x, neg, pos) {
    nm <- -match(x, neg)
    pm <- match(x, pos)
    
    if(any(is.na(nm) == is.na(pm))) stop("overlapping gene and node names")
    pm[is.na(pm)] <- nm[!is.na(nm)]
    return(pm)
  }
  merge.col1 <- get.tree.index(tree.data[,2], gene.labels, merge.labels)
  merge.col2 <- get.tree.index(tree.data[,3], gene.labels, merge.labels)

  list(merge=cbind(merge.col1, merge.col2),
       labels=pretty.labels,
       height=tree.data[,4])
}

main()

