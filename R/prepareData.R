prepareData <- function(eset1,eset2,mapping=NULL){
  x <- list(eset1,eset2)
  n <- length(x)

  ### Some sanity checks first.
  if (is.null(mapping)){
    M <- exprs(x[[1]]$data)
    nn <- rownames(M)
    for (i in 2:n){
      M <- exprs(x[[i]]$data)
      if (sum(rownames(M)%in%nn)!=length(nn)){stop("Probe ID differ between studies. Please specify a mapping.\n")}
    }
  }

  if (!is.null(mapping)){
    if (sum(is.na(mapping))>0){stop("'mapping' must not contain missing values.\n")}
    nn <- length(mapping)
    mm <- length(mapping[[1]])
    if (nn!=n){stop("Number of mappings does not match number of studies.\n")}
    for (i in 1:n){
      M <- exprs(x[[i]]$data)
      ww <- mapping[x[[i]]$name]
      ww <- unlist(ww[[1]])
      select <- which(ww%in%rownames(M))
      mapping <- mapping[select,]
    }
  }

  # modify names of studies if needed
  dataset.names <- rep("", n)
  for (i in 1:n) dataset.names[i] <- x[[i]]$name
  unique.names <- make.unique(dataset.names)
  if (!identical(dataset.names, unique.names)) {
    if (is.null(mapping)){
      for (i in 1:n) x[[i]]$name <- unique.names[i]
    } else {
      stop("Some dataset names are equal: mapping is not yet implemented")
    }
  }

  ### Merge studies into one exprSet.
  eset <- numeric()
  outcome <- character()
  dataset <- character()
  class <- character()
  paired <- character()
  
  for (i in 1:n){
    M <- exprs(x[[i]]$data)
    colnames(M) <- paste(colnames(M), i, sep=".")
    p <- pData(x[[i]]$data)
    rownames(p) <- colnames(M)
    
    if (!is.null(mapping)){
      ww <- mapping[x[[i]]$name]
      ww <- as.character(unlist(ww[[1]])      )
      M <- M[ww,]
    }
    
    out <- x[[i]]$out
    if (length(out)!=2){stop("'out' must be of length 2.\n")}

    var <- x[[i]]$var
    x1 <- rownames(p)[p[,var]%in%out[1]] ### bad outcome
    x2 <- rownames(p)[p[,var]%in%out[2]] ### good outcome

    eset <- cbind(eset,M[,x1],M[,x2])
    outcome <- c(outcome,rep(out[1],length(x1)),rep(out[2],length(x2)))
    class <- c(class,rep("bad",length(x1)),rep("good",length(x2)))
    dataset <- c(dataset,rep(x[[i]]$name,length(x1)+length(x2)))
    paired <- c(paired,rep(x[[i]]$paired,length(x1)+length(x2)))
  }
  if (!is.null(mapping)) {
    rownames(eset) <- as.character(apply(mapping,1,paste,collapse="/"))
  }

  p <- data.frame(outcome=outcome,dataset=dataset,class=class,paired=paired)
  rownames(p) <- colnames(eset)
  desc <- list("outcome","dataset","class","paired")
  names(desc) <- names(p)

  pdata <- new("phenoData", pData=p, varLabels=desc) 
  eset  <- new("exprSet", exprs=eset, phenoData=pdata)

  return(eset)
}
