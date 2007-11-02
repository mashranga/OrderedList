prepareData <- function(eset1,eset2,mapping=NULL){

	if(!is(eset1$data,"eSet"))
	{
	stop("Object '",eset1$name, "' needs to inherit from 'eSet', but is of class '",class(eset1$data),"' instead !\n")
	}
	if(!is(eset2$data,"eSet"))
	{
	stop("Object '",eset2$name, "' needs to inherit from 'eSet', but is of class '",class(eset2$data),"' instead !\n")
	}


  x <- list(eset1,eset2) 
  n <- length(x) 

  ### Some sanity checks first.
  if (is.null(mapping)){
    M <- exprs(x[[1]]$data) # methodCall 'exprs()' of class:ExpressionSet
    nn <- rownames(M) # get rownames of M
    for (i in 2:n){
      M <- exprs(x[[i]]$data)
      if (sum(rownames(M)%in%nn)!=length(nn)){stop("Probe ID differ between studies. Please specify a mapping.\n")}
    }
  }

  if (!is.null(mapping)){
    if (sum(is.na(mapping))>0){stop("'mapping' must not contain missing values.\n")} ### any NA-values?
    nn <- length(mapping) 	### nr of elements
    mm <- length(mapping[[1]]) 	### value 'mm' is never used !!
    if (nn!=n){stop("Number of mappings does not match number of studies.\n")}
    for (i in 1:n){ 	### n= nr of studies, here: 2
      M <- exprs(x[[i]]$data)
      ww <- mapping[x[[i]]$name] 	### get element of mapping by index 'name' of x
      ww <- unlist(ww[[1]]) # convert to vector to gain atomic values
      select <- which(ww%in%rownames(M)) 	### get index of elements that are TRUE
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

  ### Merge studies into one ExpressionSet (deprecated: exprSet).
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

    eset <- cbind(eset,M[,x1],M[,x2]) # merge datasets into one ExpressionSet instance
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

   

   desc <- data.frame(labelDescription=c("outcome","dataset","class","paired"))
   
   
  rownames(desc) <- colnames(p)
  

  pdata <- new("AnnotatedDataFrame", data=p, varMetadata=desc) ### data:=data.frame ; varMetadata:=data.frame
  eset  <- new("ExpressionSet", exprs=eset, phenoData=pdata) ### exprs:= matrix ; phenoData:=AnnotatedDataFrame 


  return(eset)

}
