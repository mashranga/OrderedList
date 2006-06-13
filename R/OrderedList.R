### compute overlap
overlap <- function(x1,x2,n){
  r <- match(x1[1:n],x2[1:n])
  select <- r <= n
  overlapRanks <- pmax(r[select], 1:n)
  tmp <- table(overlapRanks)
  x <- integer(n)
  x[as.integer(names(tmp))] <- tmp
  x <- cumsum(x)
  return(x)
}

# helper functions for test statistics on matrices
check.test.args <- function(m, cl, paired) {
  lbls <-unique(cl)
  if (length(lbls) != 2)
    stop("Exactly two labels needed, received instead: ", lbls)
  if (ncol(m) != length(cl))
    stop("Exactly one class label needed per column, received ", length(c), 
         " labels for ", ncol(m), " columns instead")
  if (paired) {
    if (sum(cl==lbls[1]) != sum(cl==lbls[2])) 
      stop("Groups must have same size in paired test, recieved ", sum(cl==lbls[1]), " ",
           lbls[1], " samples and ", sum(cl==lbls[2]), " ", lbls[2], " samples instead")
  }
  return(lbls)
}

test.fc <- function(m, cl, paired=FALSE) {
  lbls <- check.test.args(m, cl, paired)
  if (paired) fc <- rowMeans(m[,cl==lbls[1]] - m[,cl==lbls[2]])
  else        fc <- rowMeans(m[,cl==lbls[1]]) - rowMeans(m[,cl==lbls[2]])
  return(fc)
}

test.t <- function(m, cl, paired=FALSE) {
  lbls <- check.test.args(m, cl, paired)
  cll <- integer(length(cl))
  cll[cl==cl[1]] <- 1
  t.stat <- twilight.teststat(xin=m, yin=cll, method="t", paired=paired)$observed
  names(t.stat) <- rownames(m)
  return(t.stat)
}

test.z <- function(m, cl, paired=FALSE) {
  lbls <- check.test.args(m, cl, paired)
  cll <- integer(length(cl))
  cll[cl==cl[1]] <- 1
  z.stat <- twilight.teststat(xin=m, yin=cll, method="z", paired=paired)$observed
  names(z.stat) <- rownames(m)
  return(z.stat)
}

# score rankings by test-statistic on sampled data
scoreOrderComparison <- function(exprs1, labels1, paired1,  exprs2, labels2, paired2, 
                                 test.method=test.z, nn, bases, two.sided) {
  # computes scores for matched direction only
  
  r1 <- rank(test.method(exprs1, labels1, paired1))
  r2 <- rank(test.method(exprs2, labels2, paired2))
  return(scoreRankings(r1, r2, nn, bases, two.sided))
}

scoreOrderComparisonBoth <- function(exprs1, labels1, paired1,  exprs2, labels2, paired2, 
                                 test.method=test.z, nn, bases, two.sided) {
  # computes scores for both, direct and reversed orders.
  
  n <- nrow(exprs1)
  r1 <- rank(test.method(exprs1, labels1, paired1))
  r2 <- rank(test.method(exprs2, labels2, paired2))
  return(c(scoreRankings(r1, r2, nn, bases, two.sided),
           scoreRankings(r1, n+1-r2, nn, bases, two.sided)))
}

# compare lists with resampling
# -----------------------------
preparePermutations <- function(ids, paired, B, sample.ratio=0.8) {
  # determine number of needed subsamples
  nbad <- floor(sum(ids)*sample.ratio)
  ngood <- floor(sum(1-ids)*sample.ratio)

  # prepare matrices
  if (paired){
    yperm <- twilight.permute.pair(ids,B+1,bal=FALSE)
    ysubs <- matrix(nrow=B,ncol=(2*nbad))
    for (i in 1:B){
      x <- sample(1:(sum(ids)),nbad)
      ysubs[i,] <- c(x,sum(ids)+x)
    }  
  } else { # not paired
    yperm <- twilight.permute.unpair(ids,B+1,bal=FALSE)
    ysubs <- matrix(nrow=B,ncol=(nbad+ngood))
    for (i in 1:B){
      x <- sample(1:(sum(ids)),nbad)
      y <- sample((sum(ids)+1):(length(ids)),ngood)
      ysubs[i,] <- c(x,y)
    }  
  }
  return(list(yperm=yperm, ysubs=ysubs))
}

OrderedList <- function(eset, B=1000, test="z", beta=1, percent=0.95, verbose=TRUE,
                        alpha=NULL, min.weight=1e-5){

  ### check arguments

  if (class(eset)!="exprSet"){
    stop("Please prepare the input expression set with function 'prepareData'.")
  }
  
  if (sum(names(pData(eset))%in%c("outcome","dataset","class","paired"))!=4){
    stop("Please prepare the input expression set with function 'prepareData'.")
  }

  if(test%in%c("fc","t","z")){
    if (test == "fc") current.test <- test.fc
    if (test == "t" ) current.test <- test.t
    if (test == "z" ) current.test <- test.z
  } else {
    stop("'test' can only bet set to 'fc','t' or 'z'.")
  }

  if (!beta%in%c(1,0.5)){
    stop("'beta' can only be set to 1 or 0.5.")
  }

  if ((percent<=0)|(percent>1)){
    stop("'percent' has to be in the interval (0,1].")
  }

  if (!is.null(alpha)) {
    if (any(!is.numeric(alpha)))
      stop("Values for alpha must be numeric, received ", alpha[!is.numeric(alpha)], 
           " instead")
  }

  ### store comparison label

  pdata <- pData(eset)
  x <- levels(pdata$dataset)
  label <- paste(x[1],x[2],sep="~")

  ### prepare input data

  data1 <- pdata$dataset==x[1]
  data2 <- pdata$dataset==x[2]

  eset1 <- exprs(eset)[,data1]
  eset2 <- exprs(eset)[,data2]
  
  ngenes <- nrow(eset1)
  probes <- rownames(eset1)
  
  id1 <- as.character(pdata[data1,]$class)
  id2 <- as.character(pdata[data2,]$class)
  
  paired1 <- as.logical(unique(pdata[data1,]$paired))
  paired2 <- as.logical(unique(pdata[data2,]$paired))
  
  id1[id1=="bad"] <- 1
  id2[id2=="bad"] <- 1
  id1[id1=="good"] <- 0
  id2[id2=="good"] <- 0
  
  id1 <- as.numeric(id1)
  id2 <- as.numeric(id2)
  
  ### prepare permutation matrices

  samplings1 <- preparePermutations(id1, paired1, B)
  samplings2 <- preparePermutations(id2, paired2, B)
  
 
  ### prepare weights

  if (is.null(alpha)) {
    nn <- c(100, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000, 2500)
    alpha <- -log(min.weight)/nn
  } else {
    nn <- floor(-log(min.weight)/alpha)
  }
  select <- nn < ngenes
  alpha <- alpha[select]
  nn <- nn[select]
  bases <- exp(-alpha)
  nalpha <- length(alpha)
  nmax <- max(nn)

  if (beta == 0.5) scoreOrderComparison <- scoreOrderComparisonBoth

  ### score the actual orders

  SIM.obs <- scoreOrderComparison(eset1, id1, paired1, eset2, id2, paired2,
                                  current.test, nn, bases, TRUE)

  ### compute distributions of alternative and null-scores

  samplings1$yperm <- samplings1$yperm[-1,]
  samplings2$yperm <- samplings2$yperm[-1,]

  if (verbose) {
    cat("\nSimulating score distributions...\n")
    if (B < 50) cat("          0%", rep(".",B-6), "100%\n", sep="")
    else cat("          0%.......:.........:.........:.........:......100%\n")
    cat("  Random: ")
    dotStep <- ceiling(B/50)
  } else dotStep=Inf

  indirectScore <- function(i,yin1,xin1,paired1,yin2,xin2,paired2,
                            current.test, nn, bases, dotStep=Inf) {
    if ((i %% dotStep)==0) cat("-")
    scores <- scoreOrderComparison(xin1,yin1[i,], paired1, xin2,yin2[i,], paired2,
                                   current.test, nn, bases, TRUE)
    return(scores)
  }
  SIM.random <- t(sapply(1:B,indirectScore,samplings1$yperm,eset1,paired1,
                                               samplings2$yperm,eset2,paired2,
                                               current.test,nn,bases,dotStep))

  if (verbose) cat("\nObserved: ")
  subsetScore <- function(i,index1,yin1,xin1,paired1,index2,yin2,xin2,paired2,
                          current.test, nn, bases, dotStep=Inf) {
    if ((i %% dotStep)==0) cat("-")
    scores <- scoreOrderComparison(xin1[,index1[i,]],yin1[index1[i,]], paired1, 
                                   xin2[,index2[i,]],yin2[index2[i,]], paired2,
                                   current.test, nn, bases, TRUE)
    return(scores)
  }
  SIM.observed <- t(sapply(1:B,subsetScore,samplings1$ysubs,id1,eset1,paired1,
                                               samplings2$ysubs,id2,eset2,paired2,
                                               current.test,nn,bases,dotStep))
  if (verbose) cat("\n")
  
  ### compute pAUC scores of separability

  pauc <- function(x,A,B){
#    u <- sort(unique(A),decreasing=TRUE)
#    t <- numeric(length(u))
#    r <- numeric(length(u))
#    for (i in 1:length(u)){
#      t[i] <- sum(A[B==0]>u[i])/sum(1-B)
#      r[i] <- sum(A[B==1]>u[i])/sum(B)
#    }

    n <- length(A)
    o <- order(A, decreasing=TRUE)
    r <- c(0, cumsum(B[o[-n]]))
    t <- (0:(n-1)) - r
    r <- r/sum(B)
    t <- t/sum(1-B)

    roc <- numeric(length(x))
    for (i in 1:length(x)){
      z <- which(t<=x[i])
      z <- z[length(z)]
      roc[i] <- r[z]
    }
    
    return(roc)
  }

  pAUC <- numeric(length(SIM.obs))
  y <- c(rep(0,B),rep(1,B))
  for (i in 1:length(SIM.obs)){
    X <- c(SIM.random[,i],SIM.observed[,i])
    pAUC[i] <- integrate(pauc,0,0.1,A=X,B=y,stop.on.error=FALSE)$value
  }
  
  ### select optimal alpha and direction

  xx <- which.max(pAUC)
  alpha.opt <- alpha[ifelse(xx > nalpha, xx-nalpha, xx)]
  direction <- ifelse(xx > nalpha, -1, 1)

  ### prepare output

  res <- list()
  res$SIM.observed    <- SIM.obs[xx]
  res$SIM.alternative <- SIM.observed[,xx]
  res$SIM.random      <- SIM.random[,xx]
  res$subSample       <- TRUE

  
  pauc <- matrix(pAUC,ncol=nalpha, byrow=TRUE)
  colnames(pauc) <- alpha
  if (beta == 1) rownames(pauc) <- "direct"
  else rownames(pauc) <- c("direct","reversed")

  x1 <- current.test(eset1, id1, paired1)
  x2 <- current.test(eset2, id2, paired2)
  scores <- cbind(x1,x2)
  colnames(scores) <- unlist(strsplit(label,"~"))
  rownames(scores) <- rownames(eset1)

  order1 <- rownames(eset1)[order(x1, decreasing=TRUE)]
  order2 <- rownames(eset2)[order(x2, decreasing=TRUE)]
  if (direction < 0) order2 <- rev(order2)

#  pp <- c(res$SIM.observed,res$SIM.random)
  pp <- c(median(res$SIM.alternative),res$SIM.random)
  pp <- rank(pp)
  pvalue <- (length(pp)-pp[1]+1)/length(pp)
  
  ### retrieve intersecting probes

  nn <- ifelse (direction > 0, nn[xx], nn[xx-nalpha])
  ww <- exp(-alpha.opt*(1:nn))
  rr <- ww*(overlap(order1, order2, nn) + overlap(rev(order1), rev(order2), nn))
  y <- min(which(cumsum(rr) >= percent*sum(rr)-(1e-5)))

  y1 <- intersect(order1[1:y], order2[1:y])
  y2 <- intersect(order1[(ngenes-y+1):ngenes], order2[(ngenes-y+1):ngenes])

  ### prepare final output

  x <- list(n=ngenes,
            label=label,
            p=pvalue,
            intersect=sort(c(y1,y2)),
            alpha=alpha.opt,
            direction=direction,
            scores=scores,
            sim.scores=res,
            pauc=pauc,
            call=list(B=B,test=test,beta=beta,percent=percent,
                      alpha=alpha,min.weight=min.weight,
                      labels1=as.character(unique(pdata[data1,"outcome"])),
                      labels2=as.character(unique(pdata[data2,"outcome"])))
            )
  
  class(x) <- "OrderedList" 
  return(x)
}
