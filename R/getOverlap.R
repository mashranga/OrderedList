### get the overlapping genes from a listComparison object
getOverlap <- function(x,max.rank=100,percent=0.95){

  if (class(x)!="listComparison"){stop("Input must be a 'listComparison' object.")}

  if (!max.rank%in%x$nn){stop("Choose 'max.rank' from the evaluated gene numbers.")}

  nn <- max.rank
  NN <- length(x$overlaps)
  ngenes <- x$n

  i <- which(x$nn==nn)
  
  alpha <- x$call$alphas[i]
  ww <- exp(-alpha*(1:nn))

  res <- list()
  res$n <- x$n
  res$call <- x$call
  res$call$alphas <- alpha
  res$nn <- nn
  
  overlaps <- x$overlaps[1:nn]
  revOverlaps <- x$revOverlaps[1:nn]

  score <- ww*(overlaps[1:nn])
  revScore <- ww*(revOverlaps[1:nn])

  if (x$call$two.sided){
    overlaps <- c(overlaps, x$overlaps[(NN-nn+1):NN])
    revOverlaps <- c(revOverlaps, x$revOverlaps[(NN-nn+1):NN])

    score <- ww*(overlaps[1:nn]+rev(overlaps[(nn+1):(2*nn)]))
    revScore <- ww*(revOverlaps[1:nn]+rev(revOverlaps[(nn+1):(2*nn)]))
  }
  
  
  res$score <- x$scores[i]
  res$pvalue <- x$pvalues[i]
  res$overlaps <- overlaps
  res$randomScores <- x$randomScores[,i]
  res$direction <- 1
  
  ## which direction?
  if (sum(score)>sum(revScore)){

    y <- min(which(cumsum(score) >= percent*sum(score)-(1e-5)))
    y1 <- intersect(x$ID.List1[1:y], x$ID.List2[1:y])

    if (x$call$two.sided){
      y1 <- c(y1, intersect(x$ID.List1[(ngenes-y+1):ngenes], x$ID.List2[(ngenes-y+1):ngenes]))
    }
  }
  else {
    res$score <- x$revScores[i]
    res$pvalue <- x$revPvalues[i]
    res$overlaps <- revOverlaps
    res$direction <- -1
    
    y <- min(which(cumsum(revScore) >= percent*sum(revScore)-(1e-5)))    
    y1 <- intersect(x$ID.List1[1:y], rev(x$ID.List2)[1:y])

    if (x$call$two.sided){
      y1 <- c(y1, intersect(x$ID.List1[(ngenes-y+1):ngenes], rev(x$ID.List2)[(ngenes-y+1):ngenes]))
    }
  } 
  
  res$call$percent <- percent
  res$intersect <- sort(y1)
  class(res) <- "listComparisonOverlap"
  
  return(res)
}


plot.listComparisonOverlap <- function(x, which="overlap",...) {
  if (which=="overlap"){
    overlaps <- x$overlaps
    N <- x$n
    NN <- length(overlaps)
    nn <- x$nn
    ylims <- range(overlaps)
    xlims <- c(1,nn)
    if (x$call$two.sided){xlims[2] <- 2*nn}

    mu <- (1:nn)^2/N
    upper <- qhyper(0.975,1:nn,N-(1:nn),1:nn)
    lower <- qhyper(0.025,1:nn,N-(1:nn),1:nn)
    
    polyX <- c(1:nn,nn:1)
    polyY <- c(upper,rev(lower))
    polyX <- rep(polyX,rep(2,length(polyX)))
    polyY <- rep(polyY,rep(2,length(polyY)))
    polyX <- polyX[-1]
    polyY <- polyY[-length(polyY)]
    
    plot(0,type="n",xaxt="n", xlim=xlims, ylim=ylims, 
         xlab="",ylab="size of overlap",
         main=paste("Overlap for alpha:",round(x$call$alphas,3)))

    par(fg="orange")
    polygon(polyX,polyY,angle=90,density=30,border=NA)
    par(fg="black")

    lines(1:nn,mu,col="orange",t="s")
    lines(1:nn,overlaps[1:nn],t="s")
    
    if (x$call$two.sided) {
      polyX <- c((nn+1):(2*nn),rev((nn+1):(2*nn)))
      polyY <- c(rev(upper),lower)
      polyX <- rep(polyX,rep(2,length(polyX)))
      polyY <- rep(polyY,rep(2,length(polyY)))
      polyX <- polyX[-1]
      polyY <- polyY[-length(polyY)]

      par(fg="orange")
      polygon(polyX,polyY,angle=90,density=30,border=NA)
      par(fg="black")

      lines((nn+1):(2*nn),rev(mu),col="orange",t="s")
      lines((nn+1):(2*nn),overlaps[(NN-nn+1):NN],t="s")
      
      axStep <- floor(nn / 45) * 10
      axis(side=1, at=c(axStep*(0:4), 2*nn-axStep*(0:4)), 
           labels=c(1, axStep*(1:4), 1, axStep*(1:4)))
      mtext("top ranks", side=1, at=0.5*nn, line=2.5)
      mtext("bottom ranks", side=1, at=1.5*nn, line=2.5)
      abline(v=nn+0.5,lty=2)
    }
    else { # one sided
      axis(side=1)
      mtext("top ranks", side=1, at=nn/2, line=2.5)
    }
    legend(1, ylims[2], legend=c("observed", "expected"), col=c("black", "orange"), lty=1, bty="n")
    return(invisible())
  }
  if (which=="density"){
    d <- density(x$randomScores, from=0)
    plot(d, main=paste("Random Scores for alpha:", round(x$call$alphas,3)), xlab="Similarity score")
    abline(v=x$score)
    rug(x$randomScores)
    mtext("observed", side=3, at=x$score, line=0.1)
    return(invisible())
  }
  stop("Value for parameter 'which' not recognized: ", which)
}


print.listComparisonOverlap <- function(x, ...) {
  cat("List comparison")
  cat("\n  Assessing similarity of  :", ifelse(x$call$two.sided, "top and bottom ranks", "top ranks"))
  cat("\n  Length of lists          :", x$n)
  cat("\n  Number of random samples :", x$call$B)
  cat("\n--------------------------------------")
  cat("\n  Chosen regularization parameter: alpha =", x$call$alphas)
  cat("\n  Lists are more alike in", ifelse(x$direction > 0,
                                            "direct", "reversed"), "order")
  cat("\n  Weighted overlap score: ", x$score)
  cat("\n  Significance of similarity: p-value =", x$pvalue)
  cat("\n  Number of genes contributing", x$call$percent * 100,
      "% to similarity score:", length(x$intersect))
  cat("\n")
}

