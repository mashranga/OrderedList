### get the overlapping genes from a listComparison object
getOverlap <- function(x,max.rank=NULL,percent=0.95){

  if (class(x)!="listComparison"){stop("Input must be a 'listComparison' object.")}

  if (is.null(max.rank)) {
    # find the rank where the best p-value is achieved
#    tmp <- which.min(c(x$pvalue, x$revPvalue)) # which.min = R !
    tmpx <- c(x$pvalue, x$revPvalue)
    tmp <- max(which(tmpx == min(tmpx)))
    max.rank <- x$nn[(tmp-1)%%length(x$nn)+1]
  }
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
  
  res$direction <- ifelse(x$scores[i] > x$revScores[i], 1, -1)
  res$score <- max(x$scores[i], x$revScores[i])
  res$pvalue <- min(x$pvalues[i], x$revPvalues[i])
  res$randomScores <- x$randomScores[,i]
  
  ## which direction?
  if (res$direction > 0) {
    res$overlaps <- x$overlaps

    y <- min(which(cumsum(score) >= percent*sum(score)-(1e-5)))
    y1 <- intersect(x$ID.List1[1:y], x$ID.List2[1:y])

    if (x$call$two.sided){
      y1 <- c(y1, intersect(x$ID.List1[(ngenes-y+1):ngenes], 
                            x$ID.List2[(ngenes-y+1):ngenes]))
    }
  }
  else {
    res$overlaps <- x$revOverlaps
    
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


plot.listComparisonOverlap <- function(x, which="overlap", no.title=FALSE, no.legend=FALSE,
                                       list.name1="List 1", list.name2="List 2", ...) {
  if (no.title) par(mai=par("mai")[c(1,2,4,4)])
  if (which=="overlap"){
    overlaps <- x$overlaps
    N <- trunc((1-x$call$invar.q) * x$n)
    NN <- length(overlaps)
    nn <- x$nn

    mu <- (1:nn)^2/N
    Ns <- N-(1:nn); Ns <- ifelse(Ns > 0, Ns, 0)
    upper <- qhyper(0.975,1:nn,Ns,1:nn)
    lower <- qhyper(0.025,1:nn,Ns,1:nn)

    if (x$call$two.sided){
      ylims <- range(c(overlaps[c(1:nn, (NN-nn+1):NN)], lower, upper))
      xlims <- c(1,2*nn)
    } else {
      ylims <- range(c(overlaps[1:nn], lower, upper))
      xlims <- c(1,nn)
    }
    
    polyX <- c(1:nn,nn:1)
    polyY <- c(upper,rev(lower))
    polyX <- rep(polyX,rep(2,length(polyX)))
    polyY <- rep(polyY,rep(2,length(polyY)))
    polyX <- polyX[-1]
    polyY <- polyY[-length(polyY)]
    
    plot(0,type="n",xaxt="n", xlim=xlims, ylim=ylims, 
         xlab="",ylab="size of overlap", ...)
    if (!no.title) title(main=paste("Overlap for alpha:",round(x$call$alphas,3)))
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
    }
    if (!no.legend) {
      if (x$direction > 0) {
        text(1, ylims[2], adj=c(0,1),
             label=paste("top in ", list.name1, "\ntop in ", list.name2, sep=""))
        text(2*nn, ylims[2], adj=c(1,1),
             label=paste("bottom in ", list.name1, "\nbottmo in ", list.name2, sep=""))
      } else {
        text(1, ylims[2], adj=c(0,1),
             label=paste("top in ", list.name1, "\nbottom in ", list.name2, sep=""))
        text(2*nn, ylims[2], adj=c(1,1),
             label=paste("bottom in ", list.name1, "\ntop in ", list.name2, sep=""))
      }
      legend(nn+1, ylims[1], yjust=0, lty=1, bty="n",
             legend=c("observed", "expected"), col=c("black", "orange"))
    }
    return(invisible())
  }
  if (which=="scores"){
    d <- density(x$randomScores, from=0)
    xmax <- max(c(x$randomScores, x$score))
    plot(d, xlab="similarity score", main="", xlim=c(0,xmax), ...)
    if (!no.title) title(main="Distribution of Random Scores")
    abline(v=x$score)
    rug(x$randomScores)
    mtext("observed", side=3, at=x$score, line=0.1)
    if (!no.legend) {
      text(0, max(d$y), adj=c(0,1),
           label=paste("direction: ", ifelse(x$direction>0, "direct", "reversed"), 
                       "\np-value: ", x$pvalue, 
                       "\nalpha: ",   round(x$call$alphas,3),
                       "\nlist size: ", x$n,
                       "\ninvariants: ", round(x$call$invar.q*100,0), "%", sep=""))
    }
    return(invisible())
  }
  stop("Value for parameter 'which' not recognized: ", which)
}


print.listComparisonOverlap <- function(x, ...) {
  cat("List comparison")
  cat("\n  Assessing similarity of               :", ifelse(x$call$two.sided, "top and bottom ranks", "top ranks"))
  cat("\n  Length of lists                       :", x$n)
  cat("\n  Number of random samples              :", x$call$B)
  cat("\n----------------------------------------------------------")
  cat("\n  Lists are more alike in", ifelse(x$direction > 0,
                                            "direct", "reversed"), "order")
  cat("\n  Chosen regularization parameter       : alpha =", round(x$call$alphas,3),
      "(", x$nn, "genes)")
  cat("\n  Weighted overlap score                :", x$score)
  cat("\n  Significance of similarity            : p-value =", x$pvalue)
  cat("\n  Score percentage for common entries   :", x$call$percent*100)
  cat("\n  Entries contributing score percentage :", length(x$intersect))
  cat("\n")
}

