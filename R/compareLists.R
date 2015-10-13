scoreRankings <- function(r1, r2, nn, bases, two.sided=TRUE) {
  n <- length(r1)
  if (length(r2)!=n) stop("Rankings have unequal length")
  nbases <- length(bases)

  nmax <- max(nn)
  select <- r1 <= nmax & r2 <= nmax
  overlapRanks <- pmax(r1[select], r2[select])

  if (two.sided) {
    nmin <- n - nmax
    select <- r1 >= nmin & r2 >= nmin
    overlapRanks <- c(overlapRanks, n+1 - pmin(r1[select], r2[select]))
  }

  scores <- numeric(nbases)
  for (i in 1:nbases) {
    select <- overlapRanks <= nn[i]
    nselect <- sum(select)
    scores[i] <- sum(bases[i]^overlapRanks[select], -nselect*bases[i]^(nn[i]+1))
  }
  scores <- scores/(1-bases)

  return(scores)
}

shuffledRandomScores <- function(n, nn, bases, B=1000, two.sided=TRUE) {
  randomScores <- matrix(0, nrow=B, ncol=length(bases))
  dotStep <- ceiling(B/50)
  cat("  Simulating random scores...\n")
  if (B < 50) cat("  0%", rep(".",B-6), "100%\n  ", sep="")
  else cat("  0%.......:.........:.........:.........:......100%\n  ")
  for (i in 1:B) {
    rList <- sample(n)
    randomScores[i,] <- scoreRankings(rList, 1:n, nn, bases, two.sided)
    if ((i %% dotStep)==0) cat("-")
  }
  cat("\n")
  colnames(randomScores) <- bases
  class(randomScores) <- "shuffledRandomScores"
  return(randomScores)
}

print.shuffledRandomScores <- function(x, ...) {
  cat("Random scores from shuffled lists:\n")
  cat("  alphas  :", round(-log(as.numeric(colnames(x))), 3), "\n")
  cat("  nalphs ::", ncol(x), "\n")
  cat("  samples :", nrow(x), "\n")
}

plot.shuffledRandomScores <- function(x, observeds=NULL, revObserveds=NULL, ...) {
  for (i in 1:ncol(x)) {
    d <- density(x[,i], from=0)
    plot(d, main=paste("Random Scores for alpha:", 
              round(-log(as.numeric(colnames(x)[i])), 3)), xlab="Similarity score")
    rug(x[,i])
    if (!is.null(observeds)) {
      abline(v=observeds[i])
      mtext("observed", side=3, at=observeds[i], line=0.1)
    }
    if (!is.null(revObserveds)) {
      abline(v=revObserveds[i])
      mtext("revObserved", side=3, at=revObserveds[i], line=0.1)
    }
    readline("<hit enter for next alpha>")
  }
}


compareLists <- function(ID.List1, ID.List2, mapping=NULL, 
                         two.sided=TRUE, B=1000, alphas=NULL, invar.q=0.5,
                         min.weight=1e-5, no.reverse=FALSE) {

  res <- list()

  # checkout mapping
  if (!is.null(mapping)) {
    res$mapping <- list()
    tmp <- ID.List1 %in% mapping[,1]
    if (any(!tmp)) 
      cat(sum(!tmp)," of ",length(tmp)," elements in first list not found in mapping\n")
    res$mapping$missed1 <- sum(!tmp)
    res$mapping$rawlen1 <- length(tmp)
    tmp <- ID.List2 %in% mapping[,2]
    if (any(!tmp)) 
      cat(sum(!tmp)," of ",length(tmp), " elements in second list not found in mapping\n")
    res$mapping$missed2 <- sum(!tmp)
    res$mapping$rawlen2 <- length(tmp)

    # restrict map
    mapping <- mapping[mapping[,1] %in% ID.List1, ,drop=FALSE]
    mapping <- mapping[mapping[,2] %in% ID.List2, ,drop=FALSE]
    n <- nrow(mapping)

    # convert list IDs into combined mapping-IDs
    mapIDs <- apply(mapping, 1, paste, collapse="/")
    oo <- order(match(mapping[,1], ID.List1))
    ID.List1 <- mapIDs[oo]
    oo <- order(match(mapping[,2], ID.List2))
    ID.List2 <- mapIDs[oo]
  }

  # check list consistency
  tmp <- sum(!(ID.List1 %in% ID.List2))
  if (tmp > 0) stop(tmp, " element(s) of first list not found in second")
  tmp <- sum(!(ID.List2 %in% ID.List1))
  if (tmp > 0) stop(tmp, " element(s) of second list not found in first")

  # determine ranks
  n <- length(ID.List2)
  Ranks.List1 <- match(ID.List1, ID.List2)
  Ranks.List2 <- 1:n

  # initialize
  if (is.null(alphas)) {
    nn <- c(100, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000, 2500)
    alphas <- -log(min.weight)/nn
  } else {
    nn <- floor(-log(min.weight)/alphas)
  }
  select <- nn < n
  alphas <- alphas[select]
  nn <- nn[select]
  nalphas <- length(alphas)

  res$n          <- n
  res$call       <- list(B=B,alphas=alphas,invar.q=invar.q,
                         two.sided=two.sided,min.weight=min.weight,
                         no.reverse=no.reverse)
  res$nn         <- nn
  res$scores     <- numeric(nalphas)
  res$revScores  <- numeric(nalphas)
  res$pvalues    <- numeric(nalphas)
  res$revPvalues <- numeric(nalphas)
  class(res) <- "listComparison"

  res$overlaps <- overlap(ID.List1, ID.List2, max(nn))
  if (no.reverse) { res$revOverlaps <- rep(0, length(res$overlaps)) } 
  else { res$revOverlaps <- overlap(ID.List1, rev(ID.List2), max(nn)) }
  if (two.sided) {
    res$overlaps <- c(res$overlaps, rev(overlap(rev(ID.List1), rev(ID.List2), max(nn))))
    if (no.reverse) { res$revOverlaps <- rep(res$revOverlaps, 2) }
    else { res$revOverlaps <- c(res$revOverlaps, rev(overlap(rev(ID.List1), ID.List2, max(nn)))) }
  }

  # compute observed
  bases <- exp(-alphas)
  res$scores <- scoreRankings(Ranks.List1, Ranks.List2, nn, bases, two.sided)
  if (no.reverse) { res$revScores <- rep(0, length(res$scores)) }
  else { res$revScores <- scoreRankings(Ranks.List1, n+1-Ranks.List2, nn, bases, two.sided) }

  # compute random distributions and p-values
  res$randomScores <- shuffledRandomScores(trunc(n*(1-invar.q)), nn, bases, B, two.sided)
  for (i in 1:nalphas) {
    res$pvalues[i] <- sum(res$randomScores[,i] > res$scores[i])/B
    if (no.reverse) { res$revPvalues[i] <- 1 }
    else { res$revPvalues[i] <- sum(res$randomScores[,i] > res$revScores[i])/B }
  }

  res$ID.List1 <- ID.List1
  res$ID.List2 <- ID.List2
  
  return(res)
}

print.listComparison <- function(x, ...) {
  cat("List comparison")
  cat("\n  Assessing similarity of     :", ifelse(x$call$two.sided, "top and bottom ranks", "top ranks"))
  if (is.null(x$mapping)) {
    cat("\n  Length of lists             :", x$n)
  } else {
    cat("\n  Mapped entries from list 1  :", 
        x$mapping$rawlen1-x$mapping$missed1, "of", x$mapping$rawlen1)
    cat("\n  Mapped entries from list 2  :", 
        x$mapping$rawlen2-x$mapping$missed2, "of", x$mapping$rawlen2)
    cat("\n  Number of mapping pairs     :", x$n)
  }
  cat("\n  Quantile of invariant genes :", x$call$invar.q)
  cat("\n  Number of random samples    :", x$call$B)
  cat("\n--------------------------------------")
  m <- data.frame(Genes=x$nn, 
                  Scores=x$scores,
                  p.values=x$pvalues,
                  Rev.Scores=x$revScores,
                  Rev.p.values=x$revPvalues)
  rownames(m) <- round(x$call$alphas,3)
  cat("\n")
  print(m)
}

plot.listComparison <- function(x, which="overlap",...) {
  if (which == "overlap") {
    for (i in 1:length(x$nn)){
      N <- x$n
      NN <- length(x$overlaps)
      nn <- x$nn[i]
      ylims <- range(x$overlaps, x$revOverlaps)
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
           main=paste("Overlap for alpha:",round(x$call$alphas[i],3)))

      par(fg="orange")
      polygon(polyX,polyY,angle=90,density=30,border=NA)
      par(fg="black")

      lines(1:nn,mu, col="orange", t="s")
      lines(1:nn,x$overlaps[1:nn], col="blue", t="s")
      lines(1:nn,x$revOverlaps[1:nn], col="green", t="s")
      
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

        lines((nn+1):(2*nn),rev(mu), col="orange", t="s")
        lines((nn+1):(2*nn),x$overlaps[(NN-nn+1):NN], col="blue", t="s")
        lines((nn+1):(2*nn),x$revOverlaps[(NN-nn+1):NN], col="green", t="s")
        
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
      legend(1, ylims[2], legend=c("direct", "reversed", "expected"), col=c("blue", "green", "orange"),
             lty=1, bty="n")
      readline("<hit enter for next alpha>")
    }
    return(invisible())
  }
  if (which == "density") {
    plot(x$randomScores, x$scores, x$revScores)
    return(invisible())
  }
  stop("Value for parameter 'which' not recognized: ", which)
}


