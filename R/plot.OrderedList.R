plot.OrderedList <- function(x, which=NULL, no.title=FALSE, ...){

  y <- x$sim.score
  if (no.title) par(mai=par("mai")[c(1,2,4,4)])

  if (!is.null(which)){
    if (which=="scores"){
      rr <- range(y$SIM.alternative,y$SIM.random,y$SIM.observed)
      dalt <- density(y$SIM.alternative)
      dnul <- density(y$SIM.random)
      ylims <- range(c(dalt$y, dnul$y))
      tmp <- range(dnul$y)
      tmp <- ylims[1] + (tmp[2]-tmp[1])*1.2
      ylims[2] <- max(ylims[2], tmp)
      plot(dalt,col="red",xlim=rr,ylim=ylims,xlab="Similarity score",main="",...)
      lines(dnul)
      abline(v=y$SIM.observed,col="red")
      mtext("observed",side=3,at=y$SIM.observed)
      rug(y$SIM.alternative,col="red")
      rug(y$SIM.random)
      box()
      legend(x=min(rr),y=ylims[2], xjust=0,yjust=1, legend=c("resampled","random"),
             col=c("red","black"),lty=1, bty="n")
      text(max(rr), ylims[2], adj=c(1,1),
           label=paste(paste("direction:", ifelse(x$direction<0,"reversed","direct")),
                       paste("p-value:", round(x$p,3)),
                       paste("alpha:", round(x$alpha,3)),sep="\n"))
    }  
    
    if (which=="overlap"){
      order1 <- rownames(x$scores)[order(x$scores[,1],decreasing=TRUE)]
      order2 <- rownames(x$scores)[order(x$scores[,2],decreasing=TRUE)]
      if (x$direction < 0) {
        order2 <- rev(order2)
        tmp <- x$call$labels2[1]
        x$call$labels2[1] <- x$call$labels2[2]
        x$call$labels2[2] <- tmp
      }
      nn <- ceiling(-log(x$call$min.weight)/x$alpha)
      N <- x$n
      up <- overlap(order1, order2, nn)
      dn <- overlap(rev(order1), rev(order2), nn)
      overlapnum <- c(up, rev(dn))
      ylims <- range(overlapnum)

      if (!is.null(x$empirical)){
        median <- x$empirical$top[,2]           
        lower  <- x$empirical$top[,1]           
        upper  <- x$empirical$top[,3]           
      }
      else {
        median <- qhyper(0.5,1:nn,N-(1:nn),1:nn)
        upper  <- qhyper(0.975,1:nn,N-(1:nn),1:nn)
        lower  <- qhyper(0.025,1:nn,N-(1:nn),1:nn)
      }


      plot(0,xaxt="n",type="n",xlab="",ylab="size of overlap",xlim=c(1,2*nn), 
           ylim=ylims)

      polyX <- c(1:nn,nn:1)
      polyY <- c(upper,rev(lower))
      polyX <- rep(polyX,rep(2,length(polyX)))
      polyY <- rep(polyY,rep(2,length(polyY)))
      polyX <- polyX[-1]
      polyY <- polyY[-length(polyY)]
      
      par(fg="orange")
      polygon(polyX,polyY,angle=90,density=30,border=NA)
      lines(1:nn, median, col="orange", t="s")

      if (!is.null(x$empirical)){
        median <- x$empirical$bottom[,2]           
        lower  <- x$empirical$bottom[,1]           
        upper  <- x$empirical$bottom[,3]           
      }

      polyX <- c((nn+1):(2*nn),rev((nn+1):(2*nn)))
      polyY <- c(rev(upper),lower)
      polyX <- rep(polyX,rep(2,length(polyX)))
      polyY <- rep(polyY,rep(2,length(polyY)))
      polyX <- polyX[-1]
      polyY <- polyY[-length(polyY)]

      polygon(polyX,polyY,angle=90,density=30,border=NA)
      lines((nn+1):(2*nn),rev(median), col="orange", t="s")
      par(fg="black")

      lines(1:nn,overlapnum[1:nn], t="s")
      lines((nn+1):(2*nn),overlapnum[(nn+1):(2*nn)], t="s")
      
      axStep <- floor(nn / 45) * 10
      axis(side=1, at=c(axStep*(0:4), 2*nn-axStep*(0:4)), 
           labels=c(1, axStep*(1:4), 1, axStep*(1:4)))
      mtext("top ranks", side=1, at=0.5*nn, line=2.5)
      mtext("bottom ranks", side=1, at=1.5*nn, line=2.5)
      abline(v=nn+0.5,lty=2)
      text(1, max(overlapnum), adj=c(0,1),
           label=paste("Upregulated:",
                       paste(x$call$labels1[1], "in", colnames(x$scores)[1]),
                       paste(x$call$labels2[1], "in", colnames(x$scores)[2]),sep="\n"))
      text(2*nn, max(overlapnum), adj=c(1,1),
           label=paste("Upregulated:",
                       paste(x$call$labels1[2], "in", colnames(x$scores)[1]),
                       paste(x$call$labels2[2], "in", colnames(x$scores)[2]),sep="\n"))
      legend(nn+1, ylims[1], lty=1, bty="n", yjust=0, 
             legend=c("observed", "expected"), col=c("black", "orange"))
 
    }
    
    if (which=="pauc"){
      ylims <- range(x$pauc)
      alphas <- as.numeric(colnames(x$pauc))
      plot(alphas, x$pauc[1,], type="l", col="blue", xlim=range(alphas), ylim=ylims,
           main=, xlab=expression(alpha),ylab="pAUC score", log="x")
      abline(v=x$alpha,col="red")
      mtext(expression(alpha[opt]), side=3, at=x$alpha)
      if (nrow(x$pauc) == 2) {
        lines(alphas,x$pauc[2,], col="green")
        legend(x=min(alphas), y=ylims[2], xjust=0,yjust=1, legend=c("direct","reversed"),
               col=c("blue","green"),lty=1, bty="n")
      
      }
    }
    if (!no.title) title(main=paste("Comparison: ", x$label,sep=""))
  }
  
  if (is.null(which)){
#    X11(width=8.1,height=2.5)
#    par(mfrow=c(1,3), ps=6)
    par(mfrow=c(1,3))
    plot(x,"scores",no.title=TRUE)
    plot(x,"overlap",no.title=TRUE)
    plot(x,"pauc",no.title=TRUE)
  }

  return(invisible())
}
