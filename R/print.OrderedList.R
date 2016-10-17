print.OrderedList <- function(x, ...){

  cat("Similarity of Ordered Gene Lists")
  cat("\n Comparison          :",x$label)
  cat("\n Number of genes     :", x$n)
  cat("\n Test statistic      :",x$call$test)
  cat("\n Number of subsamples:",x$call$B)
  if (x$call$beta==1){
    cat("\n beta = 1 -> corresponding labels could be matched in different studies")
  }
  if (x$call$beta==0.5){
    cat("\n beta = 0.5 -> corresponding labels could not be matched in different studies")
  }

  cat("\n--------------------------------------")
  cat("\n Optimal regularization parameter: alpha =",x$alpha)
  cat("\n Lists are more alike in", 
      ifelse(x$direction > 0, "direct", "reversed"), "order")
  cat("\n Weighted overlap score: ", x$sim.scores$SIM.obs)
  cat("\n Significance of similarity: p-value =",x$p)  
  cat("\n Number of genes contributing",x$call$percent*100,"% to similarity score:",length(x$intersect))
  cat("\n")
}
