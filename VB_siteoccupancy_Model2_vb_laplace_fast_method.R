vb_la<-function(W, X, y, alpha_0, beta_0, Sigma_alpha_0, Sigma_beta_0, LargeSample=TRUE, epsilon=1e-5) {UseMethod("vb_la")}

vb_la.default<-function(W, X, y, alpha_0, beta_0, Sigma_alpha_0, Sigma_beta_0, LargeSample=TRUE, epsilon=1e-5)
{
  
  out<-vb_model2_la(W, X, y, alpha_0, beta_0, Sigma_alpha_0, Sigma_beta_0, LargeSample, epsilon)
  out$alpha <- as.vector(out$alpha)
  out$beta <- as.vector(out$beta)

  out$call <- match.call()
  class(out) <- "vb_la"
  out
}
  

print.vb_la <- function(out)
{
  cat("Call:\n")
  print(out$call)
  cat("\nOccurrence Coefficients:\n")
  cat("------------------------ \n")
  print(out$alpha)
  cat("\nDetection Coefficients:\n")
  cat("----------------------- \n")
  print(out$beta)
}



summary.vb_la <- function(object, ...)
{
  se1 <- sqrt(diag(object$Sigma_alpha))
  se2 <- sqrt(diag(object$Sigma_beta))
  
  TAB <- cbind(Estimate=c(object$alpha,object$beta), StdErr=c(se1,se2))
  rownames(TAB) <-NULL
  
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.vb_la"
  res
}


print.summary.vb_la <- function(x)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$coefficients)
}

#save(y, file="y.rda")
#save(X, file="X.rda")
#save(W, file="W.rda")

#vb1<-vb_la(W=W, X=X, y=y, alpha_0=alpha_0, beta_0=beta_0, Sigma_alpha_0=Sigma_alpha_0, Sigma_beta_0=Sigma_beta_0,LargeSample=TRUE, epsilon=1e-5)
