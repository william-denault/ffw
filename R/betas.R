#'@title Compute Beta value from Bayesian linear regression
#'@description  Compute Beta value from Bayesian linear regression
#'@param y phenotype vector/variable of interest has to be numeric.
#'@param x numerical, variable to regress on.
#'@param confounder the confounding matrix with the same sample order as Y. The intercept should not be included. If missing will generate an intercept matrix.
#'@param sigma_b the parameter of the NIG prior used for the Bayes Factor computation. We advised setting this value between 0.1 and 0.2
#'@param all logical, if set as TRUE return all the Beta value (including the ones from the confounding factors). If set as FALSE, only return the estimate for x, set as FALSE if missing.
#'@details The Wavelet_screaming function computes the Likelihood ratio used for testing the significance of a genetic region. In addition, it computes
#'the proportion of wavelets coefficients associated by the level of resolution, and the Bayes factor used for this estimation. All the details
#'of the computation can be found in our paper "Detecting differentially methylated regions using a fast wavelet-based approach to functional association analysis" Denault and Jugessur
#'@return A named vector. First position, the estimated value of the Lambda statistics, then the proportion of association per level of resolution then the computed Bayes Factor per wavelet coefficient.
#'@examples \dontrun{
#'
#'x <- rnorm(1000)
#'y <- x+3
#'sigma_b <- 0.2
#'betas(y=y,x=x,sigma_b = sigma_b,all=TRUE)
#'}




betas <- function(y,x,confounder, sigma_b,all)
{
  if(missing(confounder) )
  {
    confounder <- data.frame( confounder=rep(1,length(y)) )
  }
  else{
    confounder <- data.frame(confounder)
  }
  if(missing(all))
  {
    all=FALSE
  }

  pc <- dim(confounder)[2]
  Dmat <- cbind(confounder,x)
  Dmat <- as.matrix(Dmat)
  index <- pc+1
  betas <- solve(t(Dmat) %*% Dmat + diag(1/sigma_b/sigma_b,index)) %*% t(Dmat)%*% y
  colnames(betas) <- "Beta"
  #postvar <- diag( solve(t(Dmat) %*% Dmat + diag(1/sigma_b/sigma_b,2)))*1/sigma_b/sigma_b
  if(all==FALSE)
  {
    return(betas[index,1])
  }
  if(all==TRUE)
  {
    return(betas)
  }

}





