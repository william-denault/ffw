#'@title Analytical pvalue for loci with only one Bayes Factor over a threshold
#'@title Analytical pvalue for loci with only one Bayes Factor over a threshold
#'@description  Compute the analytical pvalue for loci with only one Bayes Factor over a threshold
#'@param nsimu number of simulations to perform.
#'@param lambda the lambda parameter of the null distribution of the Bayes Factor.
#'@param lev_res the maximum level of resolution needed.
#'@param ncp the lambda parameter of the Bayes Factor's null distribution, if not specified set at 0.
#'@details Use the theoretical null distribution of the generated Bayes factor to perform simulation of null distribution of the Wavelet screaming procedure's test statistics.


analytical_p <- function(sub,lambda,ncp,thresh )
{
  if(missing(thresh))
  {
    thresh <-1
  }
  if(missing(ncp))
  {
    ncp <-0
  }


  if(length( which( as.numeric(sub)>thresh ) )>1 )
  {
    print( "More than one Bayes Factor above the choosen threshold (set as one per deault), no analytical formula available. Please use the function adaptative_Simu_Lambda_null")
    break
  }

  BF <- as.numeric(sub[which( as.numeric(sub)>thresh ) ])

  quant  <- (2*log(BF)-log(1 - lambda))/lambda


  pval <- pchisq(quant, df=1, ncp =  ncp, lower.tail = FALSE)

  return(pval)

}
