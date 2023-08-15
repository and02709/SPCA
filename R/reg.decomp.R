#' This function decomposes the learned matrix
#' 
#' @param X nxp matrix of predictors
#' @param L nxn response kernel
#' @param npc number of desired principal components
#' @param n an integer storing the number response observations
#' @param p an integer storing the number of predictors 
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

reg.decomp <- function(X, L, npc, n, p){
  H <- diag(x=1, nrow=n) - 1/n*rep(1,n)%*%t(rep(1,n))
  M <- Rfast::Crossprod(X, Rfast::mat.mult(H, Rfast::mat.mult(L, Rfast::mat.mult(H, X))))
  Md <- eigen(M)
  return(list(vectors=Md$vectors[1:npc,], values=Md$values[1:npc]))
}