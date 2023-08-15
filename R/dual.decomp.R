#' This function converts the principal component decomposition to the dual
#' supervised form
#' @param X nxp matrix of predictors
#' @param L nxn response kernel
#' @param npc number of desired principal components
#' @param n an integer storing the number response observations
#' @param p an integer storing the number of predictors 
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

dual.decomp <- function(X, L, npc, n, p){
  
  Psi <- get.psi(X=X,L=L,n=n,p=p)
  # Generate Dual supervised matrix  
  #   This matrix is nxn instead of pxp
  M <- Rfast::Tcrossprod(Psi, Psi)
  # Decomposition of Dual supervised form
  Md <- eigen(M)
  if(npc==1){
    Sigma_Meval <- sqrt(Md$values[1])
    U <- Rfast::mat.mult(solve(Sigma_Meval),Rfast::Crossprod(as.matrix(Md$vectors)[,1], Psi))
    #U <- t(t(Psi)%*%Md$vectors%*%t(solve(Sigma_Meval)))
    return(list(vectors=t(U), values=Md$values[1]))
  }
  else{
    Sigma_Meval <- diag(sqrt(Md$values[1:npc]))
    #Ut <- t(t(Psi)%*%Uvectors%*%t(solve(Siggy)))
    U <- Rfast::mat.mult(solve(Sigma_Meval),Rfast::Crossprod(as.matrix(Md$vectors)[,1:npc], Psi))
    return(list(vectors=t(U), values=Md$values[1:npc]))
  }
}