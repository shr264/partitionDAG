## this is the block where we implement

require(Rcpp)
require('pdagDFS')

#' Estimates an adjacencey matrix for a DAG based on l1 penalized negative likelihood minimization given a partitioning of the nodes into two groups
#'
#' @param B a matrix of size p by p 
#' @param S a matrix of size p by p 
#' @param l penalization parameter
#' @param m1 the node at which the partition occurs
#' @param m2 not relevant
#' @param p maximum number of variables
#'
#' @return graph adjacency B
#'
#' @examples
#' diagonalBLock(B,S,l, m1 = 0, m2, p)
#'
#' @export
diagonalBLock <- function(B,S,l, m1 = 0, m2, p){
  cat('Processing diagonal block!\n')
  if(m2 == m1 + 1){
    i = m2
    km1_ = 1:p
    km1_ = km1_[!km1_ == i]
    sumterm = sum(B[km1_,i] * S[km1_,i])
    B[i,i] = (-sumterm + sqrt(1*sumterm^2 + 4*S[i, i]))/(2*S[i, i])
    return(B)
    }
  else{
    for(i in (m1+1):m2){
      km1_ = 1:p
      km1_ = km1_[!km1_ == i]
      sumterm = sum(B[km1_,i] * S[km1_,i])
      B[i,i] = (-sumterm + sqrt(1*sumterm^2 + 4*S[i, i]))/(2*S[i, i])
    }
    for(i in (m1+2):m2){
      for(k in (m1+1):(i-1)){
        Bupdateki = Bki(S,B,k,i,l)
        Bupdateik = Bik(S,B,i,k,l)
        hascycleki = 0
        hascycleik = 0
        if(abs(Bupdateki)>0){
          graphki = abs(B)>0
          diag(graphki) = 0
          graphki[k,i] = 1
          if(dfs2(graphki,p)==1){
            hascycleki = 1
          }
        }
  
        if(abs(Bupdateik)>0){
          graphik = abs(B)>0
          diag(graphik) = 0
          graphik[i,k] = 1
          if(dfs2(graphik,p)==1){
            hascycleik = 1
          }
        }
        if(hascycleik==1){
          Bupdateik = 0.0
        } else if(hascycleki==1){
          Bupdateki = 0.0
        } else {
          qik = qBik(S,Bupdateik,B,i,k,l)
          qki = qBki(S,Bupdateki,B,k,i,l)
          if(qik$q>qki$q){
            Bupdateki = qki$b
            Bupdateik = 0
          } else if (qik$q<qki$q) {
            Bupdateik = qik$b
            Bupdateki = 0
          } else {
            u = runif(n=1)
            if(u<0.5){
              Bupdateki = qki$b
              Bupdateik = 0
            } else {
              Bupdateik = qik$b
              Bupdateki = 0
            }
          }
        }
        B[i,k] = Bupdateik
        B[k,i] = Bupdateki
      }
    }
    return(B)
  }
}

#' Estimates an adjacencey matrix for a DAG based on l1 penalized negative likelihood minimization given a partitioning of the nodes into two groups
#'
#' @param B a matrix of size p by p 
#' @param S a matrix of size p by p 
#' @param l penalization parameter
#' @param m1 the node at which the partition occurs
#' @param m2 not relevant
#' @param p maximum number of variables
#'
#' @return graph adjacency B
#'
#' @examples
#' diagonalBLock_ll(B,S,l, m1 = 0, m2, p)
#'
#' @export
diagonalBLock_ll <- function(B,S,l, m1 = 0, m2, p){
  cat('Processing diagonal block!\n')
  for(i in (m1+2):m2){
    for(k in (m1+1):(i-1)){
      Bupdateki = Bki(S,B,k,i,l)
      Bupdateik = Bik(S,B,i,k,l)
      hascycleki = 0
      hascycleik = 0
      if(abs(Bupdateki)>0){
        graphki = abs(B)>0
        diag(graphki) = 0
        graphki[k,i] = 1
        if(dfs2(graphki,p)==1){
          hascycleki = 1
        }
      }

      if(abs(Bupdateik)>0){
        graphik = abs(B)>0
        diag(graphik) = 0
        graphik[i,k] = 1
        if(dfs2(graphik,p)==1){
          hascycleik = 1
        }
      }
      if(hascycleik==1){
        Bupdateik = 0.0
      } else if(hascycleki==1){
        Bupdateki = 0.0
      } else {
        qik = qBik_ll(S,Bupdateik,B,i,k,l)
        qki = qBki_ll(S,Bupdateki,B,k,i,l)
        if(qik$q>qki$q){
          Bupdateki = qki$b
          Bupdateik = 0
        } else if (qik$q<qki$q) {
          Bupdateik = qik$b
          Bupdateki = 0
        } else {
          u = runif(n=1)
          if(u<0.5){
            Bupdateki = qki$b
            Bupdateik = 0
          } else {
            Bupdateik = qik$b
            Bupdateki = 0
          }
        }
      }
      B[i,k] = Bupdateik
      B[k,i] = Bupdateki
    }
  }
  return(B)
}
