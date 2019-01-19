## this is the block where we implement

require(Rcpp)
require('pdagDFS')

diagonalBLock <- function(B,S,l, m1 = 0, m2, p){
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
