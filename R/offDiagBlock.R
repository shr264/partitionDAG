
offDiagBlock <- function(B,S,l,m1,m2,p){
  cat('Processing OFF-diagonal block!')
  for(i in (m1+1):m2){
    for(k in 1:m1){
      B[k,i] = 0
      B[i,k] = Bik(S,B,i,k,l)
    }
  }
  return(B)
}
