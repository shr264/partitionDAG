softhresh <- function(x,l){
  sign(x)*max(abs(x)-l,0)
}

Bii <- function(S,B,i){
  temp = sum(S[i,-i]*B[i,-i])
  out = (-2*temp + sqrt(4*temp^2+8*S[i,i]))/(4*S[i,i])
  out
}

Bki <- function(S,B,k,i,l){
  out = softhresh(-2*sum(S[i,-i]*B[k,-i])/(4*S[i,i]),l/(4*S[i,i]))
  out
}

qBki <- function(S,Bki,B,k,i,l){
  out1 = S[i,i]*Bki^2 + 2*Bki*sum(S[i,-i]*B[k,-i]) + l*abs(Bki)
  out2 = 0
  if(out1<out2){return(list(q = out1, b = Bki))
  } else {return(list(q = out2, b = 0))}
}

penloglik <- function(B,S,l){
  Omega = t(B)%*%B
  out = sum(diag(Omega%*%S)) - log(det(Omega)) + l*sum(abs(B))
  out
}

qBki_ll <- function(S,Bki,B,k,i,l){
  B[k,i] = 0
  old_ll = penloglik(B,S,l)
  B[k,i] = Bki
  B[i,k] = 0
  new_ll = penloglik(B,S,l)
  if(new_ll<old_ll){return(list(q = new_ll, b = Bki))
  } else {return(list(q = old_ll, b = 0))}
}

Bik <- function(S,B,i,k,l){
  out = softhresh(-2*sum(S[k,-k]*B[i,-k])/(4*S[k,k]),l/(4*S[k,k]))
  out
}

qBik <- function(S,Bik,B,i,k,l){
  out1 = S[k,k]*Bik^2 + 2*Bik*sum(S[k,-k]*B[i,-k]) + l*abs(Bik)
  out2 = 0
  if(out1<out2){return(list(q = out1, b = Bik))
  } else {return(list(q = out2, b = 0))}
}

qBik_ll <- function(S,Bik,B,i,k,l){
  B[i,k] = 0
  old_ll = penloglik(B,S,l)
  B[i,k] = Bik
  B[k,i] = 0
  new_ll = penloglik(B,S,l)
  if(new_ll<old_ll){return(list(q = new_ll, b = Bik))
  } else {return(list(q = old_ll, b = 0))}
}
