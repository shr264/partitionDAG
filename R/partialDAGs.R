
ccdr_custom <- function(X,l,a=NULL,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,eps = 10^(-4),maxitr = 100, init=NULL){
    (n = dim(X)[1])
    (p = dim(X)[2])
    (S = (1/n)*t(X)%*%X)
    if(is.null(init)){
      B11 = diag(p)
      #+ matrix(rnorm(p^2,0,0.01)*ifelse(runif(p^2)<0.3,  1, 0), nrow = p, ncol = p)
      B = Bold = B11
    } else {B = Bold = init}
    diff = 1
    itr = 1
    ptm <- rep(0,maxitr)
    while((diff>eps)&(itr<maxitr)){
      ptm[itr] = proc.time()[3]
      cat('itr:', itr, '...')
      ## first all the diagonals
      for(i in 1:p){
        B[i,i] = Bii(S,B,i)
      }

      B = diagonalBLock(B,S,l,0,p,p)
      ptm[itr] = proc.time()[3] - ptm[itr]
      diff = max(B-Bold)
      Bold = B

      itr = itr +1
    }
    time <- sum(ptm)
    return(list(B=B,itr = itr, time = time))
}

ccdr_custom_ll <- function(X,l,a=NULL,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  (S = (1/n)*t(X)%*%X)
  if(is.null(init)){
    B11 = diag(p)
    #+ matrix(rnorm(p^2,0,0.01)*ifelse(runif(p^2)<0.3,  1, 0), nrow = p, ncol = p)
    B = Bold = B11
  } else {B = Bold = init}
  diff = 1
  itr = 1
  ptm <- rep(0,maxitr)
  while((diff>eps)&(itr<maxitr)){
    ptm[itr] = proc.time()[3]
    cat('itr:', itr, '...')
    ## first all the diagonals
    for(i in 1:p){
      B[i,i] = Bii(S,B,i)
    }

    B = diagonalBLock_ll(B,S,l,0,p,p)
    ptm[itr] = proc.time()[3] - ptm[itr]
    diff = max(B-Bold)
    Bold = B

    itr = itr +1
  }
  time <- sum(ptm)
  return(list(B=B,itr = itr, time = time))
}

partial2 <- function(X,l,a=NULL,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  (S = (1/n)*t(X)%*%X)
  if(is.null(init)){
    B11 = diag(m1)
    B22 = diag(p-m1)
    B = Bold = as.matrix(bdiag(B11,B22))
  } else {B = Bold = init}
  diff = 1
  itr = 1
  ptm <- rep(0,maxitr)
  while((diff>eps)&(itr<maxitr)){
    time = proc.time()[3]
    cat('itr:', itr, '...')
    ## first all the diagonals
    for(i in 1:p){
      B[i,i] = Bii(S,B,i)
    }
    time = proc.time()[3] - time

    B = diagonalBLock(B,S,l,0,m1,p)

    time2 = proc.time()[3]
    B = offDiagBlock(B,S,l,m1,p,p)

    B = diagonalBLock(B,S,l,m1,p,p)
    time2 = proc.time()[3] - time2

    diff = max(B-Bold)
    Bold = B
    ptm[itr] = time + time2
    itr = itr +1
  }
  time <- sum(ptm)
  return(list(B=B,itr = itr, time = time))
}

partial3 <- function(X,l,a=NULL,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  (S = (1/n)*t(X)%*%X)
  if(is.null(init)){
    B11 = diag(m1)
    B22 = diag(m2-m1)
    B33 = diag(p-m2)
    B = Bold = as.matrix(bdiag(B11,B22,B33))
  } else {B = Bold = init}
  diff = 1
  itr = 1
  ptm <- rep(0,maxitr)
  while((diff>eps)&(itr<maxitr)){
    time = proc.time()[3]
    cat('itr:', itr, '...')
    ## first all the diagonals
    for(i in 1:p){
      B[i,i] = Bii(S,B,i)
    }
    time = proc.time()[3] - time

    B = diagonalBLock(B,S,l,0,m1,p)

    B = offDiagBlock(B,S,l,m1,m2,p)

    B = diagonalBLock(B,S,l,m1,m2,p)

    time2 = proc.time()[3]
    B = offDiagBlock(B,S,l,m2,p,p)

    B = diagonalBLock(B,S,l,m2,p,p)
    time2 = proc.time()[3] - time2

    diff = max(B-Bold)
    Bold = B
    ptm[itr] = time + time2
    itr = itr +1
  }
  time <- sum(ptm)
  return(list(B=B,itr = itr, time = time))
}

partial4 <- function(X,l,a=NULL,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  (S = (1/n)*t(X)%*%X)
  if(is.null(init)){
    B11 = diag(m1)
    B22 = diag(m2-m1)
    B33 = diag(p-m2)
    B = Bold = as.matrix(bdiag(B11,B22,B33))
  } else {B = Bold = init}
  diff = 1
  itr = 1
  ptm <- rep(0,maxitr)
  while((diff>eps)&(itr<maxitr)){
    time = proc.time()[3]
    cat('itr:', itr, '...')
    ## first all the diagonals
    for(i in 1:p){
      B[i,i] = Bii(S,B,i)
    }
    time = proc.time()[3] - time

    B = diagonalBLock(B,S,l,0,m1,p)

    B = offDiagBlock(B,S,l,m1,m2,p)

    B = diagonalBLock(B,S,l,m1,m2,p)

    B = offDiagBlock(B,S,l,m2,m3,p)

    B = diagonalBLock(B,S,l,m2,m3,p)

    time2 = proc.time()[3]
    B = offDiagBlock(B,S,l,m3,p,p)

    B = diagonalBLock(B,S,l,m3,p,p)
    time2 = proc.time()[3] - time2

    diff = max(B-Bold)
    Bold = B
    ptm[itr] = time + time2
    itr = itr +1
  }
  time <- sum(ptm)
  return(list(B=B,itr = itr, time = time))
}

partial5 <- function(X,l,a=NULL,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  (S = (1/n)*t(X)%*%X)
  if(is.null(init)){
    B11 = diag(m1)
    B22 = diag(m2-m1)
    B33 = diag(p-m2)
    B = Bold = as.matrix(bdiag(B11,B22,B33))
  } else {B = Bold = init}
  diff = 1
  itr = 1
  while((diff>eps)&(itr<maxitr)){
    cat('itr:', itr, '...')
    ## first all the diagonals
    for(i in 1:p){
      B[i,i] = Bii(S,B,i)
    }

    B = diagonalBLock(B,S,l,0,m1,p)

    B = offDiagBlock(B,S,l,m1,m2,p)

    B = diagonalBLock(B,S,l,m1,m2,p)

    B = offDiagBlock(B,S,l,m2,m3,p)

    B = diagonalBLock(B,S,l,m2,m3,p)

    B = offDiagBlock(B,S,l,m3,m4,p)

    B = diagonalBLock(B,S,l,m3,m4,p)

    B = offDiagBlock(B,S,l,m4,p,p)

    B = diagonalBLock(B,S,l,m4,p,p)

    diff = max(B-Bold)
    Bold = B
    itr = itr +1
  }
  return(list(B=B,itr = itr))
}

partial8 <- function(X,l,a=NULL,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  (S = (1/n)*t(X)%*%X)
  if(is.null(init)){
    B11 = diag(m1)
    B22 = diag(m2-m1)
    B33 = diag(p-m2)
    B = Bold = as.matrix(bdiag(B11,B22,B33))
  } else {B = Bold = init}
  diff = 1
  itr = 1
  while((diff>eps)&(itr<maxitr)){
    cat('itr:', itr, '...')
    ## first all the diagonals
    for(i in 1:p){
      B[i,i] = Bii(S,B,i)
    }

    B = diagonalBLock(B,S,1*l,0,m1,p)

    B = offDiagBlock(B,S,1*l,m1,m2,p)

    B = diagonalBLock(B,S,1*l,m1,m2,p)

    B = offDiagBlock(B,S,1*l,m2,m3,p)

    B = diagonalBLock(B,S,1*l,m2,m3,p)

    B = offDiagBlock(B,S,1*l,m3,m4,p)

    B = diagonalBLock(B,S,1*l,m3,m4,p)

    B = offDiagBlock(B,S,1*l,m4,m5,p)

    B = diagonalBLock(B,S,1*l,m4,m5,p)

    B = offDiagBlock(B,S,1*l,m5,m6,p)

    B = diagonalBLock(B,S,1*l,m5,m6,p)

    B = offDiagBlock(B,S,1*l,m6,m7,p)

    B = diagonalBLock(B,S,1*l,m6,m7,p)

    B = offDiagBlock(B,S,1*l,m7,p,p)

    B = diagonalBLock(B,S,1*l,m7,p,p)

    diff = max(B-Bold)
    Bold = B
    itr = itr +1
  }
  return(list(B=B,itr = itr))
}

partial10 <- function(X,l,a=NULL,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  (S = (1/n)*t(X)%*%X)
  if(is.null(init)){
    B11 = diag(m1)
    B22 = diag(m2-m1)
    B33 = diag(p-m2)
    B = Bold = as.matrix(bdiag(B11,B22,B33))
  } else {B = Bold = init}
  diff = 1
  itr = 1
  while((diff>eps)&(itr<maxitr)){
    cat('itr:', itr, '...')
    ## first all the diagonals
    for(i in 1:p){
      B[i,i] = Bii(S,B,i)
    }

    B = diagonalBLock(B,S,1*l,0,m1,p)

    B = offDiagBlock(B,S,1*l,m1,m2,p)

    B = diagonalBLock(B,S,1*l,m1,m2,p)

    B = offDiagBlock(B,S,1*l,m2,m3,p)

    B = diagonalBLock(B,S,1*l,m2,m3,p)

    B = offDiagBlock(B,S,1*l,m3,m4,p)

    B = diagonalBLock(B,S,1*l,m3,m4,p)

    B = offDiagBlock(B,S,1*l,m4,m5,p)

    B = diagonalBLock(B,S,1*l,m4,m5,p)

    B = offDiagBlock(B,S,1*l,m5,m6,p)

    B = diagonalBLock(B,S,1*l,m5,m6,p)

    B = offDiagBlock(B,S,1*l,m6,m7,p)

    B = diagonalBLock(B,S,1*l,m6,m7,p)

    B = offDiagBlock(B,S,1*l,m7,m8,p)

    B = diagonalBLock(B,S,1*l,m7,m8,p)

    B = offDiagBlock(B,S,1*l,m8,m9,p)

    B = diagonalBLock(B,S,1*l,m8,m9,p)

    B = offDiagBlock(B,S,1*l,m9,p,p)

    B = diagonalBLock(B,S,1*l,m9,p,p)

    diff = max(B-Bold)
    Bold = B
    itr = itr +1
  }
  return(list(B=B,itr = itr))
}

partial8_inc <- function(X , l, m1, m2, m3, m4, m5, m6, m7, y = 0, z = 0,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  (S = (1/n)*t(X)%*%X)
  if(is.null(init)){
    B11 = diag(m1)
    B22 = diag(m2-m1)
    B33 = diag(p-m2)
    B = Bold = as.matrix(bdiag(B11,B22,B33))
  } else {B = Bold = init}
  diff = 1
  itr = 1
  while((diff>eps)&(itr<maxitr)){
    cat('itr:', itr, '...')
    ## first all the diagonals
    for(i in 1:p){
      B[i,i] = Bii(S,B,i)
    }

    B = diagonalBLock(B,S,(1+1*y)^(1*z)*l,0,m1,p)

    B = offDiagBlock(B,S,(1+1*y)^(1*z)*l,m1,m2,p)

    B = diagonalBLock(B,S,(1+2*y)^(2*z)*l,m1,m2,p)

    B = offDiagBlock(B,S,(1+2*y)^(2*z)*l,m2,m3,p)

    B = diagonalBLock(B,S,(1+3*y)^(3*z)*l,m2,m3,p)

    B = offDiagBlock(B,S,(1+3*y)^(3*z)*l,m3,m4,p)

    B = diagonalBLock(B,S,(1+4*y)^(4*z)*l,m3,m4,p)

    B = offDiagBlock(B,S,(1+4*y)^(4*z)*l,m4,m5,p)

    B = diagonalBLock(B,S,(1+5*y)^(5*z)*l,m4,m5,p)

    B = offDiagBlock(B,S,(1+5*y)^(5*z)*l,m5,m6,p)

    B = diagonalBLock(B,S,(1+6*y)^(6*z)*l,m5,m6,p)

    B = offDiagBlock(B,S,(1+6*y)^(6*z)*l,m6,m7,p)

    B = diagonalBLock(B,S,(1+7*y)^(7*z)*l,m6,m7,p)

    B = offDiagBlock(B,S,(1+8*y)^(8*z)*l,m7,p,p)

    B = diagonalBLock(B,S,(1+8*y)^(8*z)*l,m7,p,p)

    diff = max(B-Bold)
    Bold = B
    itr = itr +1
  }
  return(list(B=B,itr = itr))
}

partial8_dec <- function(X,l,m1,m2,m3,m4,m5,m6,m7,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  (S = (1/n)*t(X)%*%X)
  if(is.null(init)){
    B11 = diag(m1)
    B22 = diag(m2-m1)
    B33 = diag(p-m2)
    B = Bold = as.matrix(bdiag(B11,B22,B33))
  } else {B = Bold = init}
  diff = 1
  itr = 1
  while((diff>eps)&(itr<maxitr)){
    cat('itr:', itr, '...')
    ## first all the diagonals
    for(i in 1:p){
      B[i,i] = Bii(S,B,i)
    }

    B = diagonalBLock(B,S,8*l,0,m1,p)

    B = offDiagBlock(B,S,7*l,m1,m2,p)

    B = diagonalBLock(B,S,7*l,m1,m2,p)

    B = offDiagBlock(B,S,6*l,m2,m3,p)

    B = diagonalBLock(B,S,6*l,m2,m3,p)

    B = offDiagBlock(B,S,5*l,m3,m4,p)

    B = diagonalBLock(B,S,5*l,m3,m4,p)

    B = offDiagBlock(B,S,4*l,m4,m5,p)

    B = diagonalBLock(B,S,4*l,m4,m5,p)

    B = offDiagBlock(B,S,3*l,m5,m6,p)

    B = diagonalBLock(B,S,3*l,m5,m6,p)

    B = offDiagBlock(B,S,2*l,m6,m7,p)

    B = diagonalBLock(B,S,2*l,m6,m7,p)

    B = offDiagBlock(B,S,1*l,m7,p,p)

    B = diagonalBLock(B,S,1*l,m7,p,p)

    diff = max(B-Bold)
    Bold = B
    itr = itr +1
  }
  return(list(B=B,itr = itr))
}
