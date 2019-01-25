library(MASS)

generate_random_B <- function(p = 10, a = 0.3, b = 0.7, m = 2, z = 0.05, s = 0.5){
  set.seed(12345)
  D = rep(1,p)
  plower = p*(p-1)/2
  ## off-diagonals
  T = diag(p)
  T[upper.tri(T)] = 0
  T[lower.tri(T)] = (ifelse(runif(plower)<s, -1, 1) *
                       ifelse(runif(plower)<z,  1, 0) *
                       runif(plower, a, b))
  L = diag(1.0/sqrt(D)) %*% T   # cholesky factor
  if(is.null(m)||(m==1)){
    Border = sample(1:p)
  } else {
    M = floor(p/m)
    Border = c()
    pos = 1
    while((pos+M+1)<p){
      Border = append(Border,sample(pos:(pos+M-1),M))
      pos = pos+M
    }
    Border = append(Border,sample(pos:p,p-pos+1))
  }
  B = L[Border,Border]
  adj = (B!=0)
  return(list(adj = adj,B = B))
}

#' Generates random multivariate normal data from a random B or provided adjacency
#'
#' @param n number of samples to generate
#' @param p number of dimesnions of data
#' @param adj adjacency matrix
#'
#' @return graph adjacency B
#'
#' @examples
#' generate_data(100,20)
#'
#' @export
generate_data <- function(n,p,B=NULL){
  if(is.null(B)){B = generate_random_B()$B}
  p = dim(B)[1]
  omega = t(B) %*% B            # omega
  sigma = solve(omega)
  sigma[abs(sigma)<1e-10] = 0   # set numerical error to zero
  set.seed(3689) #seed for generating data
  X = mvrnorm(n, mu=rep(0, p), Sigma=sigma) # observations
  X = scale(X,center = TRUE, scale=TRUE)
  return(X)
}
