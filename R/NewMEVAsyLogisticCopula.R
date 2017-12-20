#' Creates a multivariate asymmetric logistic copula
#'
#' Creates an instance of the multivariate asymmetric copula with parameters \eqn{\theta} and r.
#'
#' If \code{theta} has entries \eqn{\theta_{ij}} and \code{r} has entries \eqn{r_j} (\eqn{i=1,\dots,k} and \eqn{j=1,\dots,d}), then the following parameterisation of the copula is used:
#' \deqn{C(u_1,\dots,u_d) = \exp\left(- \sum_{i=1}^k \left\{ \sum_{j=1}^d (\theta_{ij} \bar u_j)^{r_i} \right\}^{1/r_i}  \right)}{C(u_1,\dots,u_d) = exp(- \sum_i { \sum_j (\theta_{ij} v_j)^{r_i} }^{1/r_i} )}
#' where \eqn{\bar u_j = -\log(u_j)}{v_j = -log(u_j)}, \eqn{j=1,\dots,d}.
#'
#' Necessary and sufficient conditions for the parameters are
#' \itemize{
#'   \item all entries in \code{theta} have to be non-negative.
#'   \item each column of \code{theta} has to add to one.
#'   \item each row of \code{theta} must have a unique pattern of non-zero values, including the pattern that has no zeros in a row.
#'   \item if a row of \code{theta} has only one non-zero value, then the corresponding entry in \code{r} has to be one.
#'   \item if a row of \code{theta} has more than one non-zero value, then the corresponding entry of \code{r} must be greater than one.
#' }
#'
#' @param theta a \eqn{k \times d}{k x d} matrix of reals.
#' @param r a vector of \eqn{k} reals
#'
#' @return A function that evaluates the multivariate asymmetric logistic copula (with parameters \eqn{\theta} and r) at a given \eqn{d}-dimensional vector in the unit square.  Note that for this function the dimension of vectors at which the copula can be evaluated is determined by the input parameters.  The environment of the function also contains a function called \code{pdfCopula} that evaluates the probability density function of the multivariate asymmetric logistic copula via automatic differentation.
#'
#' @seealso \code{\link{NewBEVAsyLogisticCopula}}, \code{\link{NewMEVGumbelCopula}}
#'
#' @author  Berwin A. Turlach <berwin.turlach@gmail.com>
#'
#' @keywords distribution
#'
#' @examples
#' theta <-  rbind(c(0, 0.2, 0.8), c(1,0.8,0.2))
#' r <- c(2,3)
#' cop <- NewMEVAsyLogisticCopula(theta, r)
#'
#' ## Creates the same copula
#' theta <- 0.2
#' phi <- 0.4
#' r <- 2
#' cop1 <- NewBEVAsyLogisticCopula(r, theta, phi)
#' theta <- cbind(c(phi, 1-phi, 0), c(theta, 0, 1-theta))
#' r <- c(r, 1,  1)
#' cop2 <- NewMEVAsyLogisticCopula(theta, r)
#'
#' @export
NewMEVAsyLogisticCopula <- function(theta, r){

  dim <- NCOL(theta)

  if( !all(theta >= 0 ) )
    stop("Entries in 'theta' are not non-negative.")
  if( !isTRUE(all.equal(rep(1,dim), colSums(theta), check.attributes=FALSE)) )
    stop("Columns of 'theta' do not add to one.")
  tmp <- theta != 0
  if( any(duplicated(tmp)) )
    stop("Repeated zero patterns in rows of 'theta'.")
  ii <- rowSums(tmp) == 1
  jj <- ! ii
  if( !isTRUE(all.equal(rep(1,sum(ii)), r[ii], check.attributes=FALSE)) )
    stop("Entries in 'r' corresponding to rows with only one non-zero entry in 'theta' are not one.")
  if( !all(r[jj] > 1) )
    stop("Redundancies in parameters, entries of one in 'r' corresponding to rows in 'theta' with more than one non-zero entry.")

  rm(tmp, ii, jj)

  param <- list(theta = theta, r = r)
  type <- c("asymmetric logistic", paste0(dim,"-variate"))

  pdfCopula <- function(x){
    if( length(x) != dim )
      stop(paste("x has the wrong length, should be ", dim, ".", sep=""))
    if( any(x<0) || any(x>1) )
      return(0)

    acc1 <- lx <- vector("list", dim)
    for(i in 1:dim){
      lx[[i]] <- NewCube(0, dim=dim)
      lx[[i]]$val[1] <- -log(x[i])
      lx[[i]]$val[2^(i-1)+1] <- -1/x[i]
      acc1[[i]] <- NewCube(0, dim=dim)
    }

    res <- NewCube(0, dim=dim)
    tmp <- NewCube(0, dim=dim)
    for(j in 1:NROW(theta)){
      for(i in 1:dim){
        acc1[[i]]$val <- lx[[i]]$val * theta[j,i]
        if( theta[j,i] != 0 && r[j] != 1 ){
          k <- 2^(i-1)+1
          acc1[[i]]$val[k] <- r[j]*acc1[[i]]$val[1]^(r[j]-1)*acc1[[i]]$val[k]
          acc1[[i]]$val[1] <- acc1[[i]]$val[1]^r[j]
        }
      }
      tmp$val <- rowSums(sapply(acc1, function(x) x$val))
      if( r[j] != 1 )
        tmp <- CrossPow(tmp, 1/r[j])
      res$val <- res$val - tmp$val
    }
    res <- CrossExp(res)
    if( !isTRUE(all.equal(res$val[1], cdfCopula(x), check.attributes=FALSE)) )
      stop("Inconsistency between implementation of cdf and pdf of this copula.")
    res$val[2^dim]
  }

  res <- cdfCopula <- function(x){
    if( length(x) != dim )
      stop(paste("x has the wrong length, should be ", dim, ".", sep=""))

    if( any(x==0) ) return(0)
    if( all(x==1) ) return(1)
    tmp <- log(x)
    sumtmp <- sum(tmp)
    tmp <- tmp/sumtmp

    tmp1 <- sweep(theta, 2, tmp, "*")^r
    tmp1 <- rowSums(tmp1)^(1/r)
    exp(sumtmp*sum(tmp1))
  }
  class(res) <- "SimCop"
  res
}
