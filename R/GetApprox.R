#' Approximate a copula by a histogram density
#'
#' Approximates the \dQuote{density} of a copula by a piece-wise constant function.
#'
#' This function provides two methods for subdividing the \eqn{d}-dimensional unit cube into hyper-rectangles, with \eqn{d} being passed to the parameter \code{dim}.  As most of the functions in this package which create a new copula return a function that can be evaluated at points in arbitrary dimensions, it is necessary to specify for which dimension \eqn{d} one wishes to calculate the approximation to the copula's \dQuote{density}.
#'
#' The first method (Approximation I) determines \eqn{2^m} hyper-rectangles (where \eqn{m} is the parameter \code{depth}), each containing the same probability mass with respect to the copula.  The second method (Approximation II) dividies the unit cube into \eqn{m^d} hyper-squares.
#'
#' These approximations can be interpreted as piecewise constant approximations of the copula's probability density function if the copula is absolutely continuous.  For futher details see \sQuote{References}.
#'
#'
#' @return \code{GetApprox} returns an object of \code{\link{class}} \sQuote{CopApprox} according to its inputs. The returned object is a list containing a matrix that holds the information of the approximation, the argument \code{Cop}, which approximation was determined, and other auxiliary information.
#'
#' The only method for objects of class \sQuote{CopApprox} implemented so far are for the generic function \code{\link{plot}}, and then only for the case if \code{dim} was 2.
#'
#' @param Cop A function defining the copula.
#' @param dim The approximation should be calculated on the dim-dimensional unit cube, defaults to 2.
#' @param depth The number of hyperrectangles to be used to devide the unit cube, defaults to 10 for Approximation I and to 32 for Approximation II.
#' @param type Whether Approximation I or Approximation II should be used, defaults to one.
#' @param TOL A numerical tolerance used when calculating Approximation I.
#'
#' @seealso \code{\link{plot.CopApprox}}
#'
#' @references Tajvidi, N. and Turlach, B.A. (2017). A general approach to generate random variates for multivariate copulae, \emph{Australian & New Zealand Journal of Statistics}. Doi:10.1111/anzs.12209.
#'
#' @author  Berwin A. Turlach <berwin.turlach@gmail.com>
#'
#' @examples
#' Cop <- NewMEVGumbelCopula(3)
#' CopApprox <- GetApprox(Cop, dim=2)
#' plot(CopApprox)
#'
#' @export
GetApprox <- function(Cop, dim=2, depth=ifelse(type==1, 10, 32),
                      type=1,
                      TOL=1000*.Machine$double.eps){
  if(type==1){
    res <- GetApprox1(Cop, dim, depth, TOL)
  }else{
    res <- GetApprox2(Cop, dim, depth)
  }
  class(res) <- "CopApprox"

  res
}

GetApprox1 <- function(Cop, dim=2, depth=10,
                       TOL=1000*.Machine$double.eps){
  res <- matrix(NA, nrow=2^depth, ncol=2*dim)
  ires <- 1

  tracker <- matrix(0, nrow=depth, ncol= 2*dim + 4)
  il <- 1:dim
  ir <- il + dim
  ib <- c(il, ir)
  icand <- 2*dim + 1
  idim  <- 2*dim + 2
  icut  <- 2*dim + 3
  iprob <- 2*dim + 4
  tracker[1, ir] <- 1
  tracker[, iprob] <- 1/2^(1:depth)

  tracker[1, icand] <- 0.5
  tracker[1, idim] <- 1
  ind <- 1

  while(TRUE){

    while( ind < depth ){
      tracker[ind+1, ib] <- tracker[ind, ib]
      if(tracker[ind, icut] == 0){
        tracker[ind+1, il[tracker[ind, idim]]] <- tracker[ind, icand]
        tracker[ind, icut] <- 1
      }else{
        tracker[ind+1, ir[tracker[ind, idim]]] <- tracker[ind, icand]
        tracker[ind, icut] <- 2
      }
      ind <- ind + 1

      edgelen <- tracker[ind, ir] - tracker[ind, il]

      cutdim <- which.max(edgelen)
      CutFct <- .CreateCutFct(Cop, tracker[ind,il], tracker[ind, ir],
                             cutdim, offset=tracker[ind,iprob])

      cand <- stats::uniroot(CutFct,
                             lower = tracker[ind, il[cutdim]],
                             upper = tracker[ind, ir[cutdim]],
                             tol = TOL)$root

      tracker[ind, icand] <- cand
      tracker[ind, idim] <-  cutdim
    }

    res[ires, ib] <- res[ires+1, ib] <- tracker[ind, ib]
    res[ires,   il[tracker[ind, idim]]] <- tracker[ind, icand]
    res[ires+1, ir[tracker[ind, idim]]] <- tracker[ind, icand]
    ires <- ires + 2

    ind <- ind - 1
    while( ind > 0 && tracker[ind, icut] == 2 ){
      tracker[ind, icut] <- 0
      ind <- ind - 1
      if( ind == 0 ) break
    }
    if( ind == 0 )
      break
  }
  list(res=res, Copula=Cop, type=1)
}

GetApprox2 <- function(Cop, dim=2, depth=32){

  xgr <- seq(from=0, to=1, length=depth+1)[-(depth+1)]
  h <- 1/depth
  tmp <- vector("list", dim)
  for(i in 1:dim) tmp[[i]] <- xgr

  res <- data.matrix(expand.grid(tmp))
  dimnames(res) <- NULL

  prob <- function(x) .ProbRect(Cop, x, x+h)
  probs <- apply(res, 1, prob)

  list(res=res, Copula=Cop, type=2, probs=probs, h=h)
}

.CreateCutFct <- function(Cop, left, right, dim, offset=0){
  d <- length(right)

  tmp <- 0
  track <- rep(1, d)
  x <- right
  s <- 1
  ind <- 1
  repeat{
    tmp <- tmp + s * Cop(x)

    while( ind <= d ){
      if( ind != dim ){
        if( track[ind] == 1 ){
          track[ind] <- -1
          x[ind] <- left[ind]
          s <- -s
          ind <- 1
          break
        }
        track[ind] <- 1
        x[ind] <- right[ind]
        s <- -s
      }
      ind <- ind+1
    }
    if(ind > d)
      break
  }

  rm(track, x, s, ind)
  function(x){
    res <- tmp
    track <- rep(1,d)
    xx <- right
    xx[dim] <- x
    s <- -1
    ind <- 1
    repeat{
      res <- res + s * Cop(xx)

      while( ind <= d ){
        if(ind != dim){
          if( track[ind] == 1 ){
            track[ind] <- -1
            xx[ind] <- left[ind]
            s <- -s
            ind <- 1
            break
          }
          track[ind] <- 1
          xx[ind] <- right[ind]
          s <- -s
        }
        ind <- ind+1
      }
      if(ind > d)
        break
    }
    res-offset
  }
}

.CheckApprox <- function(obj){

  d <- NCOL(obj$res)/2
  il <- 1:d
  ir <- il+d
  prob <- function(x) .ProbRect(obj$Copula, x[il], x[ir])

  apply(obj$res, 1, prob)
}

.ProbRect <- function(Cop, left, right){
  d <- length(right)

  res <- 0
  track <- rep(1, d)
  x <- right
  s <- 1
  ind <- 1
  repeat{
    res <- res + s * Cop(x)

    while( ind <= d ){
      if( track[ind] == 1 ){
        track[ind] <- -1
        x[ind] <- left[ind]
        s <- -s
        ind <- 1
        break
      }
      track[ind] <- 1
      x[ind] <- right[ind]
      s <- -s
      ind <- ind+1
    }
    if(ind > d)
      break
  }
  if(res < 0){
    if(abs(res) > 1000*.Machine$double.eps){
      print(res)
      stop("Large negative probability calculated")
    }else{
      0
    }
  }else{
    res
  }
}

