fitGev <- function(samp, shape = NA, scale = NA, location = NA,
                   lower = c( - Inf, 0, - Inf))

{
  ext.ll.nlminb <-
    function(x, y)
  {
    gamma <- x[1]
    sigma <- x[2]
    mu <- x[3]
    n.y <- length(y)
    if(any(((1 + (gamma * (y - mu))/sigma) < 4 * .Machine$
      double.eps)))
      return(10e9)
    if(abs(gamma) < 4 * .Machine$double.eps)
      return(n.y * log(sigma) + sum((y - mu)/sigma) + sum(exp(
                                                          - (y - mu)/sigma)))
    n.y * log(sigma) + (1/gamma + 1) * sum(log(1 + (gamma * (y - mu
    ))/sigma)) + sum((1 + (gamma * (y - mu))/sigma)^(-1/
                                                     gamma))
  }
  
  n <- length(samp)
  samp.pwm <- sort(samp)
  koeff1 <- ((1:n) - 1)/(n - 1)
  koeff2 <- (koeff1 * ((1:n) - 2))/(n - 2)
  b2 <- sum(koeff2 * samp.pwm)/n
  b1 <- sum(koeff1 * samp.pwm)/n
  b0 <- sum(samp.pwm)/n
  z <- (2 * b1 - b0)/(3 * b2 - b0) - log(2)/log(3)
  shape.pwm <- 7.859 * z + 2.9554 * z^2
  scale.pwm <- ((2 * b1 - b0) * shape.pwm)/(gamma(1 +
			shape.pwm) * (1 - 2^( - shape.pwm)))
  location.pwm <- b0 + (scale.pwm * (gamma(1 + shape.pwm) -
			1))/shape.pwm

  res <- stats::nlminb(start = c(ifelse(missing(shape),  -
                                                         shape.pwm, shape),ifelse(missing(scale),
                                                                                  scale.pwm, scale), ifelse(missing(
                                                                                                         location), location.pwm, location)), objective = ext.ll.nlminb, lower = lower,
                       y = samp)
  ests <- c(location=res$par[3],scale=res$par[2],shape=res$par[1])
  list(mle=ests,pwm = c(
                    location=location.pwm, scale=scale.pwm, shape=-shape.pwm))
  
}

transformToFrechet            <- function(obs.x,verbose=F)
{
  est.res<- fitGev(obs.x)$mle
  gamma<- est.res["shape"]
  sigma<- est.res["scale"]
  mu<- est.res["location"]
  if (verbose) {
    cat("gamma =",round(gamma,3),"\n")
    cat("sigma =",round(sigma,3),"\n")
    cat("mu =",round(mu,3),"\n")
  }
  if (gamma >0){
    left.point<-mu-sigma/gamma
    valid.check<- obs.x >= left.point
    if (!all(valid.check))
    {
      rem.index<-seq(along=obs.x)[!valid.check]
      cat("Observations", obs.x[!valid.check],"with index", rem.index, "is outside the support of distribution\n\n", "Estimated left point of distribution is",left.point,"\n" )
    }
  }
  if (gamma <0){
    right.point<- mu-sigma/gamma
    valid.check<- obs.x <= right.point
    if (!all(valid.check))
    {
      rem.index<-seq(along=obs.x)[!valid.check]
      cat("Observations", obs.x[!valid.check],"with index", rem.index, "is outside the support of distribution\n\n", "Estimated right point of distribution is",right.point,"\n" )
    }
  }
  (1+(gamma*(obs.x-mu)/sigma))^(1/gamma)
}
