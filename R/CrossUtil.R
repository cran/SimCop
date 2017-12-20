.crossmult <- function(h, u, w, v, iu, iw, iv){
  if(h == 1){
    v$val[iv] <- v$val[iv] + u$val[iu]*w$val[iw]
  }else{
    h <- h/2
    Recall(h, u, w, v, iu, iw+h, iv+h)
    Recall(h, u, w, v, iu+h, iw, iv+h)
    Recall(h, u, w, v, iu, iw, iv)
  }
  invisible()
}

.crossdeconv <- function(h, u, w, v, iu, iw, iv){
  v$val[iv] <- v$val[iv] - u$val[iu]*w$val[iw]
  i <- 1
  while(i < h){
    Recall(i, u, w, v, iu+i, iw, iv+i)
    Recall(i, u, w, v, iu, iw+i, iv+i)
    i <- i*2
  }
  invisible()
}

.crossdivide <- function(h, u, w, v, iu, iw, iv){
  w0 <- w$val[iw]
  w$val[iw] <- 0

  v$val <- u$val/w0
  w$val <- w$val/w0

  i <- 1
  while(i < h){
    .crossdeconv(i, w, v, v, iw+i, 1, iv+i)
    .crossdeconv(i, w, v, v, iw, iv+i, iv+i)
    i <- i*2
  }
  w$val <- w$val*w0
  w$val[iw] <- w0
  invisible()
}
