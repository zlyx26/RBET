#' @title MAC
#' @description Maximum adjusted chi-squared statistics for testing the difference between two objects.
#'
#' @param x A list of object.
#' @param y Another list of object.
#' @param k The sub-sampling size of MAC.
#'
#' @return The statistic value of MAC.
#' @export
#'
#' @import Rcpp

MAC2 <- function(x,y,k=20)
{
  x=as.matrix(x)
  y=as.matrix(y)
  dx=dim(x)[2]
  dy=dim(y)[2]

  if(dx!= 2 || dy!= 2)
    stop('The data must be in 2 columns')

  rx=dim(x)[1]
  ry=dim(y)[1]

  # if(rx != ry)
  #   stop('The sample size of X and Y must be the same')

  ## whether do sample
  if((rx+ry)<k){
    x <- x[order(x[,1]),]
    y <- y[order(y[,1]),]
    z <- rbind(x,y)
    res=mac2(x,y,z,rx,ry,k)
  }
  else{
    x <- x[order(x[,1]),]
    y <- y[order(y[,1]),]
    z <- rbind(x,y)
    ## sample k value to cal dist
    z <- z[order(z[,1]),]
    alpha = 0.01
    t <- floor(alpha*(rx+ry)+(1-2*alpha)*(rx+ry)/(k-1)*seq(0,(k-1),1))
    z <- z[t,]

    res=mac2(x,y,z,rx,ry,k)
  }
  return(res)
}

