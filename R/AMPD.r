#' AMPD
#'
#' Implementation of the automatic multiscale-based peak detection (AMPD) algorithm.
#' AMPD enables to detect peaks in noisy periodic and quasi-periodic signals
#'
#' @param data list containing noisy data
#' @param L, maximum number of scales (corresponding tothe maximum width of the window used)
#' @param plotting default FALSE, plots useful information if TRUE
#' @return A list with different variables: LMS ($LMS), rescaled LMS ($rLMS), position of global minimum of the row-wise summation of LMS ($minLoc), position of maxima ($maximaLoc)
#' @keywords peak detection
#' @references Scholkmann et al. (2012). An Efficient Algorithm for Automatic Peak Detection in Noisy Periodic and Quasi-Periodic Signals. Algorithms, 5 (4), 588-603 http://www.mdpi.com/1999-4893/5/4/588
#' @importFrom grDevices heat.colors
#' @importFrom graphics abline image par plot points
#' @importFrom stats fitted lm runif sd
#' @export
#' @examples
#' t = seq(0,2,0.005)
#' data = sin(25*t)*sin(0.3*t)+0.4*t
#' dataNoise = jitter(data,1000)
#' returnedList = AMPD(dataNoise,200,TRUE)
#' returnedList$maximaLoc

AMPD = function(data,L,plotting=FALSE){

  # Implementation of the automatic multiscale-based peak detection (AMPD) algorithm.
  # For more information see:
  #
  #   Scholkmann et al. (2012).
  #   An Efficient Algorithm for Automatic Peak Detection in Noisy Periodic
  #   and Quasi-Periodic Signals. Algorithms, 5 (4), 588-603
  #   http://www.mdpi.com/1999-4893/5/4/588
  #
  #__________________________________________________________________________
  # Version 1.0
  # Felix Scholkmann, felix.scholkmann@gmail.com
  # Implemented in R by Oliver Sieber, oliver.sieber@gmail.com
  #
  # For plotting, plotly should be installed: install.packages("plotly")
  #__________________________________________________________________________
  #
  #
  # INPUT
  # data:       1 D time series
  # L:          Maximum number of scales (corresponding to the maximum width of the window used)
  # plotting:   z = 1: show figures. z = 0: only perform calculation
  #
  # OUTPUT
  # M:        Local maxima scalogram (LMS)
  # rescaled: Rescaled local maxima scalogram
  # minLoc:   lambda parameter
  # p:        Indices of the peaks
  #__________________________________________________________________________


  # (1) Linear detrending of the signal x
  f = fitted(lm(data~seq(1:length(data))))
  dataDetrended = data-f

  # (2) Calculate the Local Maxima Scalogram (LMS)
  M = array(1,c(length(data),L))
  j = 1;
  for (i in 1:L){
    M[,j] = localExtrema(dataDetrended,i);
    j = j+1
  }

  M = t(M)

  minLoc = which.min(colSums(t(M)))

  p = which(apply(M[1:minLoc,],2,sd)==0)
  p = p - 1
  rescaled = M[1:minLoc,]


  if (plotting){
    par(mfrow=c(2,2))
    image(seq(1:length(data)),seq(1:L),t(M),col=heat.colors(256),main="Local maxima scalogram (LMS)",xlab="Index",ylab="Scale [#]")
    image(seq(1:length(data)),seq(1:minLoc),t(rescaled),col=heat.colors(256),main="Rescaled LMS",xlab="Index",ylab="Scale [#]")
    plot(colSums(t(M)),type="l",main="Row-wise summation of the LMS",xlab="Scale [#]",ylab=expression(sigma))
    abline(v=minLoc, col="red")
    plot(data,type="l",main="Detected peaks",xlab="Index",ylab="signal")
    points(p,data[p],col="red")
    if (length(p) <= 20){
      for (i in 1:length(p)){
        abline(v=p[i], col="gray")
      }
    }
  }

  newList <- list("LMS" = M, "rLMS" = rescaled, "minPos"=minLoc, "maximaLoc"=p)
  return(newList)
}



localExtrema = function(x,k){
  alpha = 1
  N = length(x)

  # (1) Calculate the local maxima for different skales given by k
  m = array(1,c(N-2*k,1))
  j = 1;
  for (i in (k+2):(N-k+1)){
    if (x[i-1] > x[i-k-1] && x[i-1] > x[i+k-1]){
      m[j,] = c(0)
    }
    else{
      m[j,] = 1+alpha* runif(1)
    }
    j = j+1
  }
  return(array(c(1 + alpha * runif(k+1), c(t(m)) ,1 + alpha * runif(k-1)),c(N,1)))
}

