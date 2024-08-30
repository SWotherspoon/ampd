##' AMPD
##'
##' Implementation of the automatic multiscale-based peak detection
##' (AMPD) algorithm.  AMPD enables to detect peaks in noisy periodic
##' and quasi-periodic signals
##'
##' @param data list containing noisy data with quasi-periodic
##'   distribution of local maxima
##' @param L, maximum number of scales (corresponding to the maximum
##'   width of the window used)
##' @param extended default FALSE, if TRUE it uses an optimzed
##'   algorithm, which is capable to detect Maxima at the beginning
##'   and/or the end of the data. This was a bottleneck for the old
##'   algorithm. FALSE uses the traditional method.
##' @param splitting default FALSE, splits the data vector to achieve
##'   better computation performance. However, with this option, no
##'   other parameter than the positions of the maxima will be
##'   returned
##' @param splittingSize default NA, it is recommended to specify a
##'   splitting size, otherwise a splitting size of 000 is used. The
##'   splitting size is adapted anyway, so that at least 10 local
##'   maxima are included.
##' @return A list with different variables: LMS ($LMS), rescaled LMS
##'   ($rLMS), position of global minimum of the row-wise summation of
##'   LMS ($minLoc), position of maxima ($maximaLoc). Note that for
##'   big data vectors, only the positions of the local maxima will be
##'   returned. The others remain empty.
##' @keywords peak detection
##' @references Scholkmann et al. (2012). An Efficient Algorithm for
##' Automatic Peak Detection in Noisy Periodic and Quasi-Periodic Signals.
##' Algorithms, 5 (4), 588-603 http://www.mdpi.com/1999-4893/5/4/588
##' @importFrom grDevices heat.colors
##' @importFrom graphics abline image par plot points
##' @importFrom stats fitted lm runif sd
##' @export
##' @examples
##' t <- seq(0,2,0.005)
##' data <- sin(25*t)*sin(0.3*t)+0.4*t
##' dataNoise <- jitter(data,1000)
##' result <- AMPD(dataNoise)
##' result2 <- AMPD(dataNoise, extended=TRUE)
##' par(mfrow=c(1,2))
##' plot(dataNoise, main="traditional algorithm", type="l")
##' points(result$maximaLoc, dataNoise[result$maximaLoc],col="red")
##' plot(dataNoise, main="extended algorithm", type="l")
##' points(result2$maximaLoc, dataNoise[result2$maximaLoc],col="red")

AMPD <- function(data, L=NA, extended=FALSE, splitting=FALSE, splittingSize=NA) {

  ## (1) Linear detrending of the signal x
  N <- length(data)
  f <- fitted(lm(data~seq(1:N)))
  dataDetrended <- data-f

  ## if no splittingSize is given, the algorithm starts with 1000, but
  ## will automatically adapt it to a size which will involve at least
  ## 10 maxima
  if(splitting) {
    if(is.na(splittingSize))
      splittingSize <- min(1000,length(dataDetrended))
    p <- 0
    ScanWindow <- ceiling(splittingSize/2)-1
    start <- 1
    end <- start+splittingSize-1
    while(splitting) {
      if(end>=N) {
        end <- N
        start <- end - splittingSize + 1
        splitting <- FALSE
      }
      maximas <- doCalculation(dataDetrended[start:end],splittingSize,ScanWindow,extended=TRUE, LMS=FALSE)$maximaLoc + start - 1
      if(!is.element(p[length(p)],maximas)) p <- p[-length(p)]
      p <- append(p,maximas)
      if(start==1) {
        goodWindowLength <- ceiling(mean(diff(p))*10)
        splittingSize <- goodWindowLength
        ScanWindow <- ceiling(splittingSize/2)-1
      }

      start <- start + splittingSize-floor(splittingSize/10*1.5)
      end <- start+splittingSize-1
    }
    p <- unique(p)
    list("LMS" = NULL, "rLMS" = NULL, "minPos"=NULL, "maximaLoc"=p)
  } else {
    if(is.na(L)) L <- ceiling(N/2)-1
    doCalculation(dataDetrended,N,L,extended,LMS=TRUE)
  }
}

doCalculation <- function(dataDetrended,N,L,extended,LMS) {
  ## Calculates the Local Maxima Scalogram (LMS)
  M <- array(1,c(N,L))
  if(extended)
    M <- sapply(seq_len(L), function(k) localExtremaAdvanced(dataDetrended,N,k,min(dataDetrended)))
  else
    M <- sapply(seq_len(L), function(k) localExtrema(dataDetrended,N,k))

  if (extended) {
    ## weight of each row decreases with increasing row number this
    ## makes the algorithm more stable.
    colSumM <- colSums(M)
    f <- fitted(lm(colSumM~seq_along(colSumM)))
    colSumM <- c((colSumM-f))
    minLoc <- which.min(colSumM)
  } else {
    ## the old simple search for the global minimum
    minLoc <- which.min(colSums(M))
  }
  M <- t(M)
  ## calculate rescaled LMS
  rescaled <- M[1:minLoc,,drop=FALSE]

  ## search for local maximum, which means for rows with standard deviation = 0
  p <- which(apply(rescaled,2,sd)==0)

  if(extended && length(p)>0) {
    ## make sure that the local maximum is not at the edges of the dataset
    if(p[length(p)]==length(p)) p <- p[-length(p)]
    if(p[1]==1) p <- p[-1]
  }

  if(LMS)
    list("LMS" = M, "rLMS" = rescaled, "minPos"=minLoc, "maximaLoc"=p)
  else
    list("LMS" = NULL, "rLMS" = NULL, "minPos"=NULL, "maximaLoc"=p)
}

localExtremaAdvanced <- function(x, N, k, dataMin) {
  ## min(data) will be added to the data vector at the beginning and at the end
  alpha <- 1
  m <- runif(N)*alpha+1
  x <- c(rep(dataMin,k),x,rep(dataMin,k))
  ## checks if data point is max around +- k
  locMax <- ((x[(k+1):(N+k)]>x[1:N]) & (x[(k+1):(N+k)]>x[(k+1+k):(N+k+k)]))
  m[locMax] <- 0
  return(m)
}

localExtrema <- function(x, N, k) {
  alpha <- 1
  locMax <- logical(N)
  m <- runif(N) * alpha + 1
  locMax[(k+1):(N-k)] <- ((x[(k+1):(N-k)]>x[(1):(N-2*k)]) & (x[(k+1):(N-k)]>x[(k+1+k):(N)]))
  m[locMax] <- 0
  return(m)
}
