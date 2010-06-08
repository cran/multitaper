##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2010 Karim Rahim 

##     This file is part of the multitaper package for R.

##     The multitaper package is free software: you can redistribute it and
##     /or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 2 of the License, or
##     any later version.

##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.

##     You should have received a copy of the GNU General Public License
##     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

##     If you wish to report bugs please contact the author. 
##     karim.rahim@gmail.com
##     112 Jeffery Hall, Queen's University, Kingston Ontario
##     Canada, K7L 3N6


spec.mtm <- function(timeSeries,
                     k,
                     nw,
                     nFFT="default", 
                     centreWithSlepians=TRUE,
                     dpssIN=NULL,
                     returnZeroFreq=TRUE,
                     Ftest=FALSE,
                     jackknife=FALSE,
                     jkCIProb=.95,
                     maxAdaptiveIterations=100,
                     plot=TRUE,
                     na.action=na.fail,
                     returnInternals=FALSE,
                     ...) {

    series <- deparse(substitute(timeSeries))
    stopifnot(nw >= 0.5, k >= 1, nw <= 500, k <= 1.5+2*nw) 
    timeSeries <- as.ts(timeSeries)
    na.action(timeSeries)
    deltaT <- deltat(timeSeries)
    if(k==1) {
        Ftest=FALSE
        jackknife=FALSE
    }
    timeSeries <- as.double(timeSeries)
    n <- length(timeSeries)
    
    stopifnot(n > 8)

    sigma2 <- var(timeSeries) * (n-1)/n
 
    dw <- NULL
    ev <- NULL
    receivedDW <- TRUE
    if(!.is.dpss(dpssIN)) {
        receivedDW <- FALSE
        dpssIN <- dpss(n, k, nw=nw,
                       returnEigenvalues=TRUE)
        dw <- dpssIN$v*sqrt(deltaT)
        ev <- dpssIN$eigen
    } else {
        dw <- .dpssV(dpss)
        ev <- .dpssEigen(dpss)
        if(ev == NULL) {
            ev <- dpssToEigenvalues(dw, nw)
        }
        dw <- dw*sqrt(deltaT)
    }

    if(nFFT == "default") {
        nFFT <- 2* 2^ceiling(log2(n))
    } else {
        stopifnot(is.numeric(nFFT))
    }
    stopifnot(nFFT >= n)
    
    nFreqs <- nFFT %/% 2 + as.numeric(returnZeroFreq)
    offSet <- if(returnZeroFreq) 0 else 1 
    scaleFreq <- 1 / as.double(nFFT * deltaT)
    
    swz <- NULL ## Percival and Walden H0
    ssqswz <- NULL
    
    if(centreWithSlepians) {
        swz <- apply(dw, 2, sum)
        ## zero swz where theoretically zero
        swz[1:(k/2)*2] <- 0.0
        ssqswz <- sum(swz**2)
        res <- .mweave(timeSeries, dw, swz,
                       n, k, ssqswz, deltaT)
        timeSeries <- timeSeries - res$cntr
    }
     
    taperedData <- dw*timeSeries
    
    nPadLen <- nFFT - n
    paddedTaperedData <- rbind(taperedData, matrix(0, nPadLen, k))
    cft <- mvfft(paddedTaperedData)
    cft <- cft[(1+offSet):(nFreqs+offSet),]
    sa <- abs(cft)**2
    
    resultFreqs <- ((0+offSet):(nFreqs+offSet-1))*scaleFreq 

    adaptive <-  NULL
    jk <- NULL
    PWdofs <- NULL
    if(!jackknife) {
        adaptive <- .mw2wta(sa, nFreqs, k,
                            sigma2, deltaT, ev)
    } else {
        stopifnot(jkCIProb < 1, jkCIProb > .5)
        adaptive <- .mw2jkw(sa, nFreqs, k,
                            sigma2, deltaT, ev)
        scl <- exp(qt(jkCIProb,adaptive$dofs)*
                   sqrt(adaptive$varjk))
        upperCI <- adaptive$s*scl
        lowerCI <- adaptive$s/scl
        minVal = min(lowerCI)
        maxVal = max(upperCI)
        jk <- list(varjk=adaptive$varjk,
                   bcjk=adaptive$bcjk,
                   sjk=adaptive$sjk,
                   upperCI=upperCI,
                   lowerCI=lowerCI,
                   maxVal=maxVal,
                   minVal=minVal)
    }

    ftestRes <- NULL

    if(Ftest) {
        if(is.null(swz)) {
            swz <- apply(dw, 2, sum)
        }
        ftestRes <- .HF4mp1(cft,
                            swz,
                            k,
                            ssqswz)
    }

    eigenCoef1 <- NULL
    wtCoef1 <- NULL
    
    if(returnInternals) {
        eigenCoef1 <- cft
        wtCoef1 <- sqrt(adaptive$wt)
    }
    auxiliary <- list(dpss=dpssIN,
                      ##eigenSpectra=sa,
                      eigenCoefs=eigenCoef1,
                      eigenCoefWt=wtCoef1,
                      nfreqs=nFreqs,
                      nFFT=nFFT,
                      jk=jk,
                      Ftest=ftestRes$F,
                      dofs=adaptive$dofs,
                      nw=nw,
                      k=k)
    
    spec.out <- list(origin.n=n,
                     method="Multitaper Spectral Estimate",
                     pad= nFFT - n,
                     spec=adaptive$s,
                     freq=resultFreqs,
                     series=series,
                     mtm=auxiliary)

    class(spec.out) <- c("mtm", "spec")
    
    if(Ftest) {
        class(spec.out) <- c("mtm", "Ftest", "spec")
    }
    
    if(plot) {
        plot.mtm(spec.out, jackknife=jackknife, ...)
        return(invisible(spec.out))
    } else {
        return(spec.out)
    }

}



centre <- function(x, nw=NULL, k=NULL, deltaT=NULL, trim=0) {
    na.fail(x)
    res <- NULL
    if(is.null(nw) && is.null(k) ) {
        res <- x - mean(x, trim=trim)
    } else {
        if(trim != 0) {
            warning(paste("Ignoring trim =", trim))
        }
        
        stopifnot(nw >= 0.5, k >= 1, nw <= 500, k <= 1.5+2*nw)
        if(is.null(deltaT)) {
            if(is.ts(x)) {
                deltaT <- deltat(ts)
            } else {
                warning("deltaT not sepecified using deltaT=1.")
                deltaT <- 1
            }
        }
        n <- length(x)
        dpssRes <- dpss(n, k=k, nw=nw,
                        returnEigenvalues=TRUE)
        dw <- dpssRes$v*sqrt(deltaT)
        ev <- dpssRes$eigen
        swz <- apply(dw, 2, sum)
        ## zero swz where theoretically zero
        swz[1:(k/2)*2] <- 0.0
        ssqswz <- sum(swz**2)
        res <- .mweave(x, dw, swz,
                       n, k, ssqswz, deltaT)
        res <- x - res$cntr
    }
    return(res)
}


