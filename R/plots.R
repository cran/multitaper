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

plot.mtm <- function(x, jackknife=FALSE, Ftest=FALSE, ftbase=1.01,
                     siglines=NULL, ...) {
    nw <- x$mtm$nw
    k <- x$mtm$k
    sub <- paste("NW = ", nw, " K = ", k, sep="")

    if(Ftest) {
        .plotFtest(x)
    } else {
        log <- match.call(expand.dots = )$log
        dBPlot <- FALSE
        if(!is.null(log) && log== "dB" ) {
            dBPlot <- TRUE
        }
        
        if(jackknife && !is.null(x$mtm$jk)) {
            if(dBPlot) {
                upperCI <- 10*log10(x$mtm$jk$upperCI)
                lowerCI <- 10*log10(x$mtm$jk$lowerCI)
                minVal <- 10*log10(x$mtm$jk$minVal)
                maxVal <- 10*log10(x$mtm$jk$maxVal)
            } else {
                upperCI <- x$mtm$jk$upperCI
                lowerCI <- x$mtm$jk$lowerCI
                minVal <- x$mtm$jk$minVal
                maxVal <- x$mtm$jk$maxVal
        }
            yRange <- c(minVal, maxVal)
            plot.spec(x, sub=sub, ylim=yRange, ...)
            lines(x$freq, upperCI, lty=2, col=2)
            lines(x$freq, lowerCI, lty=2, col=3)
        } else {
            plot.spec(x, sub=sub, ...)
        }
    }
}

