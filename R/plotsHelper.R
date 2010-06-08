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

.plotFtest <- function(x, ftbase=1.01, siglines=NULL, ...) {

    stopifnot(!is.null(x$mtm$Ftest) ||
              "Ftest" %in% class(mtm))
    
    log <- match.call(expand.dots = )$log
    xlab <- match.call(expand.dots = )$xlab
    ylab <- match.call(expand.dots = )$ylab

    if(is.null(xlab)) xlab <- "Frequency"
    if(is.null(ylab)) ylab <- "F Test"
    
    ylog = "n"
    if(is.null(log) || log == "yes") {
        ylog = "y"
    }

    ftestVals = x$mtm$Ftest
    ftestVals[ftestVals < ftbase] <- ftbase
    plot(x$freq, ftestVals, log=ylog, ylab=ylab, xlab=xlab,
         type="l",...)
}
