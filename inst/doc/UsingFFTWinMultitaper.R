### R code from vignette source 'UsingFFTWinMultitaper.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 60)
multitaper1 <- packageDescription("multitaper")


###################################################
### code chunk number 2: <setupFunction
###################################################
require("fftwtools")


useFFTWGE <- 2^17 ## about 130 thousand values with padding.
useMVFFTWGE <- 2^16 ## Based on two tests; one of 5 columns and one of 10 columns of data


fft <- function(z, inverse = FALSE) {
    if(length(z) >= useFFTWGE) {
        fftwtools::fftw(z, inverse=inverse)
    } else {
        stats::fft(z, inverse=inverse)
    }
}

## do the same for now
mvfft <- function(z, inverse=FALSE) {
    if(dim(z)[1] >= useFFTWGE) {
        fftwtools::mvfftw(z, inverse=inverse)
    } else {
        stats::mvfft(z, inverse=inverse)
    }
}



###################################################
### code chunk number 3: <restore
###################################################
rm(fft, mvfft)


