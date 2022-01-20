\name{decomp2d}
\alias{decomp2d}
\title{Decomposition of an Image}
\description{
This function decomposes an image into two components according to frequency.
}
\usage{
decomp2d(z, method="wavelet",
    dct.frequency=NULL,
    emd.sm=FALSE, emd.spar=NULL, emd.tol=0.1^2, emd.maxiter=20,
    ept.tau=NULL, ept.process=c("average", "average"),
    ept.tol=0.1^2, ept.maxiter=20,
    pca.freqcomp=NULL,
    ssa.L=NULL, ssa.freqcomp=NULL,
    wavelet.highlevel=NULL)
}
\arguments{
\item{z}{matrix of an image.}
\item{method}{decomposition method of \code{"dct"} for discrete cosine transform, \code{"emd"} for bidimensional empirical mode decomposition, \code{"ept"} for ensemble pactch transform, \code{"pca"} for two-dimensional principal component analysis, \code{"ssa"} for two-dimensional singular spectrum analysis, and \code{"wavelet"} for wavelet transform.}
\item{dct.frequency}{threshold of frequencies for \code{"dct"}.}
\item{emd.sm}{specifies whether envelop is constructed by interpolation or local polynomial smoothing for \code{"emd"}. Use \code{FALSE} for interpolation or \code{TRUE} for local polynomial smoothing.}
\item{emd.spar}{specifies user-supplied smoothing parameter of local polynomial smoothing for \code{"emd"}.}
\item{emd.tol}{tolerance for stopping rule of sifting for \code{"emd"}.}
\item{emd.maxiter}{the maximum number of sifting for \code{"emd"}.}
\item{ept.tau}{a size parameter for \code{"ept"}: \code{tau[1]} for horizontal size and \code{tau[2]} for vertical size of a two-dimensional patch. When   \code{length(tau) is 1}, the horizontal and vertical size are the same.}
\item{ept.process}{specifies transform types for \code{"ept"}: \code{process[1]} for patch process and \code{process[2]} for ensemble process.}
\item{ept.tol}{tolerance for stopping rule of sifting for \code{"ept"}.}
\item{ept.maxiter}{the maximum number of sifting for \code{"ept"}.}
\item{pca.freqcomp}{numeric vectors of frequency components for \code{"pca"}.}
\item{ssa.L}{numeric vector with length 2 of window length for \code{"ssa"}.}
\item{ssa.freqcomp}{numeric vectors of frequency components for \code{"ssa"}.}
\item{wavelet.highlevel}{specifies resolution level of high-frequency component for \code{"wavelet"}.}
}
\details{
This function decomposes an image into frequency component and residue of two-dimensional image.
}
\value{
\item{fc}{high-frequency component decomposed from an image \code{z}.}
\item{residue}{residue image decomposed from an image \code{z}.}
}
\seealso{
\code{\link{empperiod}}.
}
\examples{
#### example : composite of two components having different frequencies
nr <- nc <- 128; x <- seq(0, 1, length=nr); y <- seq(0, 1, length=nc)

coscomp1 <- outer(cos(20 * pi * x), cos(20 * pi * y))
coscomp2 <- outer(cos(5* pi * x), cos(5 * pi * y))
cosmeanf <- coscomp1 + coscomp2

op <- par(mfcol=c(3,1), mar=c(0,0.5,2,0.5))
image(cosmeanf, xlab="", ylab="", main="a composite image",
    col=gray(0:100/100), axes=FALSE)
image(coscomp1, xlab="", ylab="", main="high-frequency component",
    col=gray(0:100/100), axes=FALSE)
image(coscomp2, xlab="", ylab="", main="low-frequency component",
    col=gray(0:100/100), axes=FALSE)

#### Decomposition by Wavelet Transform
outcoswr3 <- decomp2d(cosmeanf, method="wavelet", wavelet.highlevel=3)

par(mfcol=c(2,1), mar=rep(0.1, 4), oma=c(0,1.35,1.35,0))
image(outcoswr3$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE)
mtext("high-frequency component", side = 2, line = 0.3, cex=0.85, font=2)
mtext("level 3", side=3, line=0.1, cex=0.85, font=2)
image(outcoswr3$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE)
mtext("low-frequency component", side=2, line=0.3, cex=0.85, font=2)
par(op)
}
\keyword{nonparametric}
