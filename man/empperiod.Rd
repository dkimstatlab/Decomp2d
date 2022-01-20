\name{empperiod}
\alias{empperiod}
\title{Calculating Empirical Period of an Image}
\description{
This function calculates empirical period of an image.
}
\usage{
empperiod(z)
}
\arguments{
\item{z}{matrix of an input image}
}
\details{
This function calculates empirical period of an image.
}
\value{
\item{rowperiod}{vector of empirical period between row-wise local maxima.}
\item{colperiod}{vector of empirical period between column-wise local maxima.}
}
\seealso{
\code{\link{decomp2d}}.
}
\examples{
nr <- 128; x <- y <- seq(0, 1, length=nr)

coscomp1 <- outer(cos(20 * pi * x), cos(20 * pi * y))
coscomp2 <- outer(cos(5* pi * x), cos(5 * pi * y))
cosmeanf <- coscomp1 + coscomp2

op <- par(mfrow=c(1,2), mar=c(2,2,1,1))
hist(empperiod(cosmeanf)$rowperiod, xaxt = "n", breaks=seq(4, 55, by=3), freq=FALSE,
    main="empirical period of vertical direction", xlab="")
axis(1, seq(4, 55, by=3), seq(4, 55, by=3))
hist(empperiod(cosmeanf)$colperiod, xaxt = "n", breaks=seq(4, 55, by=3), freq=FALSE,
    main="empirical period of horizontal direction", xlab="")
axis(1, seq(4, 55, by=3), seq(4, 55, by=3))
par(op)
}
\keyword{nonparametric}
