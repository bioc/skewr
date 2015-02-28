\name{getMethyLumiSet}
\alias{getMethyLumiSet}

\title{
Read idat files and create a \code{MethyLumiSet} object
}
\description{
This a wrapper function for \code{methylumIDAT} that does not require a vector of barcodes to be provided.
}
\usage{
getMethyLumiSet(path = getwd(), barcodes = NULL,
                norm = c("none", "illumina", "SWAN", "dasen"),
                bg.corr = TRUE)
}

\arguments{
  \item{path}{
The path to the directory containing the idat files.
}
  \item{barcodes}{
A vector of barcodes specifying which idat's to read.
}
  \item{norm}{
Should normalization be done on the resulting \code{MethyLumiSet}. The default is "none"
}
  \item{bg.corr}{
logical; if \code{TRUE}, an Illumina style background subtraction will be performed only if \code{norm} is set to \code{'illumina'}. Otherwise, it is ignored. If background subtraction without any normalization is desired, the \code{preprocess} method must be used.
}
}
\details{
If only \code{path} is provided, all idat's found in the given folder will be pulled. If only \code{barcodes} is given, corresponding idat's will be pulled from the current working directory. Both \code{path} and \code{barcodes} may be passed for finer control. The default is to pull all idat's found in the current working directory.
}
\value{
A \code{MethyLumiSet} object
}
\references{
Davis S, Du P, Bilke S, Triche T, Jr. and Bootwalla M (2014). methylumi: Handle Illumina methylation data. R package version 2.12.0.
}
\author{
Ryan Putney \email{ryanputney@gmail.com}
}
\note{
One would probably not normally want to use the preprocess option at this stage. It is more likely that a \code{MethyLumiSet} of the raw data will be desired. Then the \code{preprocess} method may be used to normalize the raw data or use background subtraction only on the raw data. See the vignette for example workflow.
}

\seealso{
\code{\link{getBarcodes}}, \code{\link{preprocess}}
}
\examples{
\donttest{
if(require('minfiData')) {
  path <- system.file("extdata/5723646052", package="minfiData")
  barcodes <- getBarcodes(path = path)
  methylumiset.raw <- getMethyLumiSet(path = path, barcodes = barcodes[1])
}
}
}
\keyword{ IO }
