\name{preprocess}
\alias{preprocess}

\title{
Normalize a MethyLumiSet object using some popular choices
}
\description{
This is a wrapper function that allows normalizing of a MethyLumiSet using either a BeadStudio approximation, SWAN, or dasen. If desired, background correction only may be performed on the raw data.
}
\usage{
preprocess(MethyLumiSet, norm = c("none", "illumina", "SWAN", "dasen"),
           bg.corr = TRUE)
}

\arguments{
  \item{MethyLumiSet}{
A \code{MethyLumiSet} object
}
  \item{norm}{
The normalization method to be used
}
  \item{bg.corr}{
If TRUE, background subtarction using negative controls is performed. Ignored unless \code{norm} equals \code{'illumina'} or \code{'none'}
}
}
\details{
Both Illumina style normalization via controls and the background correct method are handled by \code{methylumi}. The \code{SWAN} and \code{dasen} normalization methods are both performed by \code{wateRmelon}
}
\value{
A \code{MethyLumiSet}
}
\references{
Davis S, Du P, Bilke S, Triche T, Jr. and Bootwalla M (2014). methylumi: Handle Illumina methylation data. R package version 2.12.0.

Maksimovic J, Gordon L, Oshlack A (2012). SWAN: Subset Quantile Within-Array Normalization for Illumina Infinium HumanMethylation450 BeadChips. Genome Biology, 13:R44.

Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalwyk LC(2013). A data-driven approach to preprocessing Illumina 450k methylation array data. BMC Genomics, 14:293.

Schalkwyk LC, Pidsley R, Wong CC, Touleimat N, Defrance M, Teschendorff A and Maksimovic J (2013). wateRmelon: Illumina 450 methylation array normalization and metrics. R package version 1.5.1.


}
\author{
Ryan Putney \email{ryanputney@gmail.com}
}

\seealso{
\code{\link{getMethyLumiSet}}
}
\examples{
\donttest{
if(require('wateRmelon')) {
  data(melon)
  melon.dasen <- preprocess(melon, norm = 'dasen')
}
}
}
\keyword{ manip }
