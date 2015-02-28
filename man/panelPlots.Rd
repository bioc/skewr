\name{panelPlots}
\alias{panelPlots}

\title{
Plot the distributions of the probe intensities and the components of the skew-normal mixture model
}
\description{
Creates a panel of nine plots. Six of the plots represent the density of either the methylated intensity or the unmethylated intensity given by one of three subsets of the 485,577 total probes. These subsets include Type I-red, Type I-green, and Type II.The remaining three distributions give the density of the beta-values for these same three subsets. Each of the nine plots optionally displays the distributions of the "rs" SNP probes and the probes associated with imprinted genes(Pidsley,2013) as a series of 'tick' marks located above the x-axis.
}
\usage{
panelPlots(MethyLumiSet, typeIRedModels, typeIGreenModels, typeIIModels,
           plot = c("panel", "frames"), samp.num = NULL, frame.nums = 1:9,
           norm = "", idmr = TRUE, snps = TRUE)
}

\arguments{
  \item{MethyLumiSet}{
The \code{MethyLumiSet} object from which the mixture models were derived
}
  \item{typeIRedModels}{
A \code{list} of the Type I-red mixture models listed in the following order: methylated models followed by unmethylated models
}
  \item{typeIGreenModels}{
A \code{list} of the Type I-green mixture models listed in the following order: methylated models followed by unmethylated models
}
  \item{typeIIModels}{
A \code{list} of the Type II mixture models listed in the following order: methylated models followed by unmethylated models
}
  \item{plot}{
Should the output consist of panel plots--one panel per sample or a single panel if \code{samp.num} is specified; or should the function output separate plots corresponding to the frames, given by \code{frame.nums}, for a single sample. The default is \code{"panel"}. If set to \code{"frames"}, \code{samp.num} must be specified
}
  \item{samp.num}{
If plotting for a single sample is desired, for which sample. The number given simply refers to the \code{MethyLumiSet} column that corresponds to the sample of interest
}
  \item{frame.nums}{
if \code{plot} is set to \code{"frame"}, then frame.nums is a vector that specifies which frames of the panel to plot. The default is to plot all nine frames. The frames are numbered from 1 to 9 in column-major order starting with the top left. For example, to plot the four corners, use \code{frame.nums=c(1,3,7,9)}
}
  \item{norm}{
A character string which will be displayed as part of the main title for each plot. Useful in indicated which normalization method was used for the modeled and plotted data
}
  \item{idmr}{
logical; should the intensities of the idmr probes be plotted as a series of tick-marks above the x-axis. The default is \code{TRUE}
}
  \item{snps}{
logical; should the intensities of the rs probes be plotted as a series of tick-marks above the x-axis. The default is \code{TRUE}
}
}

\value{
No return value. Only plots are generated.
}
\references{
Prates MO, Cabral CRB, Lachos VH (2013).mixsmsn: Fitting Finite Mixture of Scale Mixture of Skew-Normal Distributions. Journal of Statistical Software, 54(12), 1-20. \url{http://www.jstatsoft.org/v54/i12/}
}
\author{
Ryan Putney \email{ryanputney@gmail.com}
}
\note{
Please refer to the vignette for an example workflow.
}



\seealso{
\code{\link{getSNparams}}
}
\examples{
\donttest{
if(require('minfiData')) {
  path <- system.file("extdata/5723646052", package="minfiData")
  methylumiset.raw <- getMethyLumiSet(path = path)
  mixes.raw.meth.I.red <- getSNparams(methylumiset.raw, 'M', 'I-red')
  mixes.raw.meth.I.green <- getSNparams(methylumiset.raw, 'M', 'I-green')
  mixes.raw.meth.II <- getSNparams(methylumiset.raw, 'M', 'II')
  mixes.raw.unmeth.I.red <- getSNparams(methylumiset.raw, 'U', 'I-red')
  mixes.raw.unmeth.I.green <- getSNparams(methylumiset.raw, 'U', 'I-green')
  mixes.raw.unmeth.II <- getSNparams(methylumiset.raw, 'U', 'II')
  mixes.I.red <- list(mixes.raw.meth.I.red, mixes.raw.unmeth.I.red)
  mixes.I.green <- list(mixes.raw.meth.I.green, mixes.raw.unmeth.I.green)
  mixes.II <- list(mixes.raw.meth.II, mixes.raw.unmeth.II)
  panelPlots(methylumiset.raw, mixes.I.red, mixes.I.green, mixes.II)
}
}
}

\keyword{ aplot }
