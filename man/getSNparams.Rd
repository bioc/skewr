\name{getSNparams}
\alias{getSNparams}

\title{
Estimate parameters for finite mixture of Skew-Normal distributions
}
\description{
Utilizes \code{smsn.mix} from the \code{mixsmsn} package to find the parameters for a finite mixture of skew normal distributions to model the overall distribution of signal intensities for a subset of probes on the Illumina Infinium HumanMethylation450. The probes may be subset by type and methylated or unmethylated. It can also be specified whether the SNP(rs), imprinted(idmr), or ch probes should be included or filtered out prior to parameter estimation.
}
\usage{
getSNparams(MethyLumiSet, allele = c('M', 'U'),
            type = c('I-red', 'I-green', 'II'),
            snps = TRUE, idmr = TRUE, ch = FALSE)
}

\arguments{
  \item{MethyLumiSet}{
A \code{MethyLumiSet} object
}
  \item{allele}{
Should parameter estimation be done on the methylated or unmethylated signal intensities
}
  \item{type}{
Use the signal intensities for which probe type
}
  \item{snps}{
logical; should the rs probes be included in the dataset. The default is \code{TRUE}
}
  \item{idmr}{
logical; should the probes of imprinted gene loci be included in the dataset. The default is \code{TRUE}
}
  \item{ch}{
logical; should the ch probes be included in the dataset. The default is \code{FALSE}
}

}

\value{
A \code{Skew.normal} object as returned by \code{smsn.mix} from the \code{mixsmsn} package with the means and modes of the components added.
}
\references{
Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalwyk LC(2013). A data-driven approach to preprocessing Illumina 450k methylation array data. BMC Genomics, 14:293.

Prates MO, Cabral CRB, Lachos VH (2013).mixsmsn: Fitting Finite Mixture of Scale Mixture of Skew-Normal Distributions. Journal of Statistical Software, 54(12), 1-20. \url{http://www.jstatsoft.org/v54/i12/}
}
\author{

Ryan Putney \email{ryanputney@gmail.com}
}

\seealso{
\code{\link{subsetProbes}}
}
\examples{
\donttest{
if(require('wateRmelon')) {
  data(melon)
  mixes.raw.meth.II <- getSNparams(melon[,1], 'M', 'II')
}
}
}
\keyword{ models }

