\name{SimData}
\alias{SimData}
\docType{data}
\title{An example of simulated data that could be used with the function 'vb_model2_la'.}
\description{Simulated Data for the "vb_model2_la" function.}
\usage{data(SimData)}

\details{See SomeSimulations.}

\value{
  \item{y}{ - n by J matrix of detection/nondetection data where n=the number of locations and J = the number of visits to the sites. It is assumed that each site is visited a maximum of J times (although missing values are allowed). i.e. equal sampling effort is not assumed.}
  \item{X}{ - n by p occurrence design matrix.}
  \item{W.eg.l1}{ - A named list containing dataframes for the different site covariates (detection covariate). The dataframes are of dimension n by J where each row is associated with a site and each column represents a site visit. 
  
  e.g. W.eg.l1$W1[1, ] = the covariate values for site 1 for all of the visits. Note that some of the entries might be "NA" meaning that no visit took place at those occasions.}
  \item{W.eg.l2}{ - A named list containing dataframes for the different site covariates (detection covariate). The data for each site have been stacked one below the other either as a dataframe or as a matrix. e.g. W.eg.l2[[2]] is the data for site 2.
  
  CURRENTLY 'vb_model2_la' REQUIRES 'W' TO BE SPECIFIED AS 'W.eg.l1'.}
  \item{W_vb}{ - The site covariate (detection covariate) design matrix used internally in the function "vb_model2_la". Notice that dataframes contained in W.eg.l2 are row-bind to produce W_vb.}
}

\examples{
## Load the package as well as the data
##-------------------------------------
require(vboccupancy)
data(SimData)

SimData$y
SimData$X
SimData$W.eg.l1

}
\keyword{datasets}
