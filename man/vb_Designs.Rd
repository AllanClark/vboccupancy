\name{vb_Designs}
\alias{vb_Designs}

\title{vb_Designs}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
vb_Designs(W, X, y)
}
\arguments{
  \item{W }{ - A named list of data.frames of covariates that vary within sites. i.e. The dataframes are of dimension n by J where each row is associated with a site and each column represents a site visit. 
  
  e.g. Suppose W contained three dataframes W1, W2 and W3; W$W1[1, ] = the covariate values for site 1 for all of the visits. Note that some of the entries might be 'NA' meaning that no visit took place at those occasions. }
  \item{X }{ - A named data.frame that varies at site level.}
  \item{y }{ - An n by J matrix of the detection, non-detection data, where n is the number of sites, J is the maximum number of sampling periods per site.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Load the package as well as the data
##-------------------------------------
require(vboccupancy)
data(SimData)

SimData$W.eg.l1

##At the moment we do not use W's of the type SimData$W.eg.l2
vb_Designs(W=SimData$W.eg.l1, X=SimData$X, y=SimData$y)

##You should get the following result
##-----------------------------------
# $Y
#    V1
# 1   0
# 2   0
# 3   1
# 4   0
# 5   0
# 6   0
# 7   1
# 8   0
# 9   0
# 10  0
# 
# $X
#        x1         x2
# 1 -0.5802233 -1.0949354
# 2  1.0468539  1.3166989
# 3 -1.3879419  0.7589495
# 4  0.7897816 -0.5628716
# 5  0.1315297 -0.4178415
# 
# $W
#         W1          W2         W3
# 1  -1.1475341  0.08380091  0.6508060
# 2   0.3516124 -1.61591051  0.6021102
# 3   1.4671572 -1.71228685  0.5024673
# 4   0.2061681 -0.03591530  0.6489031
# 5   0.8607333  0.06444337 -0.0692492
# 6  -1.1047705  1.39241296 -1.5268438
# 7  -1.1515037  0.13920676 -0.2837397
# 8  -1.0836617  0.29944258  1.4945716
# 9   0.6777174 -0.10518068 -0.3372867
# 10  0.5505103 -0.73756554  0.2584308
# 
# $Names
# [1] "x1" "x2" "W1" "W2" "W3"
# 
# $nvisits
# [1] 1 3 1 3 2
# 
# $pres_abs
# [1] 0 1 0 1 0


}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
