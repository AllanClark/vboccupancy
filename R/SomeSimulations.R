#A simple example of how to construct y, X and W; the 
#absence-presence data, site covariates and observation covariates
#-----------------------------------------------------------------

require(MASS)
set.seed(1000)

beta.param = c(-1.85, 1.5, -0.5)
n = 5

#create 2 site covariates used to model occupancy
x1 = runif(n, -2,2) 
x1 = (x1 - mean(x1)) / sd(x1)
x2 = runif(n, -5,5) 
x2 = (x2 - mean(x2)) / sd(x2) 
X = cbind(rep(1,n), x1, x2)
psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used 
z = rbinom(n, size=1, prob=psi)

J = 3 #the maximum number of surveys (some sites might have fewer visits)

#three observation covariates used to model the detection probs
alpha.param = c(-1.35, 1.0, 0.5, -.25) 
w1 = runif(n*J, -5,5)
w1 = (w1 - mean(w1)) / sd(w1)
w2 = runif(n*J, -1,1)
w2 = (w2 - mean(w2)) / sd(w2)
w3 = runif(n*J, 0,5)
w3 = (w3 - mean(w3)) / sd(w3)
W = array(dim=c(n,J,4))
W[,,1] = 1
W[,,2] = w1
W[,,3] = w2
W[,,4] = w3

p = matrix(nrow=n, ncol=J)
y = matrix(nrow=n, ncol=J)
for (j in 1:J)
{
  p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
  y[, j] = rbinom(n, size=1, prob=z*p[, j]) 
}
#-----------------------------------------------------------------

#Now lets simulate the number of visits to each of the sites
#i.e. we need to set some of the y and W entries equal to 
#NA
nvisits<-sample(1:J, n, replace=T)
empty.sites<-which(nvisits!= J)

for (i in 1:length(empty.sites))
{
  #adds NA to sites with visits less than J
  y[ empty.sites[i], (nvisits[empty.sites[i]]+1):J ] <- NA
  
  #adds NA to W entries with visits less than J
  W[ empty.sites[i], (nvisits[empty.sites[i]]+1):J, ] <- NA 
}

#Note W[i,,] are the covariate values for site i
#each row is for a specific visit
#-----------------------------------------------------------------

#An nxJ matrix of the observed measured data,
#where n is the number of sites and J is the 
#maximum number of observations per site.
Y.eg<-y
#-----------------------------------------------------------------

#siteCovs  
#A data.frame of covariates that vary at the site level. 
#This should have n rows and one column per covariate
X.eg=as.data.frame(cbind(x1,x2)) 
#-----------------------------------------------------------------

#obsCovs  
#the obsCovs matrix is constructed as per the 'unmarked' package
#i.e. W.eg.l1 is a named list of data.frames of covariates that 
#vary within sites.
#i.e. The dataframes are of dimension n by J
#where each row is associated with a site
#and each column represents a site visit. 
#e.g. W.eg.l1$W1[1, ] = the covariate values for site 1 for all of the
#visits. Note that some of the entries might be 'NA'
#meaning that no visit took place at those occasions.

W1=matrix(NA,nrow=n, ncol=J)
W2=matrix(NA, nrow=n, ncol=J)
W3=matrix(NA, nrow=n, ncol=J)
for (i in 1:n)
{ 
  W1[i,]<- W[i,,2]
  W2[i,]<- W[i,,3]
  W3[i,]<- W[i,,4] 
}

#colnames(W1)<-paste("W1.",1:J,sep="")
#colnames(W2)<-paste("W2.",1:J,sep="")
#colnames(W3)<-paste("W3.",1:J,sep="")

W.eg.l1<-list(W1=W1, W2=W2, W3=W3)
W.eg.l1
#-----------------------------------------------------------------

#An alternate way of 'viewing' the site covariates is as follows:
#Create a list element; one for each site, where the data
#for each site have been stacked one below the other either as 
#a dataframe or as a matrix. e.g.
#W.eg.l2[[2]] is the data for site 2.

W.eg.l2=list(list())

for (i in 1:n)
{
  if (nvisits[i]!=1)
  {
    dframe<-as.data.frame(W[i,1:nvisits[i],][,-1])
  }else
  {
    dframe<-as.data.frame(matrix(W[i,1:nvisits[i],][-1],nrow=1))
  }
  
  names(dframe)<-c("w1","w2","w3")
  W.eg.l2[[i]]<-dframe
}
#-----------------------------------------------------------------

#Two different ways of representing the observation covariates
W.eg.l1
W.eg.l2
#-----------------------------------------------------------------

#If the site covariates are provided as per W.eg.l1
#then we can construct W.eg as follows 
#(here W.eg is the way in which 'vb_model2_la')
#creates the site covariate matrix W)
#We assume that all sites are visited at least once
#although all might not be visited J times
#We further assume that there are no missing covariate
#values for those occasions sites are visited

W.temp<-NULL
nv<-length(W.eg.l1)
for (i in 1:nv)
{ 
  W.temp<-cbind(W.temp, W.eg.l1[[i]])
}
W.temp

nvisits<-apply(W.eg.l1[[1]],1,function(x){length(na.omit(x))})
nvisits

W.eg<-NULL
n<-length(nvisits)
for (i in 1:n)
{
  #print(matrix( c(na.omit(W.temp[i,])), nrow=nvisits[i]))
  W.eg<-rbind(W.eg, matrix( c(na.omit(W.temp[i,])), nrow=nvisits[i]) )
}
W.eg
colnames(W.eg)<-c("W1","W2","W3")

#-----------------------------------------------------------------

#If the site covariates are provided as per W.eg.l2
#then we can construct W.eg as follows
W.eg<-NULL
n <-length(W.eg.l2)
for (i in 1:n)
{
  W.eg<- rbind(W.eg, W.eg.l2[[i]])
}
W.eg
#-----------------------------------------------------------------


SimData<-list(y=Y.eg, X=X.eg, W.eg.l1=W.eg.l1, W.eg.l2=W.eg.l2, W_vb=W.eg)
#setwd("C:/Users/01369661/Dropbox/PhD/Year1/ISEC2014/OccupancyModel2/ISECflash/Results/Package/vboccupancy/data")
#save(SimData,file="SimData.rda")
#setwd("C:/Users/01369661/Dropbox/PhD/Year1/ISEC2014/OccupancyModel2/ISECflash/Results/Package/vboccupancy/")

