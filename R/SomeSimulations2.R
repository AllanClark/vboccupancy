#A simple example of how to construct y, X and W; the 
#absence-presence data, site covariates and observation covariates
#-----------------------------------------------------------------

require(MASS)
set.seed(1000)

#int, beta1, beta2
beta.param = c(-1.85, 1.5, -0.5)
n = 500

#create 2 site covariates used to model occupancy
x1 = runif(n, -2,2) 
x1 = (x1 - mean(x1)) / sd(x1)
x2 = runif(n, -5,5) 
x2 = (x2 - mean(x2)) / sd(x2) 
X = cbind(rep(1,n), x1, x2)
psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used 
z = rbinom(n, size=1, prob=psi)

J = 5 #the maximum number of surveys (some sites might have fewer visits)

#three observation covariates used to model the detection probs
#int, alpha1, alpha2, alpha3
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
# nvisits<-sample(1:J, n, replace=T)
# empty.sites<-which(nvisits!= J)
# 
# for (i in 1:length(empty.sites))
# {
#   #adds NA to sites with visits less than J
#   y[ empty.sites[i], (nvisits[empty.sites[i]]+1):J ] <- NA
#   
#   #adds NA to W entries with visits less than J
#   W[ empty.sites[i], (nvisits[empty.sites[i]]+1):J, ] <- NA 
# }

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
colnames(X.eg)<-c("X1","X2")
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
#W.eg.l1
#-----------------------------------------------------------------

#An alternate way of 'viewing' the site covariates is as follows:
#Create a list element; one for each site, where the data
#for each site have been stacked one below the other either as 
#a dataframe or as a matrix. e.g.
#W.eg.l2[[2]] is the data for site 2.

W.eg.l2=list(list())
#Here we assume that we do not have any missing values
nvisits<-rep(J,n) #so equal visits to each site

for (i in 1:n)
{
  if (nvisits[i]!=1)
  {
    dframe<-as.data.frame(W[i,1:nvisits[i],][,-1])
  }else
  {
    dframe<-as.data.frame(matrix(W[i,1:nvisits[i],][-1],nrow=1))
  }
  
  names(dframe)<-c("W1","W2","W3")
  W.eg.l2[[i]]<-dframe
}
#-----------------------------------------------------------------

#Two different ways of representing the observation covariates
#not run
#W.eg.l1
#W.eg.l2
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
#W.temp

nvisits<-apply(W.eg.l1[[1]],1,function(x){length(na.omit(x))})
#nvisits

W.eg<-NULL
nv<-length(nvisits)
for (i in 1:nv)
{
  #print(matrix( c(na.omit(W.temp[i,])), nrow=nvisits[i]))
  W.eg<-rbind(W.eg, matrix( c(na.omit(W.temp[i,])), nrow=nvisits[i]) )
}
#W.eg
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

SimData2<-list(y=Y.eg, X=X.eg, W.eg.l1=W.eg.l1, W.eg.l2=W.eg.l2, W_vb=W.eg)
#setwd("C:/Users/01369661/Dropbox/PhD/Year1/ISEC2014/OccupancyModel2/ISECflash/Results/Package/vboccupancy/data")
#save(SimData2,file="SimData2.rda")
#setwd("C:/Users/01369661/Dropbox/PhD/Year1/ISEC2014/OccupancyModel2/ISECflash/Results/Package/vboccupancy/")

#data(SimData2)
#head(SimData2$y)
#head(SimData2$X)
#length(SimData2$W.eg.l1)

#Assume that the formula used will be of the following form: formula1<- y~X1+X2~W1+W2+W3
#The occupancy model uses 2 covariates, X1 and X2; while
#the detection model uses 3 covariates W1, W2 and W3
#Intercepts are included in both models
#The function does not allow one to repress the intercept term

#Coefficients in the detection model
#alpha_0 <- matrix(0, ncol=1, nrow=4)
#Covariance matrix of the coefficients in the detection model
#Sigma_alpha_0 <- diag(4)*1000
#Coefficients in the occupancy process
#beta_0 <- matrix(0, ncol=1, nrow=3)
#Covariance matrix of the coefficients in the occupancy model
#Sigma_beta_0 <- diag(3)*1000

#design_mats<-vb_Designs(W=SimData2$W.eg.l1, X=SimData2$X, y=SimData2$y)
#vb_fit<-vb_model2_la(y~X1+X2~W1+W2+W3, design_mats=design_mats, alpha_0=alpha_0, beta_0=beta_0, Sigma_alpha_0=Sigma_alpha_0, Sigma_beta_0=Sigma_beta_0,LargeSample=TRUE, epsilon=1e-5)


#Here we use the large sample approximation and run the VB algorithm
#-------------------------------------------------------------------
#vb_1<-vb_model2_la(W=W, X=X, y=y, alpha_0=alpha_0, beta_0=beta_0, Sigma_alpha_0=Sigma_alpha_0, Sigma_beta_0=Sigma_beta_0,LargeSample=TRUE, epsilon=1e-5)

#dim(W)
#vb_model2_la<-function(W, X, y, 

# vb_Designs<-function(W, X, y)
# {
#   #create the required 'response' and 'regressor matrices'
#   #using all of the X and W data
#   #the output is stored as a named list 
#   
#   #create the Y matrix that will be used
#   Y<-matrix(na.omit(matrix(t(y), ncol=1)))
#   #col.names(Y)<-names(y)
#   pres_abs <- apply(y,1,max) #check if this will work for NA's
#   
#   #create the W matrix
#   W.temp<-NULL
#   nv<-length(W)
#   
#   for (i in 1:nv){W.temp<-cbind(W.temp, W[[i]])}
#   
#   nvisits<-apply(W[[1]],1,function(x){length(na.omit(x))})
#   n<-length(nvisits)
#   
#   W.out<-NULL
#   for (i in 1:n){W.out<-rbind(W.out, matrix( c(na.omit(W.temp[i,])), nrow=nvisits[i]) )}
#   colnames(W.out)<-names(W)
#   
#   list(Y=as.data.frame(Y), X=as.data.frame(X), W=as.data.frame(W.out), 
#        Names=c( colnames(X), colnames(W.out)), nvisits=nvisits,
#        pres_abs=pres_abs)
# }
# 
# design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=y)
# vb1<-vb_model2_la_2(y~X1+X2~W1+W2+W3, design_mats=design_mats, alpha_0=alpha_0, beta_0=beta_0, Sigma_alpha_0=Sigma_alpha_0, Sigma_beta_0=Sigma_beta_0,LargeSample=TRUE, epsilon=1e-5)

