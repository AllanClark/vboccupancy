vb_model2_la_orig<-function(W, X, y, alpha_0, beta_0, Sigma_alpha_0, Sigma_beta_0, LargeSample=TRUE, epsilon=1e-5)
{
  #Date: 10 September 2014
  #Allan E Clark
  #Multi-visit site occupancy model
  #Estimation is undertaken using variational bayes and Laplace approximation
  #----------------------------------------------------------------------------------
  
  #----------------------------------------------------------------------------------
  #Arguments
  #----------------------------------------------------------------------------------
  #y <- n by J matrix of presence absence data
  #n <- the number of locations
  #J <- the number of visits to the sites
  #Assumed that each site is visited J times
  #X <- occurence design matrix - add column of 1's for an intercept if required
  #W <- n by J by q array where q indicates how many covariates are added
  #to model the detection probabilities
  #alpha_0 <- starting value of detection covariate coefficients
  #beta_0 <- starting value of occurence covariate coefficients 
  #Sigma_alpha_0, Sigma_beta_0 - the variance covariance matrix of alpha_0 and beta_0
  #LargeSample<-TRUE - indicates that the number of sites is 'large' and that an 
  #approximation to B(mu, sigma^2) is used instead of integrations
  #-----------------------------------------------------------------------------------
  
  #load the required functions
  bx<-function(x){log(1+exp(x))}
  
  b1<-function(x)
  {
    #the first deriv of bx
    1/(1+exp(-x))
  }
  
  b2<-function(x)
  {
    #the second deriv of bx
    exp(-x)/( (1+exp(-x))^2 )
  }
  
  Ebx <- function(x, mu, sigma2)
  {
    #bx(mu+sigma*x)*dnorm(x)
    argx<-mu+sqrt(sigma2)*x
    log(1+exp(argx))*dnorm(x)
  }
  
  B<-function(mu,sigma2)
  {
    #This does Brute force integrations
    #Note that the range is limited to -2000 to 2000
    integrate( Ebx, lower=-200, upper =200, mu=mu, sigma2=sigma2)$value
  }
  
  B2<-function(mu, sigma2)
  {
    #Approximation to B(mu,sigma2) for large sample sizes
    bx(mu)+b2(mu)*sigma2
  }
  
  B0.approx <- function(mu,sigma2)
  {
    #John Ormerod code
    sigma <- sqrt(sigma2)
    vB0 <- mu*pnorm(mu/sigma) + sigma2*dnorm(mu,sd<-sigma)
    vB0 <- vB0 + (0.6931472*exp(0.3298137*sigma2-0.8121745*mu))*pnorm( mu/sigma - 0.8121745*sigma)
    vB0 <- vB0 + (0.6931472*exp(0.3298137*sigma2+0.8121745*mu))*pnorm(-mu/sigma - 0.8121745*sigma)
    return(vB0)
  }
  
  B1.approx <- function(mu,sigma2)
  {
    #John Ormerod code
    sigma <- sqrt(sigma2)
    muOnSigma <- mu/sigma
    a1sigma <- 0.8121745*sigma
    a1mu <- 0.8121745*mu
    halfa1sqsigma2 <- 0.3298137*sigma2
    vB1 <- pnorm(muOnSigma)
    vB1 <- vB1 + 0.5629565*(exp(halfa1sqsigma2+a1mu)*pnorm(-muOnSigma - a1sigma) - exp(halfa1sqsigma2-a1mu)*pnorm( muOnSigma - a1sigma))
    return(vB1)
  }
  
  B2.approx <- function(mu,sigma2)
  {
    #John Ormerod code
    sigma <- sqrt(sigma2)
    muOnSigma <- mu/sigma
    a1sigma <- 0.8121745*sigma
    a1mu <- 0.8121745*mu
    halfa1sqsigma2 <- 0.3298137*sigma2
    vB2 <- -0.1259130*dnorm(mu,sd<-sigma)
    vB2 <- vB2 + 0.4572189*((exp(halfa1sqsigma2-a1mu))*pnorm(muOnSigma - a1sigma) + (exp(halfa1sqsigma2+a1mu))*pnorm(-muOnSigma - a1sigma))
    return(vB2)
  }
  
  XWX.fun2 <- function(X,w) 
  {
    #return( t(X)%*%diag(as.vector(w)))%*%X ) # Much slower
    #return( t(X*as.vector(w))%*%X  )         # Much faster
    return( crossprod(X*as.vector(w),X)  )    # faster than second version
  }
  
  diagSet <- function(d) 
  {
    #John Ormerod code
    return(  d*((1:d)-1) + (1:d) )
  }
  
  diagElements <- function(A) 
  {
    #John Ormerod code
    return( A[diagSet(nrow(A))] )
  }
  
  trace <- function(A) 
  {
    #John Ormerod code
    #return(sum(diag(A)))          # Much slower
    return(sum(diagElements(A)))  # Much faster 
  }
  
  X_row.Y_col=function(X,Y)
  {
    #does row by column matrix multiplication
    #ie row_i of matrix X %*% column_j of matrix Y
    #returns a column vector with the required elements
    #appears faster than a simple loop
    index=matrix(1:NROW(X))
    matrix(apply(index,1,function(x, m1=X,m2=Y){m1[x,]%*%m2[,x]}))
  }
  
  logp<-function(par, W, X, Y, P_tilde, p_tilde, p, alpha_0, SigmaInv_alpha_0, beta_0, SigmaInv_beta_0, Exp=0)
  {
    #log(q(alpha, beta)) when 'Exp<-0'
    
    ncol_W <- NCOL(W)
    
    alpha<- matrix(par[1:ncol_W], ncol=1)
    beta<- matrix(par[-c(1:ncol_W)] , ncol=1)
    
    alpha_x<-W%*%alpha
    alpha_diff <- alpha - alpha_0
    beta_x<-X%*%beta
    beta_diff <- beta - beta_0

    t1 <- crossprod(Y,P_tilde)%*%alpha_x - crossprod(p_tilde,bx(alpha_x) ) + crossprod(p,beta_x) - sum( bx(beta_x) )
    t2 <- 0.5*(1-Exp)*( crossprod(alpha_diff,SigmaInv_alpha_0)%*%alpha_diff + crossprod(beta_diff,SigmaInv_beta_0)%*%beta_diff  ) 
    
    Logp <- t1[1]  - t2[1] 
    return(Logp)
  }
  
  #----------------------------------------------------------------------------------
  
  #Declare certain matrices and constants
  #---------------------------------------
  #create a matrix where the elements of y are stored one vector below the other
  Y<-matrix(t(y), ncol=1) 
  X<-as.matrix(X) #X stored as a matrix
  n<-NROW(X)
  N<-length(Y)
  pn<-NCOL(X) #number of columns of X
  qn<-dim(W)[3] #number of columns of W
  
  const1 <- (pn+qn)*(log(2*pi)+1)
  const2<-determinant(Sigma_alpha_0, logarithm=T)$modulus[1] + determinant(Sigma_beta_0, logarithm=T)$modulus[1]
  const3<- 0.5*(const1+const2)
  
  SigmaInv_alpha_0 <- chol2inv(chol(Sigma_alpha_0))
  SigmaInv_beta_0 <- chol2inv(chol(Sigma_beta_0))
  SigmaInv_times_alpha_0 <- SigmaInv_alpha_0%*%alpha_0
  SigmaInv_times_beta_0 <- SigmaInv_beta_0%*%beta_0
  #---------------------------------------
  
  #a matrix containing how many times each of the sites are visited
  J<-dim(W)[2]
  V<-matrix(rep(J),n)
  
  #create the W matrix
  #Note that here were are assuming that V(i) are all equal
  W_vb <- cbind(1,c(t(W[,,qn]))) #this line is only correct if we have 1 W regressor!!!
  print(dim(W_vb))
  
  #the prior mean and covariance matrix of alpha and beta
  #------------------------------------------------------
  #the starting value for the alpha and beta vector (and covariance matrices)
  #--------------------------------------------------------------------------
  
  alpha <- alpha_0*rnorm(1)  #some randomization
  beta <- beta_0 *rnorm(1) #some randomization
  Sigma_alpha <- Sigma_alpha_0 
  Sigma_beta <- Sigma_beta_0 
  Sigma_alpha_inv <-solve(Sigma_alpha_0)
  Sigma_beta_inv <- solve(Sigma_beta_0)

  print(dim(alpha))
  
  oldparams<- c(alpha, Sigma_alpha_inv, beta, Sigma_beta_inv)
  nparams <- length(oldparams)
  
  #indentify which at which of the sites the species were observed!
  #----------------------------------------------------------------
  pres_abs <- apply(y,1,max)
  obs_vec <- which(pres_abs==1) #where they were observed
  not_obs_vec <- which(pres_abs==0) #where they were not observed
  n_obs <- length(obs_vec)
  n_not_obs <- length(not_obs_vec)
  
  #initial starting values of E(Z_tilde) and E(Z)
  #----------------------------------------------
  #if the species is observed at a site the element should be 1
  p_tilde <- matrix(1, ncol=1, nrow=N)
  p<-matrix(1, ncol=1, nrow=n)
  
  #the index of the location of the starting values of Z_tilde
  #Z_tilde <- [(z1,z1,...z1),....,(zn,.......,zn)]^T
  #only keep the ones where the species has not been observed
  starts_z <- matrix(c(1+cumsum(V)-V)[not_obs_vec] , ncol=1)
  Indices<-apply(starts_z,1,function(x,J.=J) x:(x+J.-1))
  c_Indices<-c(Indices)
  
  t1<-apply(X[not_obs_vec,],1, function(x,beta.=beta){crossprod(x,beta.)}) 
  t2<-apply( matrix( bx(W_vb[c_Indices,]%*%alpha), nrow=n_not_obs , byrow=T) , 1, sum)   
  p[not_obs_vec] <- 1/(1+ exp(-t1 + t2) ) 
  p_tilde[c_Indices]<-rep( p[not_obs_vec], each=J) 
  P_tilde<-diag(c(p_tilde))
  
  #Marginal likelihood approximation <- E(Log_p(y,z,alpha,beta))- E(q(alpha, beta, z))
  #-----------------------------------------------------------------------------------  
  Logp <- logp(c(alpha,beta), W=W_vb, X=X, Y=Y, P_tilde=P_tilde, p_tilde=p_tilde, p=p, alpha_0=alpha_0, SigmaInv_alpha_0=SigmaInv_alpha_0, beta_0=beta_0, SigmaInv_beta_0=SigmaInv_beta_0,Exp=1)
  #Logp <- Logp -.5*( const1 + log(det(Sigma_alpha_0)) + log(det(Sigma_beta_0)) )  #E(Log_p(y,z,alpha,beta))
  #Logp
  Logp<-Logp -0.5*( const1 + determinant(Sigma_alpha_0, logarithm=T)$modulus[1] + determinant(Sigma_beta_0, logarithm=T)$modulus[1] )  #E(Log_p(y,z,alpha,beta))
  #Logp
  
  #the E(q(alpha, beta, z)) part
  #------------------------------
  
  #E_log_pz <- sum(p[not_obs_vec]*log(p[not_obs_vec]/(1-p[not_obs_vec])) + log(1-p[not_obs_vec])) #gave numerical errors sometimes
  E_log_pz <-sum(log(p[not_obs_vec]^p[not_obs_vec]) + log( (1-p[not_obs_vec])^(1-p[not_obs_vec]) ))  
  E_log_p_alpha_beta <- -0.5*( const1 + determinant(Sigma_alpha, logarithm=T)$modulus[1] + determinant(Sigma_beta, logarithm=T)$modulus[1] )
  Log_ml <- Logp - E_log_p_alpha_beta[1] - E_log_pz
    
  #  cat("------------------------------------------------------------------------------------","\n")
  #  cat("STARTING VALUES","\n")
  #  cat("Log_ml <- ", Log_ml,", alpha <- ", round(alpha,digits<-3) , ", beta <- ", round(beta,digits<-3), "\n")
  #  cat("------------------------------------------------------------------------------------","\n")
  
  old_Log_ml<-Log_ml
  
  diff_Log_ml<-1
  its<-0
  
  #epsilon <- 1e-10
  
  while (diff_Log_ml > epsilon)
  {
    its<-1+its
    
    #Perform the Newton Rhaphson algortihm to calculate alpha and beta at iteration t
    #---------------------------------------------------------------------------------
    #crit <- .5 #used to assess convergence of the parameters - Rule 1
    crit<-0 #used for Rule 3
    
    #These are used in the Newton Rhaphson algortihm but don't change within the 'while' loop
    g_alpha_1<- crossprod(W_vb, P_tilde%*%Y)
    g_beta_1<- crossprod(X, p)
    
    while (crit < nparams) 
    {
      alpha_x<-W_vb%*%alpha
      g_alpha <- g_alpha_1 - crossprod(W_vb, p_tilde*b1(alpha_x)) - SigmaInv_alpha_0%*%alpha + SigmaInv_times_alpha_0
      Sigma_alpha_inv <-   XWX.fun2(W_vb, p_tilde*b2(alpha_x)) + SigmaInv_alpha_0 #Note Sigma_alpha is not stored!
      alpha <- alpha + solve(Sigma_alpha_inv,g_alpha)
      
      beta_x<-X%*%beta
      g_beta <- g_beta_1 -crossprod(X, b1(beta_x)) - SigmaInv_beta_0%*%beta + SigmaInv_times_beta_0
      Sigma_beta_inv <- XWX.fun2(X, b2(beta_x)) + SigmaInv_beta_0 #Note Sigma_beta is not stored!
      beta <- beta + solve(Sigma_beta_inv, g_beta)
  
    #This section has not been fixed as yet
    #------------------------------------------------------------------------------------------
      #Stopping method 1
      #-----------------
      #       newparams<-c(alpha, Sigma_alpha, beta, Sigma_beta)
      #       crit<-sum(abs(newparams-oldparams))
      #       oldparams<-newparams
      
      #Stopping method 2 - this method will take longer in general!
      #------------------------------------------------------------
      #Logp <- logp(c(alpha,beta), W<-W_vb, X<-X, Y<-Y, P_tilde<-P_tilde, p_tilde<-p_tilde, p<-p, alpha_0<-alpha_0, SigmaInv_alpha_0<-SigmaInv_alpha_0, beta_0<-beta_0, SigmaInv_beta_0<-SigmaInv_beta_0, Exp<-0)
      #E_log_p_alpha_beta <- -.5*( const1 + log(det(Sigma_alpha)) -log(det(Sigma_beta)) )    
      #Log_ml <- Logp - E_log_p_alpha_beta[1] - E_log_pz
      #crit <- abs(Log_ml-old_Log_ml)
      #old_Log_ml <- Log_ml
      
      #Stopping method 3
      #-----------------
      #newparams<-c(alpha, Sigma_alpha, beta, Sigma_beta)
      #crit<-sum(abs(newparams-oldparams) <= epsilon) #STRONG CONDITION!
      #oldparams<-newparams  
    #------------------------------------------------------------------------------------------
    
      #Stopping method 3
      #-----------------
      newparams<-c(alpha, Sigma_alpha_inv, beta, Sigma_beta_inv)
      crit<-sum(abs(newparams-oldparams) <= epsilon) #STRONG CONDITION!
      #print(crit)
      oldparams<-newparams  
    }
    #Now we calculate the covariance matrices
    Sigma_alpha<-solve(Sigma_alpha_inv)
    Sigma_beta<-solve(Sigma_beta_inv)
    
    #Brute force maximisation
    #starts_par<-c(alpha, beta)#*rnorm(1)
    #starts_par<-c(Alpha, Beta)#*rnorm(1)
    #fit1 <- optim(par<-starts_par, fn<-logq, W<-W_vb, X<-X, Y<-Y, P_tilde<-P_tilde, 
    #             p_tilde<-p_tilde, p<-p, alpha_0<-alpha_0, SigmaInv_alpha_0<-SigmaInv_alpha_0, beta_0<-beta_0, 
    #             SigmaInv_beta_0<-SigmaInv_beta_0, hessian<-T, control<-list(fnscale<--1), method<-"BFGS")
    #fit1$par
    #alpha <- matrix(fit1$par[1:2], ncol=1)
    #beta <- matrix(fit1$par[3:4], ncol=1)
    #Sigma_big <- solve(-fit1$hessian)
    #Sigma_alpha <- Sigma_big[1:2,1:2]
    #Sigma_beta <- Sigma_big[3:4, 3:4]
    
    #Calculate the occupancy probabilities at the different locations
    #-----------------------------------------------------------------
    #E_log_pz<-0

    t1<-apply(X[not_obs_vec,],1, function(x,beta.=beta){crossprod(x,beta.)}) #correct (not_obs_vec by 1 vector)
    #W_vb_sel<-W_vb[c_Indices,]
    ag1<-W_vb[c_Indices,]%*%alpha
    #ag2<-matrix(diag(W_vb_sel%*%Sigma_alpha%*%t(W_vb_sel)), ncol=1)
    ag2<-X_row.Y_col(W_vb[c_Indices,],Sigma_alpha%*%t(W_vb[c_Indices,]))
    
    if (LargeSample== TRUE)
    {
      t2<-apply( matrix(B2(ag1, ag2), nrow=n_not_obs , byrow=T), 1, sum)
    }else
    {
      t2<-apply( matrix(apply(cbind(ag1,ag2),1, function(x){B(x[1],x[2])}), nrow=n_not_obs , byrow=T), 1, sum) 
    }
    p[not_obs_vec] <- 1/(1+ exp(-t1 + t2) ) 
    p_tilde[c_Indices]<-rep( p[not_obs_vec], each=J) 
    P_tilde<-diag(c(p_tilde))
    
    
    #Marginal likelihood approximation <- E(Log_p(y,z,alpha,beta))- E(q(alpha, beta, z))
    #-----------------------------------------------------------------------------------
    Logp <- logp(c(alpha,beta), W=W_vb, X=X, Y=Y, P_tilde=P_tilde, p_tilde=p_tilde, p=p, alpha_0=alpha_0, SigmaInv_alpha_0=SigmaInv_alpha_0, beta_0=beta_0, SigmaInv_beta_0=SigmaInv_beta_0,Exp=1)
    Logp <- Logp -const3  #E(Log_p(y,z,alpha,beta))
    
    #the E(q(alpha, beta, z)) part
    #------------------------------
    #E_log_pz <- sum(p[not_obs_vec]*log(p[not_obs_vec]/(1-p[not_obs_vec])) + log(1-p[not_obs_vec]))
    E_log_pz <-sum(log(p[not_obs_vec]^p[not_obs_vec]) + log( (1-p[not_obs_vec])^(1-p[not_obs_vec]) ))
    #E_log_p_alpha_beta <- -0.5*( const1 + log(det(Sigma_alpha)) + log(det(Sigma_beta)) )
    E_log_p_alpha_beta <- -0.5*( const1 + determinant(Sigma_alpha, logarithm=T)$modulus[1] + determinant(Sigma_beta, logarithm=T)$modulus[1] )
    
    Log_ml <- Logp - E_log_p_alpha_beta[1] - E_log_pz
    
    diff_Log_ml <- abs(Log_ml-old_Log_ml)
    old_Log_ml <- Log_ml
    
    cat("Iteration: ", its,", Log_ml <- ", round(Log_ml,digits<-6), ", alpha <- ", round(alpha,digits<-3) , ", beta <- ", round(beta,digits<-3), "\n")
  }
  #out<- rbind(cbind(alpha ,sqrt(diag(Sigma_alpha))), cbind(beta,sqrt(diag(Sigma_beta))) )
  #list(alpha=alpha, beta=beta, Sigma_alpha=Sigma_alpha, Sigma_beta=Sigma_beta, occup_p=p, out=out, Log_ml=Log_ml)

  list(alpha=alpha, beta=beta, Sigma_alpha=Sigma_alpha, Sigma_beta=Sigma_beta, occup_p=p, Log_mla=Log_ml)
}
