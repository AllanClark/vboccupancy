vb_Designs<-function(W, X, y)
{
  #create the required 'response' and 'regressor matrices'
  #using all of the X and W data
  #the output is stored as a named list 
  
  #create the Y matrix that will be used
  Y<-matrix(na.omit(matrix(t(y), ncol=1)))
  #col.names(Y)<-names(y)
  pres_abs <- apply(y,1,max,na.rm=T) #check if this will work for NA's
  
  #create the W matrix
  W.temp<-NULL
  nv<-length(W)
  
  for (i in 1:nv){W.temp<-cbind(W.temp, W[[i]])}
  
  nvisits<-apply(W[[1]],1,function(x){length(na.omit(x))})
  n<-length(nvisits)
  
  W.out<-NULL
  for (i in 1:n){W.out<-rbind(W.out, matrix( c(na.omit(W.temp[i,])), nrow=nvisits[i]) )}
  colnames(W.out)<-names(W)
  
  list(Y=as.data.frame(Y), X=as.data.frame(X), W=as.data.frame(W.out), 
       Names=c( colnames(X), colnames(W.out)), nvisits=nvisits,
       pres_abs=pres_abs)
}
