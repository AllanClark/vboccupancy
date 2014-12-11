# vb_Designs_check<-function(formula, Names)
# {
#   #perform some checks on the design matrices
#   #check that the names in the formula call are in the dataframes provided
#   
#   detVars<-all.vars(formula)
#   
#   if ( (sum(detVars[-1]%in% Names)==length(detVars[-1]))!= 1)
#   {
#     stop(print("\n \n CHECK YOUR FORMULA CALL. \n MISMATCH BETWEEN CALL AND DESIGN MATRICES. \n i.e. You included objects in the call: '~occupancy variables ~ detection variables' that does not appear in the design matrices."))
#     #stop()
#   } 
# }
# 
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
vb_ReqDesigns<-function(formula, design_mats)
{
  vb_Designs_check(formula, design_mats$Names)
  
  #create the W matrix
  W<-model.matrix(as.formula(paste("~",formula[3],sep="")), data=design_mats$W)
  #print(W)
  
  #create the X matrix
  f_xmat<-paste(formula[[2]])
  X<-model.matrix(as.formula(paste(f_xmat[1],f_xmat[3],sep="")), data=design_mats$X)
  #print(dim(X))
  
  list(W=W, X=X)
}

# vb_call<-function(formula, design_mats)
# {
#   req_design_mats<-vb_ReqDesigns(formula, design_mats)
#   
#   W<-req_design_mats$W
#   X<-req_design_mats$X
#   Y<-design_mats$Y$V1
# }
