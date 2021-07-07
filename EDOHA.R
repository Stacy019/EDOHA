EDOHA=function(S,ns,lambda1,lambda2,lambda3,lambda4,lambda5,convergence=1e-10,maxiter=2000,
                start="cold",startValue=NULL,verbose=FALSE,penalty="fused",penalizeD=TRUE){
  
  
  K=length(S)
  p=nrow(S[[1]])
  #Initialize the parameters
  # rho is the parameter for scaled Lagrangian form, we set rho=2.5 as is used in HGLASSO.
  rho=2.5
  #primal variable
  if(start=="cold"){
    Theta=list()
    Theta.tilt=list()
    Z=list()
    Z.tilt=list()
    V=list()
    V.tilt=list()
    V.tilt.tilt=list()
    for(k in 1:K){
      Theta[[k]]=diag(p)
      Theta.tilt[[k]]=matrix(0,p,p)
      Z[[k]]=diag(p)
      Z.tilt[[k]]=matrix(0,p,p)
      V[[k]]=diag(p)
      V.tilt[[k]]=matrix(0,p,p)
      V.tilt.tilt[[k]]=matrix(0,p,p)
    }

    
    #dual variable
    W.theta=list()
    W.Z=list()
    W.V=list()
    W.V.tilt=list()
    for(k in 1:K){
      W.theta[[k]]=matrix(0,p,p)
      W.Z[[k]]=matrix(0,p,p)
      W.V[[k]]=matrix(0,p,p)
      W.V.tilt[[k]]=matrix(0,p,p)
      }
  }else{
    Theta=startValue$Theta
    Z=startValue$Z
    V=startValue$V
    Theta.tilt=startValue$Theta.tilt
    Z.tilt=startValue$Z.tilt
    V.tilt=startValue$V.tilt
    V.tilt.tilt=startValue$V.tilt.tilt
    W.theta=startValue$W.theta
    W.Z=startValue$W.Z
    W.V=startValue$W.V
    W.V.tilt=startValue$W.V.tilt

  }
  
  criteria = 1e10 	
  i <- 1  	
  # While loop for the iterations
  while(criteria > convergence && i <= maxiter){
    if(verbose){
      print(paste("interation ",i))
    }
    oldTheta=Theta
    
    #update Theta
    for(k in 1:K){
      Theta[[k]]=update.theta(S[[k]], Theta.tilt[[k]], W.theta[[k]],ns[[k]],rho)
    }
    
    #update Z
    Z=update.Z(Z.tilt,W.Z,lambda1,lambda2,rho,penalizeD,penalty)
    
    #update V
    for(k in 1:K){
       V[[k]]=update.V(V.tilt[[k]],V.tilt.tilt[[k]],W.V[[k]],W.V.tilt[[k]],lambda3,lambda4,2*rho)
    }
    
   
    #update V.tilt.tilt
    V.tilt.tilt=update.Z(V,W.V.tilt,0,lambda5,rho,penalizeD,penalty)
    
    #update Theta.tilt, Z.tilt, V.tilt
    res=update.tilt(Theta,W.theta,Z,W.Z,V,W.V,rho)
    Theta.tilt=res[[1]]
    Z.tilt=res[[2]]
    V.tilt=res[[3]]
    
    #update W.theta
    for(k in 1:K){
      W.theta[[k]]=W.theta[[k]]+Theta[[k]]-Theta.tilt[[k]]
    }
    
    #update W.Z
    for(k in 1:K){
      W.Z[[k]]=W.Z[[k]]+Z[[k]]-Z.tilt[[k]]
    }
    
    #update W.V
    for(k in 1:K){
      W.V[[k]]=W.V[[k]]+V[[k]]-V.tilt[[k]]
    }
    
    
    #update W.V.tilt
    for(k in 1:K){
      W.V.tilt[[k]]=W.V.tilt[[k]]+V.tilt.tilt[[k]]-V[[k]]
    }
    
    error=0
    for(k in 1:K)
    error <- error + sum((Theta[[k]]-oldTheta[[k]])^2)/sum((oldTheta[[k]])^2)
    criteria=error
    i=i+1
  }
  for(k in 1:K){
    
    Z[[k]]=ifelse(abs(Z[[k]])<1e-3, 0, Z[[k]])
    V[[k]]=ifelse(abs(V[[k]])<1e-3, 0, V[[k]])
    Theta[[k]]=Z[[k]]+V[[k]]+t(V[[k]])
  }
  if(i>maxiter){
    warning(paste("the real exit converge criteria is ",criteria, " instead of ",convergence ))
  }
 
  a=list()
  b=list()
  vv=list()
  hubs=list()
  for(k in 1:K){
    vv[[k]]=V[[k]]
    diag(vv[[k]])=0
    a[[k]]=abs(vv[[k]])
    b[[k]]=apply(a[[k]],2,sum)
    hubs[[k]]=sort(which(b[[k]]>0))
  }

  return(list(Theta=Theta,Z=Z,V=V,hubs=hubs,ns=ns,Theta.tilt=Theta.tilt, Z.tilt=Z.tilt, 
              V.tilt=V.tilt,V.tilt.tilt=V.tilt.tilt, W.theta=W.theta, W.Z=W.Z, W.V=W.V,W.V.tilt=W.V.tilt,criteria=criteria))
}

update.theta=function(S,Theta.tilt,W.theta,n,rho){
  
  A=Theta.tilt-W.theta-n/rho*S
  a=eigen(A)
  d=a$values
  D=diag((d+sqrt(d^2+4*n/rho))/2)
  U=a$vectors
  Theta <- U%*%D%*%t(U)
  return (Theta)
}


 update.Z=function(Z.tilt,W.Z,lambda1,lambda2,rho,penalizeD,penalty){
   K=length(Z.tilt)
   p=nrow(Z.tilt[[1]])
   lam1=penalty.as.matrix(lambda1,p,penalize.diagonal=penalizeD)
   lam2=penalty.as.matrix(lambda2,p,penalize.diagonal=TRUE)
   A=list()
   for(k in 1:K)
   {A[[k]]=Z.tilt[[k]]-W.Z[[k]]}
   if(penalty=="fused")
   {
     # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
     if(K==2){Z = flsa2(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
     if(K>2){Z = flsa.general(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
   }
   if(penalty=="group")
   {
     # minimize rho/2 ||Z-A||_F^2 + P(Z):
     Z = dsgl(A,rho,lam1,lam2,penalize.diagonal=TRUE)
   }
   for(k in 1:K)
   {diag(Z[[k]])=diag(Z.tilt[[k]]-W.Z[[k]])}
   return(Z)
 }
 

 
   update.V=function(V.tilt,V.tilt.tilt,W.V,W.V.tilt,lambda3,lambda4,rho){
     p=nrow(V.tilt)
     V=matrix(0,p,p)
     C=(V.tilt+V.tilt.tilt)/2-(W.V-W.V.tilt)/2
     C1=C
     diag(C1)=0
     for(j in 1:p){
       Cj=C1[,j]
       Scj=sign(Cj)*pmax(abs(Cj)-lambda3/rho,0)
       a=max(1-lambda4/(rho*sqrt(sum(Scj^2))),0)
       V[,j]=a*Scj
     }
     diag(V)=diag(C)
     return(V)
  
 }


update.tilt=function(Theta, W.theta, Z, W.Z, V, W.V, rho){
  K=length(Theta)
  p=nrow(Theta[[1]])
  Gam=list()
  Theta.tilt=list()
  Z.tilt=list()
  V.tilt=list()
  for(k in 1:K){
    Gam[[k]]=rho/6*(Theta[[k]]+W.theta[[k]]-(Z[[k]]+W.Z[[k]])-(V[[k]]+W.V[[k]])-t(V[[k]]+W.V[[k]]))
    Theta.tilt[[k]]=(Theta[[k]]+W.theta[[k]])-Gam[[k]]/rho
    Z.tilt[[k]]=(Z[[k]]+W.Z[[k]])+Gam[[k]]/rho
    V.tilt[[k]]=(V[[k]]+W.V[[k]])+2*Gam[[k]]/rho
    }
  return(list(Theta.tilt,Z.tilt,V.tilt))
}

penalty.as.matrix <-function(lambda,p,penalize.diagonal)
  {
    # for matrix penalties:  check dim and symmetry:
    if(is.matrix(lambda))
    {
      if(sum(lambda!= t(lambda))>0) {stop("error: penalty matrix is not symmetric")}
      if(sum(abs(dim(lambda)-p))!=0 ) {stop("error: penalty matrix has wrong dimension")}
    }
    # for scalar penalties: convert to matrix form:
    if(length(lambda)==1) {lambda=matrix(lambda,p,p)}
    # apply the penalize.diagonal argument:
    if(!penalize.diagonal) {diag(lambda)=0}
    return(lambda)
}  
