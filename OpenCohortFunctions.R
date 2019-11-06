#Functions for the Open cohort shiny app
library(swCRTdesign)

#Functions to generate design matrices
SWdesmat <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1))
  for(i in 1:(T-1)) {
    Xsw[i,(i+1):T] <- 1
  }
  return(Xsw)
}
CRXOdesmat<- function(T) {
  if((T-1)%%2 == 0) {
    
    Xcrxo <- matrix(data=0, ncol = T, nrow = T-1)
    Xcrxo[1:(T-1)/2, seq(1,T,2)] <- 1
    Xcrxo[((T-1)/2 + 1):(T-1), seq(2,T,2)] <- 1
    return(Xcrxo)    
  }
  
  
  if((T-1)%%2 == 1) {
    
    Xcrxo <- matrix(data=0, ncol = T, nrow = T)
    Xcrxo[1:(T)/2, seq(1,T,2)] <- 1
    Xcrxo[((T)/2 + 1):(T), seq(2,T,2)] <- 1
    return(Xcrxo)    
  }
  
}
plleldesmat <- function(T) {
  if((T-1)%%2 == 0) {
    Xpllel <- matrix(data=0, ncol = T, nrow =T -1)
    Xpllel[1:(T-1)/2,] <- 1
    return(Xpllel)
  }
  
  
  if((T-1)%%2 == 1) {
    Xpllel <- matrix(data=0, ncol = T, nrow =T)
    Xpllel[1:(T)/2,] <- 1
    return(Xpllel)
  }
  
}
pllelBLdesmat <- function(T) {
  #For now assume 50% of periods are baseline measurements.
  
  if((T-1)%%2 == 0) {
    Xpllel <- matrix(data=0, ncol = T, nrow =T -1)
    Xpllel[1:(T-1)/2,((T+1)/2):T] <- 1
    return(Xpllel)
  }
  
  
  if((T-1)%%2 == 1) {
    Xpllel <- matrix(data=0, ncol = T, nrow =T)
    Xpllel[1:(T)/2,(T/2 + 1):T] <- 1
    return(Xpllel)
  }
  
}

#########################################################################
#Functions to calculate variance of the treatment effect estimator
#This needs to be validated against other code for variance of 
#treatment effect estimator. (Karla's app?)
TreatEffVar_open <- function(expretent, totvar, m, desmat, rho, pi, tau){

#(decayrate, T, m, sig2E, sig2C, sig2eta, Xmat) {
 
  #Transforming the icc, cac, iac, and sd to variance components
  #icc=rho; cac= pi; iac=tau
  icc <- rho
  cac <- pi
  iac <- tau
  varclus <- (totvar)*(icc*cac)
  varCP <- (totvar)*(icc*(1-cac))
  varsubj <- (totvar)*(iac*(1-icc))
  varerror <- (totvar)*((1-iac)*(1-icc))
  
  Vi <- matrix(data=varclus + varsubj*(expretent)/m, nrow=ncol(desmat), ncol=ncol(desmat)) +
    diag(varCP + varerror/m + varsubj*(1-expretent)/m, nrow=ncol(desmat), ncol=ncol(desmat)) 
  
  K <- nrow(desmat) 
  T <- ncol(desmat)
  Xvec <-  as.vector(t(desmat))
  
  vartheta <-  1/(t(Xvec)%*%(diag(1,K)%x%solve(Vi))%*%Xvec -colSums(desmat)%*%solve(Vi)%*%(matrix(colSums(desmat),nrow=T, ncol=1))/K )
  return(vartheta)
}

SampleSize_open <- function(expretent, effsize, totvar, power, alpha, m, desmat, rho, pi, tau){
  #expretent is actually the expected proportion lost
  #The below line gives the same results as in Hooper et al, Stats in Med 2016, 
  #where the ceiling of the individually randomised sample size per arm is taken 
  #before calculating the number of clusters per sequence required.
  #n_ind <- 2*ceiling((2*totvar*(qnorm(1-alpha/2) + qnorm(power))^2)/(effsize^2))
  n_ind <- 2*((2*totvar*(qnorm(1-alpha/2) + qnorm(power))^2)/(effsize^2))
  deff_cluster <- 1 + (m-1)*rho
  r_OC <- (m*rho*pi + (1-rho)*tau*expretent)/(deff_cluster)
  
  K<- nrow(desmat)
  T <- ncol(desmat)
  
  X_bb <- sum(desmat)
  X_kb <- rowSums(desmat)
  X_bt <- colSums(desmat)
  
  deff_longitudinal <- (K^2)*(1-r_OC)*(1+ (T-1)*r_OC)/(4*(K*X_bb - sum(X_bt^2) + r_OC*(X_bb^2 + K*(T-1)*X_bb -(T-1)*sum(X_bt^2) - K*sum(X_kb^2))    ))
  
  return(deff_longitudinal*deff_cluster*n_ind/m)
  
  
  
  
}

Power_open <- function(expretent, effsize, totvar, alpha, m, desmat, rho, pi, tau){
  # A function to calculate power for a given design with a given expected
  #rate of retention (expretent)
  
  
  myvar <- TreatEffVar_open(expretent, totvar, m, desmat, rho, pi, tau)
  
  mypower <- pnorm( -qnorm(1-alpha/2) + sqrt(1/myvar)*effsize )
  return(mypower)
  
  
}

TreatEffVar_open_het <- function(expretent, Xmat, m, sigu2, sigv2, siguv =0, sig2E, sig2eta) {
    #expretent: expected propotion of subjects retained from one period to the next
    #Xmat: the design matrix
    #m: number if subjects in each cluster-period
    #sigu2: variance of cluster random intercepts
    #sigv2: variance of treatment effect slopes
    #siguv: covariance of random intercept and slope
    #sig2E: variance of subject-level errors
    #sig2eta: variance of the subject-level random effects (0 for cross-sectional data)
  
  sig2 <- sig2E/m
  
  T <- ncol(Xmat)
  K <- nrow(Xmat)
  
  
  Xvec <- as.vector(t(Xmat))
  stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
  Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
  
  #If type = 0 (HH or HG)
  #Variance matrix
  Vmat <- matrix(data=c(sigu2, siguv, siguv, sigv2), nrow=2, byrow=TRUE)
  
  Dtemp <- cbind(matrix(data=1, nrow=K, ncol=T), Xmat)
  #Generate the D matrix for all clusters
  # D <- kronecker(diag(1,2), Dvec)
  
  D<- blkDiag(z=array(as.vector(t(Dtemp)), c(T,2,K)))
  #blkDiag assumes square matrices so puffs up the matrices with zero columns.
  #need to remove these zero columns
  #To avoid removing those columns corresponding to clusters that are never
  #exposed, add a header row, which is later deleted.
  D<-rbind(rep(c(1,1,rep(0, T-2)), K),D)
  D<- D[,c(colSums(D, na.rm = TRUE)!=0)]
  D<-D[-1,]
  
  #Remove those rows of D which correspond to clusters not observed in 
  #particular periods
  D<- D[!is.na(Xvec),]
  
  #This is the only part that needs to accomodate the open cohort sampling structure
  varYbar <- D%*%kronecker(diag(1,K), Vmat)%*%t(D) +
   diag(sig2 + sig2eta*(1-expretent)/m, nrow(D)) +
    kronecker(diag(1,K), matrix(data=sig2eta*expretent/m, ncol(Xmat), ncol(Xmat)))

  #Variance of the treatment effect estimator is then given by:
  return(solve(t(Zmat)%*%solve(varYbar)%*%Zmat)[ncol(Zmat),ncol(Zmat)] )


}

SampleSize_open_TEhet <- function(expretent, effsize, power, alpha, Xmat, m, sigu2, sigv2, siguv =0, sig2E, sig2eta){
  #expretent is actually the expected proportion lost
  myvar <- TreatEffVar_open_het(expretent, Xmat, m, sigu2, sigv2, siguv, sig2E, sig2eta)
  clusperseq <- nrow(Xmat)*myvar*((qnorm(1-alpha/2) + qnorm(power))^2)/(effsize^2)
  
  
  return(clusperseq)
}

Power_open_TEhet <- function(expretent, effsize, alpha, Xmat, m, sigu2, sigv2, siguv, sig2E, sig2eta){
  # A function to calculate power for a given design with a given expected
  #rate of retention (expretent)
  
  myvar <- TreatEffVar_open_het(expretent, Xmat, m, sigu2, sigv2, siguv, sig2E, sig2eta)
  
  mypower <- pnorm( -qnorm(1-alpha/2) + sqrt(1/myvar)*effsize )
  return(mypower)
  
  
}


TreatEffVar_open_decays <- function(expretent, totvar, m, desmat, rho, pi, tau, clusdecay, subjdecay) {
  #expretent: expected propotion of subjects retained from one period to the next
  #Xmat: the design matrix
  #m: number if subjects in each cluster-period
  #clusdecay = 0 if no cluster decay; 1 if some cluster decay
  #subjdecay = 0 if no subj decay; 1 if some subj decay
  #icc=rho; cac= pi; iac=tau
  icc <- rho
  cac <- pi
  iac <- tau
  
  K <- nrow(desmat) 
  T <- ncol(desmat)
  Xvec <-  as.vector(t(desmat))
  
  #Generate the correct variance matrix for the set of 
  #cluster and subject decay possibilities
  if(clusdecay==0 && subjdecay==0){
    #Transforming the icc, cac, iac, and sd to variance components
    varclus <- (totvar)*(icc*cac)
    varCP <- (totvar)*(icc*(1-cac))
    varsubj <- (totvar)*(iac*(1-icc))
    varerror <- (totvar)*((1-iac)*(1-icc))
  
    Vi <- matrix(data=varclus + varsubj*(expretent)/m, nrow=ncol(desmat), ncol=ncol(desmat)) +
      diag(varCP + varerror/m + varsubj*(1-expretent)/m, nrow=ncol(desmat), ncol=ncol(desmat)) 
  }
  if(clusdecay==1 && subjdecay==0){
    sig2CP <- icc*totvar
    r <- cac 
    sig2E <- totvar*(1-iac)*(1-icc)
    sig2 <- sig2E/m
    sigindiv <- totvar*iac*(1-icc)/m
    sigindivtot <- sig2 + sigindiv
    #expretent matrix
    expmat <- matrix(data=expretent, nrow=T, ncol=T)
    diag(expmat) <- 1
    
    Vi <- sigindivtot*expmat*(diag(1-iac, T) + matrix(data=iac, nrow=T, ncol=T)) + 
      sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
  }
  if(clusdecay==1 && subjdecay==1){
    sig2CP <- icc*totvar
    r <- cac 
    sig2E <- totvar*(1-iac)*(1-icc)
    sig2 <- sig2E/m
    sigindiv <- totvar*iac*(1-icc)/m
    sigindivtot <- sig2 + sigindiv
    #expretent matrix
    expmat <- matrix(data=expretent, nrow=T, ncol=T)
    diag(expmat) <- 1
    
    Vi <- sigindivtot*expmat*(iac^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE))) + 
      sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
    
  }
  
  if(clusdecay==0 && subjdecay==1){
    sig2CP <- icc*totvar
    r <- cac 
    sig2E <- totvar*(1-iac)*(1-icc)
    sig2 <- sig2E/m
    sigindiv <- totvar*iac*(1-icc)/m
    sigindivtot <- sig2 + sigindiv
    #expretent matrix
    expmat <- matrix(data=expretent, nrow=T, ncol=T)
    diag(expmat) <- 1
    
    Vi <- sigindivtot*expmat*(iac^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE))) + 
      matrix(data=sig2CP*r, nrow=T, ncol=T) + diag((1-r)*sig2CP, T)
      
  }
  
  vartheta <-  1/(t(Xvec)%*%(diag(1,K)%x%solve(Vi))%*%Xvec -colSums(desmat)%*%solve(Vi)%*%(matrix(colSums(desmat),nrow=T, ncol=1))/K )
  return(vartheta)
  
}

Power_open_decays <- function(expretent, effsize, totvar, alpha, m, desmat, rho, pi, tau, clusdecay, subjdecay){
  # A function to calculate power for a given design with a given expected
  #rate of retention (expretent)
  
  
  myvar <- TreatEffVar_open_decays(expretent, totvar, m, desmat, rho, pi, tau, clusdecay, subjdecay)
  
  mypower <- pnorm( -qnorm(1-alpha/2) + sqrt(1/myvar)*effsize )
  return(mypower)
  
  
}

SampleSize_open_decays <- function(expretent, effsize, power, alpha, totvar, m, desmat, rho, pi, tau, clusdecay, subjdecay){
  #expretent is actually the expected proportion lost
  myvar <- TreatEffVar_open_decays(expretent, totvar, m, desmat, rho, pi, tau, clusdecay, subjdecay)
    
  clusperseq <- nrow(desmat)*myvar*((qnorm(1-alpha/2) + qnorm(power))^2)/(effsize^2)
  
  
  return(clusperseq)
}


INFORP_SampleSize_open <- function(p, effsize, power, alpha, totvar, m, desmat, rho, pi, tau, clusdecay, subjdecay){
  #p is the maximum number of periods subjects provide measurements for "in-for-p"
  
  
  myvar <- TreatEffVar_open_decays_inforp(p, totvar, m, desmat, rho, pi, tau, clusdecay, subjdecay)
  
  clusperseq <- nrow(desmat)*myvar*((qnorm(1-alpha/2) + qnorm(power))^2)/(effsize^2)
  
  
  return(clusperseq)
}

TreatEffVar_open_decays_inforp <- function(p, totvar, m, desmat, rho, pi, tau, clusdecay, subjdecay) {
  #p is the maximum number of periods subjects provide measurements for "in-for-p"
  #Xmat: the design matrix
  #m: number if subjects in each cluster-period
  #clusdecay = 0 if no cluster decay; 1 if some cluster decay
  #subjdecay = 0 if no subj decay; 1 if some subj decay
  #icc=rho; cac= pi; iac=tau
  icc <- rho
  cac <- pi
  iac <- tau
  
  K <- nrow(desmat) 
  T <- ncol(desmat)
  Xvec <-  as.vector(t(desmat))
  
  #Generate the correct variance matrix for the set of 
  #cluster and subject decay possibilities
  if(clusdecay==0 && subjdecay==0){
    #Transforming the icc, cac, iac, and sd to variance components
    varclus <- (totvar)*(icc*cac)
    varCP <- (totvar)*(icc*(1-cac))
    varsubj <- (totvar)*(iac*(1-icc))
    varerror <- (totvar)*((1-iac)*(1-icc))
    
    retentmat <-1- abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE))/p
    retentmat[which(retentmat < 0)] <- 0
    
    Vi <- diag(varCP + varerror/m, nrow=ncol(desmat), ncol=ncol(desmat)) +
         matrix(data=varclus, nrow=ncol(desmat), ncol=ncol(desmat)) +
         varsubj*retentmat/m
       
  }
  if(clusdecay==1 && subjdecay==0){
    sig2CP <- icc*totvar
    r <- cac 
    sig2E <- totvar*(1-iac)*(1-icc)
    sig2 <- sig2E/m
    sigindiv <- totvar*iac*(1-icc)/m
    sigindivtot <- sig2 + sigindiv
    #expretent matrix
    expmat <- matrix(data=expretent, nrow=T, ncol=T)
    diag(expmat) <- 1
    
    retentmat <-1- abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE))/p
    retentmat[which(retentmat < 0)] <- 0
    
    Vi <- sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE))) +
          sigindivtot*retentmat*(diag(1-iac, T) + matrix(data=iac, nrow=T, ncol=T))

  }
  if(clusdecay==1 && subjdecay==1){
    sig2CP <- icc*totvar
    r <- cac 
    sig2E <- totvar*(1-iac)*(1-icc)
    sig2 <- sig2E/m
    sigindiv <- totvar*iac*(1-icc)/m
    sigindivtot <- sig2 + sigindiv
    #expretent matrix
    expmat <- matrix(data=expretent, nrow=T, ncol=T)
    diag(expmat) <- 1
    retentmat <-1- abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE))/p
    retentmat[which(retentmat < 0)] <- 0
    
    Vi <- sigindivtot*retentmat*(iac^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE))) + 
      sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
    
  }
  
  if(clusdecay==0 && subjdecay==1){
    sig2CP <- icc*totvar
    r <- cac 
    sig2E <- totvar*(1-iac)*(1-icc)
    sig2 <- sig2E/m
    sigindiv <- totvar*iac*(1-icc)/m
    sigindivtot <- sig2 + sigindiv
    #expretent matrix
    expmat <- matrix(data=expretent, nrow=T, ncol=T)
    diag(expmat) <- 1
    retentmat <-1- abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE))/p
    retentmat[which(retentmat < 0)] <- 0
    
    Vi <- sigindivtot*retentmat*(iac^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE))) + 
      matrix(data=sig2CP*r, nrow=T, ncol=T) + diag((1-r)*sig2CP, T)
    
  }
  
  vartheta <-  1/(t(Xvec)%*%(diag(1,K)%x%solve(Vi))%*%Xvec -colSums(desmat)%*%solve(Vi)%*%(matrix(colSums(desmat),nrow=T, ncol=1))/K )
  return(vartheta)
  
}

Power_open_decaysINFORP <- function(p, effsize, totvar, alpha, m, desmat, rho, pi, tau, clusdecay, subjdecay){
  # A function to calculate power for a given design with a given expected
  #rate of retention (expretent)
  
  
  myvar <- TreatEffVar_open_decays_inforp(p, totvar, m, desmat, rho, pi, tau, clusdecay, subjdecay)
  
  mypower <- pnorm( -qnorm(1-alpha/2) + sqrt(1/myvar)*effsize )
  return(mypower)
  
  
}

#mydesmat <- rbind(SWdesmat(4),SWdesmat(4),SWdesmat(4),SWdesmat(4))
#Power_open(expretent = 1, effsize=2, totvar=25, alpha=0.05, 
#           m=10, desmat = mydesmat, rho=0.33, pi=0.9, tau=0.7)

#temp <- lapply(seq(0,1, 0.01), Power_open, effsize=2, totvar=25, alpha=0.05, 
#               m=10, desmat = mydesmat, rho=0.33, pi=0.9, tau=0.5)
#plot(seq(0,1, 0.01), temp)

#myvar <- TreatEffVar_open(1,25,10,SWdesmat(4),rho=0.33, pi=0.9, tau=0.7)

#mypow <- pnorm( -1.96 + sqrt(1/myvar)*2 )



#########################################################################


#SampleSize_open(expretent=1, effsize=2, totvar=25, power=0.8, alpha=0.05, m=10, desmat=SWdesmat(4),
#                rho=0.33, pi=0.9, tau=0.7)

#temp <- lapply(seq(0,1, 0.01), SampleSize_open, effsize=1, totvar=1, power=0.8, alpha=0.05, 
#               m=10, desmat=SWdesmat(4), rho=0.33, pi=0.9, tau=0.7)

#temp <- lapply(seq(0,1, 0.01), TreatEffVar_open, T=5, m=5, sig2E=0.95, sig2C=0.2, sig2eta=0.1, design=1)