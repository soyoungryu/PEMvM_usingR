##################################################################################
#Author: So Young Ryu 		
#Contact: soyoungr@unr.edu  
#PEMvM codes are based on the published codes in PEMM R package by Chen et al. 2014
#PEMvM handles varying non-ignorable missing mechanisms across experiemnts
###################################################################################
#Parameter
###########
#X: Matrix with dimension n X p; Protein abundances
#phi: a vector for the missing data mechanisms. The length of vector is the number of experiment. (beta in Ryu 2017+ manuscript). The same value needs to be used for the same group. 
#lambda:  the tuning parameter in the Inverse-Wishart penalty function (default 5)
#K:  the second tuning parameter in the Inverse-Wishart penalty function (detault 5) 
#tol: the tolerance to assess the convergence of the iterative algorithm (detaul tol=0.001)
#maxIter: the maximum number of iterations allows. (default=100)
###################################################################################
#Value
############
#mu: the estimated mean
#Sigma: the estimated covariance
#Xhat: imputed missing value
#converge: whether the algorithm was converge or not 
#lik: the log-liklihood of the last iteration
#####################################################################



PEMvM <- function(X, phi, lambda=NULL, K=NULL,   tol=0.001, maxIter=100){
	delta=NULL; pos=NULL
	if (length(phi)==1){
		warning("Use the same phi for all experiments!")
		phi=rep(phi, nrow(X))
	}
	get.llik.mPEMM <- function(X, mu, S, phi, phi0){
	    if (length(which(X<0))==0) pos<- TRUE else pos<- FALSE
	    
	    llik <- 0
	    for (j in 1:nrow(X)){
	      Xi <- X[j,]
	      idxi <- which(!is.na(Xi))
	      Xi <- Xi[idxi]
	      Si <- as.matrix(S[idxi,idxi])
	      Sii <- my.solve(Si)
	      det.Si=det(Si)
	      #if(det.Si<=0){det.Si=0.0001}
	      oo <-    - log(det.Si)  -(Xi-mu[idxi])%*%Sii%*%(Xi-mu[idxi]) 
	      oo <- 0.5*oo
	      nmis <- length(X[j,])-length(idxi)
	      if (phi[j]==0){
	          vv <- 0
	      } else {
	         if (pos){
	
		  vv1=sum(log(1-exp(-phi0-phi[j]*Xi)))
		  vv2 <- -phi0*nmis
		  if (nmis==0){
		     vv3 <- 0
		  } else {
		      Smi <- as.matrix(S[-idxi,-idxi])
		      mu.mis <- mu[-idxi]+matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%(Xi-mu[idxi])
		      S.mis <- Smi - matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%matrix(S[idxi,-idxi],ncol=nmis)
		      vv3<- -sum(mu.mis)*phi[j]+0.5*sum(S.mis)*phi[j]^2  
		  }
	          vv <- vv1+vv2+vv3
	        } else {
	            vv1=sum(log(1-exp(-phi0-phi[j]*Xi^2)))
		    vv2 <- -phi0*nmis
		    if (nmis==0){
		  	     vv3 <- 0
		    } else {
		  	  Smi <- as.matrix(S[-idxi,-idxi])
		  	  mu.mis <- mu[-idxi]+matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%(Xi-mu[idxi])
		  	  S.mis <- Smi - matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%matrix(S[idxi,-idxi],ncol=nmis)
		  	  Smis.inv <- my.solve(S.mis)
		  	  Smii <- Smis.inv
		  	  diag(Smis.inv) <- diag(Smis.inv)+2*phi[j]
		  	  A <- my.solve(Smis.inv)
	  	  
		  	  vv3<- 0.5*(log(det(A)) -log(det(S.mis)) + matrix(mu.mis,nrow=1)%*%(Smii%*%A%*%Smii-Smii)%*%matrix(mu.mis,ncol=1) ) 
		    }
	            vv <- vv1+vv2+vv3
	
	        }
	      }
	      llik <- llik+ oo+vv
	    }  
    
	    pllik <- llik
	    return(llik) 
	}

	my.solve <- function(X){
	   if (!is.matrix(X))  X <- matrix(X, nrow=sqrt(length(X)))
	   ss <- svd(X)
	   Xinv <- ss$u%*%diag(1/ss$d, nrow=nrow(X), ncol=nrow(X))%*%t(ss$v)
	   return(Xinv)
	}

	gets1 <- function(sigma, phi){
	  p <- nrow(sigma)
	  s1 <- my.solve(sigma)+diag(2*phi, p, p)
	  return(s1)
	}

	get.bg <- function(sigma, mu, phi){
	    p <- length(mu)
	    s1 <- gets1(sigma, phi)
	    s1inv <- my.solve(s1)
    
	    ccen <- my.solve(sigma)
	    A <- s1inv%*%ccen
	    beta <- A%*%mu

	    gamma <- s1inv 
    
	    return(list(beta=beta, gamma=gamma))
	}

     	find.lambda <- function(Sigma,N,p,K, delta=delta){
	    ffL <- function(lambda, Sigma, N, p, K){
	        Sigma.new <- N/(N+K)*Sigma*(N-1)/N + lambda/(N+K)*diag(1, p, p)
	        return(abs(min(as.double(eigen(N*Sigma.new)$value))))
	    }

	    Sigma2 = Sigma
	    while (!is.double(eigen(N*Sigma2)$value)){
	        delta = delta+1
	        Sigma2 = N/(N+K)*Sigma*(N-1)/N + delta/(N+K)*diag(1, p, p)
	    }
	    Sigma=Sigma2
	
	    oo <- -min(as.double(eigen(N*Sigma)$value))
	    if (oo>0){
	        lambda <- optimize(ffL, lower=0,upper=oo, Sigma=Sigma, N=N, p=p, K=K)$minimum+delta
	    } else {
	        lambda <- delta
	    }
	    return(lambda)
       }

       if (is.null(pos)){
         if (length(which(X<0))==0) pos<- TRUE else pos<- FALSE
       }

          if (sum(abs(phi))==0) {
           phi0 <- -log(mean(is.na(X))) 
       } else {
           phi0 <- 0
       }
       
       p <- ncol(X)
       N <- nrow(X)
       if (is.null(K)){
         K=5
       } 
       X.hat <- X
       
       ## initial estimates
       mu.new <- matrix(colMeans(X,na.rm=T),ncol=1)
       Sigma.new <- cov(X, use="pairwise.complete")
       Sigma.new[is.na(Sigma.new)] <- 0  

       diff <- 999
       iter <- 0
       if (is.null(delta)) {
         delta=5
        }
         Lambda <- find.lambda(Sigma.new,N=N,p=p,K=K, delta=delta)
       Sigma.new <- N/(N+K)*Sigma.new*(N-1)/N + Lambda/(N+K)*diag(1, p, p)
       illik <- 999
       while(iter<maxIter & diff>tol){
         iter <- iter+1
         mu <- mu.new
         Sigma <- Sigma.new
         
         cov.add <- matrix(0, p, p)
         for (i in 1:nrow(X)){
           ii <- which(is.na(X[i,]))
           if (length(ii)>=1){
             Soo <- as.matrix(Sigma[-ii,-ii])
             pi <- nrow(Soo)
             mu.mis <-  mu[ii]+Sigma[ii,-ii]%*%my.solve(Soo)%*%(X[i,-ii]-mu[-ii]) 
             mu.mis <- matrix(mu.mis,ncol=1)
             cov.mis <- Sigma[ii,ii] - Sigma[ii,-ii]%*%my.solve(Soo)%*%Sigma[-ii, ii]
             
             if (phi[i]!=0 & pos==TRUE){

                 X.hat[i,ii]<- mu.mis - phi[i]*cov.mis%*%matrix(1,nrow=length(mu.mis),ncol=1)
                 cov.ii <- cov.mis
                
             } else if (phi[i]!=0 & pos==FALSE){

                 oo <- get.bg(cov.mis, mu.mis, phi[i])             
                 X.hat[i,ii] <- oo$beta
                 cov.ii <- oo$gamma

             } else if (phi[i]==0) {

                 X.hat[i,ii] <- mu.mis
                 cov.ii <- cov.mis 

             }     
             cov.add[ii,ii] <- cov.add[ii,ii] + cov.ii
           }
         }          
         
         mu.new <- colMeans(X.hat)
         Sig <- cov(X.hat)*(N-1)/N+cov.add/N
         if (is.null(lambda)) {
	      Lambda <- find.lambda(Sig,N=N,p=p,K=K,delta=delta)
	 } else {
	      Lambda <- lambda
         }
         Sigma.new <- N/(N+K)*(Sig)  + Lambda/(N+K)*diag(1, p, p)

		 det.Sigma.new=det(Sigma.new)
		 #if(det.Sigma.new <= 0){det.Sigma.new=0.0001}
         illik <- c(illik, get.llik.mPEMM(X, mu.new, Sigma.new,phi=phi,phi0=phi0) - sum(diag(Lambda*my.solve(Sigma.new))) -K*log(det.Sigma.new) ) 
         diff <- abs(illik[iter+1]-illik[iter])/abs(illik[iter])
         #cat(diff, "\n")
             if (is.na(diff)) {  	
            diff <- -1  
            warning("The algorithm does not converge!")
         }
       }
       if (iter==maxIter) warning("The algorithm does not converge!")   
       if (iter<maxIter & diff>0 & diff<tol) {converge=TRUE
       	}else{converge=FALSE}
       	 
       return(list(mu=mu, Sigma=Sigma, Xhat=X.hat, converge=converge, lik=illik))
}

