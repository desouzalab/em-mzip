# R Code by Zahra Aghahosseinalishirazi

### True Z and U and Mu , Update nu###
# rm(list = ls())
# library(FDRSeg)
# library(MASS)
# library(pheatmap)

###################################################################################################
### Generate data from Zero-inflated Negative Binomial Distribution###
ZINB_r_U<- function(n,prob,theta,disp_param){
  Z<- rbinom(size=1,n=n,p=prob)
  Y<- NULL
  for(i in 1:n){
    if(Z[i]==1) Y<- c(Y,0)
    else Y<- c(Y,rnbinom(1,mu = theta,size=disp_param))
  }
  return(output=list(Y=Y,U=Z))
}
#ZINB_r_U(50,0.5,4,10)
#################################################################################################################################
#################################################################################################################################
disp_para<- function (y, mu, n = sum(w), w, limit = 1000, eps = .Machine$double.eps^0.5, trace = FALSE) {

  score <- function(n, th, mu, y, w) sum(w * (digamma(th + y) - digamma(th) + log(th) + 1 - log(th + mu) - (y + th)/(mu + th)))

  info <- function(n, th, mu, y, w) sum(w * (-trigamma(th + y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu + th)^2))

  # uses a new formula to avoid numerical errors
  score.trick <- function(n, th, mu, y, w, eps = 1e-5) {
    sum(w * ((digamma(th + y + eps) - digamma(th - eps)) / (2 * eps) + log(th) + 1 - log(th + mu) - (y + th)/(mu + th)))
  }

  info.trick <- function(n, th, mu, y, w, eps = 1e-5) {
    sum(w * ((-trigamma(th + y + eps) + trigamma(th - eps)) / (2 * eps) - 1/th + 2/(mu + th) - (y + th)/(mu + th)^2))
  }


  if (inherits(y, "lm")) {
    mu <- y$fitted.values
    y <- if (is.null(y$y))
      mu + residuals(y)
    else y$y
  }
  #if (missing(weights))
  if (missing(w))

    #  weights <- rep(1, length(y))
    w <- rep(1, length(y))

  #t0 <- n/sum(weights * (y/mu - 1)^2)
  t0 <- n/sum(w * (y/mu - 1)^2)

  it <- 0
  del <- 1
  if (trace)
    message(sprintf("theta.ml: iter %d 'theta = %f'", it,
                    signif(t0)), domain = NA)
  while ((it <- it + 1) < limit && abs(del) > eps) {
    t0 <- abs(t0)
    #del <- score(n, t0, mu, y, weights)/(i <- info(n, t0, mu, y, weights))

    del <- score(n, t0, mu, y, w)/(i <- info(n, t0, mu, y, w))

    # check using different formula, or replaces del with 0
    if (is.nan(del)) {
      del.eps <- 1e-5
      del <- score.trick(n, t0, mu, y, w, eps = del.eps) /
        (i <- info.trick(n, t0, mu, y, w, eps = del.eps))
      # print(del)

      if (is.nan(del)) {
        del <- 0
        warning("NaN obtained, replacing del with 0")
      } else {
        warning(paste("NaN obtained, used new formula with eps =", del.eps))
      }
    }

    t0 <- t0 + del

    if (trace)
      message("theta.ml: iter", it, " theta =", signif(t0))
  }
  if (t0 < 0) {
    t0 <- 0
    warning("estimate truncated at zero")
    attr(t0, "warn") <- gettext("estimate truncated at zero")
  }
  if (it == limit) {
    warning("iteration limit reached")
    attr(t0, "warn") <- gettext("iteration limit reached")
  }
  attr(t0, "SE") <- sqrt(1/i)
  t0
}
####################################################################################################################
####################################################################################################################
To_vec_transfer<- function(Y,Z,U,MU){
  N<- nrow(Y)
  G<- ncol(Y)
  K<- ncol(Z)
  UU<- 1-U

  #print(UU)
  #print(Z)

  vec_y<- as.vector(t(Y))
  ZU<- array(NA,c(N,G,K))
  for(k in 1:K){
    ZU[,,k]<- Z[,k]*UU[,,k]
  }
  #print(ZU)

  mu_mat<- matrix(NA,nrow=N*G,ncol=K)
  for(k in 1:K){
    mu_mat[,k]<- as.vector(t(MU[,,k]))
  }

  weight_mat<- matrix(NA,nrow=N*G,ncol=K)
  for(k in 1:K){
    weight_mat[,k]<- as.vector(t(ZU[,,k]))
  }

  #print(weight_mat)
  #
  # for(i in 1:nrow(weight_mat)){
  #   for(k in 1:K){
  #     if (weight_mat[i,k]<0){
  #       print(i)
  #       print(k)
  #       print(weight_mat[i,k])
  #       weight_mat[i,k]<-0
  #       warning("the element is less than zero")}
  #
  #   }
  #
  # }
  #
  #print(ZU)
  out<- list(Y_vec=vec_y,WW=weight_mat,MU_vec=mu_mat)
  return(out)

}



#################################################################################
##########################################################################################
NR_ZINB<- function(Y,totals,beta_0,rho,Z_nk,U_ng,overdisp_para,criterion=1e-5,maxiter=500){
  N<- nrow(Y)
  G<- ncol(Y)
  K<- dim(Z_nk)[2]
  nu<- overdisp_para

  a<- 1/nu
  ### Initialization###
  mu_ng<- array(data=NA,c(N,G,K))
  diff_Y_mu<- array(data=NA,c(N,G,K))
  U_diff_Y_mu<- array(data=NA,c(N,G,K))
  U_mu<- array(data=NA,c(N,G,K))

  Z_diff<- array(data=NA,c(N,G,K))
  Z_mu<- array(data=NA,c(N,G,K))
  theta<- numeric((G*(K+1)))

  ### Assign values to theta vector###
  theta[1:G]<- beta_0

  for(k in 1:K){
    theta[((k*G)+1):((k*G)+G)]<- rho[k,]
  }

  it<- 0
  diff_theta<- 1000
  while(it<maxiter & diff_theta>criterion){
    #while(it<maxiter & all(diff_theta>criterion)){

    it<- it+1
    #print(it)
    ### mu###
    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          mu_ng[n,g,k]<- exp(log(totals[n])+beta_0[g]+rho[k,g])

        }
      }
    }
    #mu_ng
    ### Y-mu###
    for(k in 1:K){
      diff_Y_mu[,,k]<- (Y-mu_ng[,,k])/(1+(a[k]*mu_ng[,,k]))
    }

    for(k in 1:K){
      U_diff_Y_mu[,,k]<- (1-U_ng[,,k])*diff_Y_mu[,,k]
      U_mu[,,k]<- (1-U_ng[,,k])*(mu_ng[,,k]*(1+a[k]*Y))/((1+a[k]*mu_ng[,,k])^2)
    }

    U_diff_Y_mu
    U_mu


    ### Z*(Y-mu)& Z*mu###
    for(k in 1:K){
      Z_diff[,,k]<- Z_nk[,k]*U_diff_Y_mu[,,k]
      Z_mu[,,k]<- Z_nk[,k]*U_mu[,,k]
    }

    #Z_diff
    #Z_mu

    ### derivatives of beta_0###
    db_0<- ddb_0<-  numeric(G)

    for(g in 1:G){
      db_0[g]<- sum(Z_diff[,g,])
      ddb_0[g]<- sum(Z_mu[,g,])
    }

    ### Derivatives of rho's###

    drho<- ddrho<- matrix(data=NA,nrow=K,ncol=G)

    for(k in 1:K){
      for(g in 1:G){
        drho[k,g]<- sum(Z_diff[,g,k])
        ddrho[k,g]<- sum(Z_mu[,g,k])
      }
    }
    ### Gradient Vecotr###
    grad<- numeric((G*(K+1)))

    grad[1:G]<- db_0
    for(k in 1:K){
      grad[((k*G)+1):((k*G)+G)]<- drho[k,]
    }
    ### Hessian Matrix###
    hessian<- matrix(0,ncol=(G*(K+1)),nrow=(G*(K+1)))
    hessian[1:G,1:G]<- hessian[1:G,1:G]-diag(ddb_0)

    for(k in 1:K){
      hessian[((k*G)+1):((k*G)+G),((k*G)+1):((k*G)+G)]<- hessian[((k*G)+1):((k*G)+G),((k*G)+1):((k*G)+G)]-diag(ddrho[k,])
      hessian[((k*G)+1):((k*G)+G),1:G]<- hessian[((k*G)+1):((k*G)+G),1:G]-diag(ddrho[k,])
      hessian[1:G,((k*G)+1):((k*G)+G)]<- hessian[1:G,((k*G)+1):((k*G)+G)]-diag(ddrho[k,])
    }

    ### Updating theta using Newton-Raphson Algorithm###
    theta_new<- theta-qr.coef(qr(hessian,tol=1e-300), grad) ### tol=1e-300

    for(i in 1:length(theta_new)){
      if(is.na(theta_new[i])) theta_new[i]<- theta[i]

    }

    diff_theta<- mean(abs(theta_new-theta))
    #diff_theta<- abs(theta_new-theta)

    #print(theta)
    #print(theta_new)

    theta<- theta_new

    #beta_0<- theta_new[1:G]
    #for(k in 1:K){
    # rho[k,]<- theta_new[((k*G)+1):((k*G)+G)]
    #}

    beta_0<- theta[1:G]
    for(k in 1:K){
      rho[k,]<- theta[((k*G)+1):((k*G)+G)]
    }



    #out<- list(beta_0=beta_0,rho=rho,thtea=theta,iteration=it)


  }

  out<- list(beta_0=beta_0,rho=rho,thtea=theta,iteration=it)

  return(out)


  #  }


}

#########################################################################################################################
#########################################################################################################################
ZINB_Mix_EM <- function(Y, totals, beta_0g, rho_gk, dispresiopn_parameter, P, Phi, eps=1e-5, maxiter_EM=1000){
  N<- nrow(Y)
  G<- ncol(Y)
  K<- length(P)
  dl<- 1+eps
  iter<- 0
  ll<- rep(0,maxiter_EM+1)

  dl <- 1 + eps


  nu<- dispresiopn_parameter

  beta_0g_cur<- beta_0g
  rho_cur<- rho_gk

  Z_nk<- matrix(data=NA,nrow=N,ncol=K)
  U_ngk<- array(data=NA,c(N,G,K))

  ZU<- array(data=NA,c(N,G,K))

  thresh<- -744


  mu<- array(data=NA,c(N,G,K))
  #print(beta_0)
  #print(rho)
  for(n in 1:N){
    for(g in 1:G){
      for(k in 1:K){
        mu[n,g,k]<- exp(log(totals[n])+beta_0g[g]+rho_gk[k,g])
      }
    }
  }
  #print(mu)



  bk<- array(NA,c(N,G,K))
  for(k in 1:K){
    for(n in 1:N){
      for(g in 1:G){
        if(Y[n,g]==0) bk[n,g,k]<- (Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[n,g,k],size=nu[k])))


        else bk[n,g,k]<- (1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[n,g,k],size=nu[k])
      }
    }
  }


  b<- array(NA,c(N,G,K))

  for(k in 1:K)
  {b[,,k]<- P[k]*bk[,,k]}

  bp_sum<- apply(b,c(1,2),sum)
  l<- sum(rowSums(log(bp_sum))) #Observed data log-likelihood
  Z_nk<- matrix(NA,ncol=K,nrow=N)
  thresh<- -744


  U_ngk<- array(data=NA,c(N,G,K))
  ZU<- array(data=NA,c(N,G,K))



  d<- array(data=NA,c(N,G,K))

  for(n in 1:N){
    for(g in 1:G){
      for(k in 1:K){
        if(Y[n,g]==0) d[n,g,k]<- P[k]*(Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[n,g,k],size=nu[k])))

        else d[n,g,k]<- 0

      }
    }
  }
  # print(d)


  ##############################################################################
  ##############################################################################
  #ratio<- 1000
  #oldl<- -1000000



  while(abs(dl)>eps & iter<maxiter_EM){
    #while(iter<maxiter_EM){

    iter<- iter+1
    ll[iter]<- l
    #pi_Update <- colSums(Zhat)/length(Y)

    ##E-step##
    ###Estimate Z_nk###
    Z_nk<- matrix(data=log(P),ncol=K,nrow=N,byrow=T)

    for(k in 1:K){
      for(g in 1:G){
        Z_nk[,k]<- Z_nk[,k]+(log(bk[,g,k]))
      }
    }

    if (K>1){
      v1<-which(apply(Z_nk,1,max)< thresh)
      v3<-1:N
      len<-length(v1)
      if(len>0){
        v2<-apply(array(Z_nk[v1,],dim=c(len,K)),1,order)[K,]
        ddd<-cbind(v1,v2)
        Z_nk[v1,]<- 0
        Z_nk[ddd]<- 1
        v3<- -v1
      }
      Z_nk[v3,]<-exp(Z_nk[v3,])}
    Z_nk<-Z_nk/apply(Z_nk,1,sum)

    epsilon <- 1e-10
    sl<-length(Z_nk[Z_nk < epsilon])
    bl<-length(Z_nk[Z_nk > 1-epsilon])
    Z_nk[Z_nk<epsilon]<-rep(epsilon,sl)
    Z_nk[Z_nk>1-epsilon]<-rep(1-epsilon,bl)
    Z_nk<-Z_nk/apply(Z_nk,1,sum)
    #print(Z_nk)
    ###############################################

    U_ngk<- array(data=NA,c(N,G,K))

    #denom<- apply(d,c(1,2),sum)
    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          #print(P[k]*Phi[k])

          if(Y[n,g]==0) U_ngk[n,g,k]<- (P[k]*Phi[k])/d[n,g,k]
          else U_ngk[n,g,k]<- 0

        }
      }
    }



    #####################################################
    ### M-step###
    for(k in 1:K){
      ZU[,,k]<- Z_nk[,k]*U_ngk[,,k]
    }
    ##Update P###
    P<- apply(Z_nk,2,mean)


    #print(P)
    #return(P)

    ### Update Phi###
    Phi<- (apply(ZU,3,sum))/(G*apply(Z_nk,2,sum))
    #print(Phi_hat)
    #return(P)
    nu<- nu
    ##########################################################################################
    ### To estimate size_MLE###
    result<- To_vec_transfer(Y=cell_matrix,Z=Z_nk,U=U_ngk,MU=mu)
    #print(result$WW)
    #print(result$MU_vec)

    size_MLE<- rep(NA,K)

    for(k in 1:K){
      size_MLE[k]<- disp_para(y=result$Y_vec, mu=result$MU_vec[,k], n = sum(result$WW[,k]), w=result$WW[,k], limit = 1000, eps = .Machine$double.eps^0.25, trace = FALSE)

    }
    nu<- size_MLE



    # print(nu)
    ################################################################################################################
    ### Updaet Parameters of lambda: rho_gk, Beta_0g###

    ### NR POisson##
    #Newton_Update<- NR_Poisson(Y=Y,T=T,beta_0=beta_0g_cur,rho=rho_cur,Z_nk=Z_nk,K=K,criterion=1e-4,maxiter=500)

    #Newton_Update<- NR_Poisson_G(Y=Y,T=T,beta_0=beta_0g_cur,rho=rho_cur,Z_nk=Z_nk,K=K,criterion=1e-5,maxiter=500)

    ### NR ZIP##
    #Newton_Update<- NR_Poisson_G_ZIP(Y=Y,T=T,beta_0=beta_0g_cur,rho=rho_cur,Z_nk=Z_nk,U_ng=U_ng,K,criterion=1e-5,maxiter=500)
    #Newton_Update<- NR_Poisson_ZIP(Y=Y,T=T,beta_0=beta_0g_cur,rho=rho_cur,Z_nk=Z_nk,U_ng=U_ng,K,criterion=1e-5,maxiter=500)

    #Newton_Update<- NR_Poisson_ZIP_1X(Y=cell_matrix,T=T,beta_0=beta_0g_cur,rho=rho_cur,beta_pg=beta_pg_cur,x=X,Z_nk=Z_nk,U_ng=U_ngk,criterion=1e-4,maxiter=500)


    #Newton_Update<- NR_G_ZINB(Y=cell_matrix,totals=T,beta_0=beta_0g_cur,rho=rho_cur,Z_nk=Z_nk,U_ng=U_ngk,overdisp_para=nu,criterion=1e-5,maxiter=500)
    Newton_Update<- NR_ZINB(Y=cell_matrix,totals=T,beta_0=beta_0g_cur,rho=rho_cur,Z_nk=Z_nk,U_ng=U_ngk,overdisp_para=nu,criterion=1e-5,maxiter=500)


    beta_0g<- Newton_Update$beta_0
    rho_gk<- Newton_Update$rho
    ##############################################################################################################
    #######################
    for(n in 1:N){
      for(g in 1:G){
        for(k in 1:K){
          mu[n,g,k]<- exp(log(totals[n])+beta_0g[g]+rho_gk[k,g])
        }
      }
    }




    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          if(Y[n,g]==0) bk[n,g,k]<- (Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[n,g,k],size=nu[k])))


          else bk[n,g,k]<- (1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[n,g,k],size=nu[k])
        }
      }
    }






    for(n in 1:N){
      for(g in 1:G){
        for(k in 1:K){
          if(Y[n,g]==0) d[n,g,k]<- P[k]*(Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[n,g,k],size=nu[k])))

          else d[n,g,k]<- 0

        }
      }
    }
    # print(d)

    beta_0g_cur<- beta_0g
    rho_cur<- rho_gk


    for(k in 1:K){
      b[,,k]<- P[k]*bk[,,k]
    }
    bp_sum<- apply(b,c(1,2),sum)
    #print(bp_sum)
    oldl<- l

    l<- sum(rowSums(log(bp_sum))) #Observed data log-likelihood
    #print(l)

    ########################

    ################################################################

    #print(l)
    dl<- l-oldl
  }
  cat("number of iterations=", iter, "\n")
  iter <- iter+1
  ll[iter] <- l
  postprobs <- Z_nk
  colnames(postprobs) <- c(paste("comp", ".", 1:K, sep = ""))


  # out <- list(P=P,phi=Phi,beta_0g=beta_0g_cur,rho_gk=rho_cur,size_est=nu,a_hat=1/nu,no_iter=iter,Postprobs=postprobs,U=U_ngk,
  #               loglik=l,all.loglik=ll[1:iter],restart=0,ft="ZINBmix")
  #   class(out)<- "mixEM"
  #   out

  return(list(
    prob   = P,
    phi    = Phi,
    rho    = rho_cur,
    beta0  = beta_0g_cur,
    size   = nu,
    U      = U_ngk,
    Z      = Z_nk,
    loglik = l
  ))

}
