# R Code by Zahra Aghahosseinalishirazi

### True Z and U and Mu , Update nu###
# rm(list = ls())
# library(FDRSeg)
# library(MASS)

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
#################################################################################################################
#################################################################################################################
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
##################################################################################################################
###################################################################################################################
# disp_para<- function (y, mu, n = sum(w), w, limit = 500, eps = .Machine$double.eps^0.25, trace = FALSE) {
#
#   score <- function(n, th, mu, y, w) sum(w * (digamma(th +y) - digamma(th) + log(th) + 1 - log(th + mu) - (y +th)/(mu + th)))
#   info <- function(n, th, mu, y, w) sum(w * (-trigamma(th +y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu +th)^2))
#   if (inherits(y, "lm")) {
#     mu <- y$fitted.values
#     y <- if (is.null(y$y))
#       mu + residuals(y)
#     else y$y
#   }
#   #if (missing(weights))
#   if (missing(w))
#
#     #  weights <- rep(1, length(y))
#     w <- rep(1, length(y))
#
#   #t0 <- n/sum(weights * (y/mu - 1)^2)
#   t0 <- n/sum(w * (y/mu - 1)^2)
#
#   it <- 0
#   del <- 1
#   if (trace)
#     message(sprintf("theta.ml: iter %d 'theta = %f'", it,
#                     signif(t0)), domain = NA)
#   while ((it <- it + 1) < limit && abs(del) > eps) {
#     t0 <- abs(t0)
#     #del <- score(n, t0, mu, y, weights)/(i <- info(n, t0, mu, y, weights))
#
#     del <- score(n, t0, mu, y, w)/(i <- info(n, t0, mu, y, w))
#     t0 <- t0 + del
#     if (trace)
#       message("theta.ml: iter", it, " theta =", signif(t0))
#   }
#   if (t0 < 0) {
#     t0 <- 0
#     warning("estimate truncated at zero")
#     attr(t0, "warn") <- gettext("estimate truncated at zero")
#   }
#   if (it == limit) {
#     warning("iteration limit reached")
#     attr(t0, "warn") <- gettext("iteration limit reached")
#   }
#   attr(t0, "SE") <- sqrt(1/i)
#   t0
# }

#cell_matrix
#as.vector(cell_matrix)
#######################################################################################################
#######################################################################################################
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
    mu_mat[,k]<- rep(MU[k,],N)
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

########################################################
##### THE EM Function for Mixture of MIZNB####
###################################################################################################
ZINB_mix_K_comp<- function(Y, P, Phi, mu, dispresiopn_parameter, eps=1e-5, maxit=1000){
  dl <- 1 + eps
  iter<-0
  ll <- rep(0, maxit+1)
  G<- ncol(Y)
  N<- nrow(Y)
  K<- length(Phi)

  nu<- dispresiopn_parameter

  bk<- array(NA,c(N,G,K))
  for(k in 1:K){
    for(n in 1:N){
      for(g in 1:G){
        if(Y[n,g]==0) bk[n,g,k]<- (Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])))


        else bk[n,g,k]<- (1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])
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
        if(Y[n,g]==0) d[n,g,k]<- P[k]*(Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])))

        else d[n,g,k]<- 0

      }
    }
  }
  # print(d)




  while (abs(dl)>eps & iter<maxit) { ## Here I use abs(dl) rather than dl#
    # cat(abs(dl), "\t", iter, "\n")
    iter<-iter+1
    ll[iter] <- l


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



    # #Z_nk#
    # Z_nk<- matrix(0,nrow=N,ncol=length(P))
    # for(i in 1:N){
    #   if(clust[i]==1)Z_nk[i,1]<- 1
    #   else Z_nk[i,2]<- 1
    # }
    #Z
    ##################################################
    ##################################################
    ## U_ng##
    # U_ng<- matrix(NA,nrow=N,ncol=G)
    # for(n in 1:N){
    #   for(g in 1:G){
    #     if(Y[n,g]==0) U_ng[n,g]<- 1
    #     else U_ng[n,g]<-0
    #   }
    # }


    #print(P)



    # U_ng<- U_ng
    #print(U_ng)

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

    #print(U_ngk)

    # for(k in 1:K){
    #   for(n in 1:N){
    #     for(g in 1:G){
    #       if(U_ngk[n,g,k]>0.7) U_ngk[n,g,k]<- 1
    #       else U_ngk[n,g,k]<- 0
    #
    #     }
    #   }
    # }
    #

    #print(U_ngk)
    #print(1-U_ngk)
    #Phi<- Phi

    for(k in 1:K){
      ZU[,,k]<- Z_nk[,k]*U_ngk[,,k]
    }
    #print(ZU)

    #print(apply(ZU,3,sum))

    P<- apply(Z_nk,2,mean)
    #print(P)
    Phi<- (apply(ZU,3,sum))/(G*apply(Z_nk,2,sum))
    #print(Phi)
    ###############################################################################################
    ################################################################################################

    mult_Z_u<- array(NA,c(N,G,K))
    mult_Z_oneminus_U<- array(NA,c(N,G,K))
    mult_Z_oneminus_U_Y<- array(NA,c(N,G,K))



    for(k in 1:K){
      mult_Z_u[,,k]<- Z_nk[,k]*U_ngk[,,k]
      mult_Z_oneminus_U[,,k]<- Z_nk[,k]*(1-U_ngk[,,k])
      mult_Z_oneminus_U_Y[,,k]<- Z_nk[,k]*(1-U_ngk[,,k])*Y

    }


    mu<- t(apply(mult_Z_oneminus_U_Y,c(2,3),sum)/apply(mult_Z_oneminus_U,c(2,3),sum))
    #print(mu)

    # mu<- mu
    ##########################################################################################
    ##########################################################################################
    ### To estimate size_MLE###
    # result<- To_vec_transfer(Y=cell_matrix,Z=Z_nk,U=U_ngk,MU=mu)
    result<- To_vec_transfer(Y, Z_nk, U_ngk, mu)
    #print(result$WW)

    size_MLE<- rep(NA,K)

    for(k in 1:K){
      size_MLE[k]<- disp_para(y=result$Y_vec, mu=result$MU_vec[,k], n = sum(result$WW[,k]), w=result$WW[,k], limit = 1000, eps = .Machine$double.eps^0.25, trace = FALSE)

    }
    nu<- size_MLE

    #nu<- theta_est
    #nu<- size_MD
     #print(nu)

    # print(nu)
    ################################################################################################################
    ################################################################################################################
    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          if(Y[n,g]==0) bk[n,g,k]<- (Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])))
          else bk[n,g,k]<- (1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])
        }
      }
    }

    ############################################################
    for(n in 1:N){
      for(g in 1:G){
        for(k in 1:K){
          if(Y[n,g]==0) d[n,g,k]<- P[k]*(Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])))
          else d[n,g,k]<- 0

        }
      }
    }

    #########################################################
    for(k in 1:K){
      b[,,k]<- P[k]*bk[,,k]
    }
    bp_sum<- apply(b,c(1,2),sum)
    #print(bp_sum)
    oldl<- l

    l<- sum(rowSums(log(bp_sum))) #Observed data log-likelihood
    #print(l)

    ################################################################

    #print(l)
    dl<- l-oldl
  }
  cat("number of iterations=", iter, "\n")
  iter <- iter+1
  ll[iter] <- l
  postprobs <- Z_nk
  colnames(postprobs) <- c(paste("comp", ".", 1:K, sep = ""))


  # out <- list(P=P,Phi=Phi,MU=mu,size_est=nu,a_hat=1/nu,postprobs=postprobs,U=U_ngk,
  #             loglik=l,all.loglik=ll[1:iter],restart=0,ft="ZINBmix")
  # class(out)<- "mixEM"
  # out

  return(list(
    prob   = P,
    phi    = Phi,
    mu     = mu,
    size   = nu,
    U      = U_ngk,
    Z      = Z_nk,
    loglik = l
  ))

}
###################################################################################################
###################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
