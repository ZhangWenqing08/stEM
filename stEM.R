########################################################################################################
#       estimation algorithm
#######################################
# Y:                   a # of subjects x # of blocks matrix of item responses
# BID:                 a # of statements x 3 matrix of item information; Columns are "Block", "Item" and "Dimensions"
# positive:            a logical vector indicating whether each statement is positive directional or not 
# M:                   # of batch (see page 3 of the manual)
# B:                   # of iterations in each batch
# a:                   initial alpha parameters in equation 2; a vector with length = # of statements
# d:                   initial beta parameters in equation 2; a vector with length = # of statements
# item.par:            initial parameters for a and d in a data frame
# sigma:               initial sigma parameters; a matrix with # of dimensions x # of dimensions
# theta:               initial theta parameters; a matrix with # of subjects x # of dimensions
# fix.sigma:           logical; TRUE if sigma is not estimated
# burnin.maxitr:       max burnin allowed
# maxitr:              max iterations allowed
# eps1:                see the manual
# eps2:                see the manual
# frac1 and frac2:     cutoffs for calculating Geweke z; see the manual
##########################################################################################################
stEM <- function(Y,BID,positive=rep(TRUE,nrow(BID)),M=10,B=20,a=NULL,d=NULL,
                 item.par=NULL,sigma=NULL,theta=NULL,fix.sigma = FALSE,burnin.maxitr=40,
                 maxitr=500,eps1=1.5,eps2=0.4,frac1=.2,frac2=.5){
  
  ##################################
  #
  #
  #
  #
  ##################################
  
  required.packages <- c("armspp", "doParallel","foreach","doRNG","coda","mvnfast","lvmcomp")
  to.be.installed <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(to.be.installed)) install.packages(to.be.installed)
  stopifnot(sapply(required.packages,require,character.only=TRUE))
  s1 <- Sys.time()
  #number of dimensions
  D <- length(unique(BID$Dim))
  N <- nrow(Y) #number of participants
  J <- max(BID$Block) # number of blocks
  Q <- nrow(BID) # number of questions
  npar <- 5*J+D*(D-1)/2
  #initial a parameters
  if(is.null(a)){
    if(!is.null(item.par)){
      a <- matrix(item.par$a,nrow = 3) #3 x J
    }else{
      a <- matrix(positive,nrow=3,ncol=J)
    }
  }
  #initial d parameters
  if(is.null(d)){
    if(!is.null(item.par)){
      d <- matrix(item.par$d,nrow = 3) #3 x J
    }else{
      d <- matrix(rnorm(3*J),nrow=3)
    }
    
  }
  #initial theta
  if(is.null(theta)){
    theta <- matrix(0,nrow = N,ncol = D)
  }
  #initial covariance matrix
  if(is.null(sigma)){
    sigm <- diag(D)
  }else{
    sigm <- sigma
  }
  
  # parallel computing settings
  
  nCPUcores = detectCores()
  if (nCPUcores < 3) {
    registerDoSEQ()
  }else{
    cl = makeCluster(nCPUcores - 1)
    registerDoParallel(cl)
  }
  plist <- list()
  
  ###########################
  #
  # burn-in phase
  #
  ###########################
  total.number.of.batch <- burn.in.size <- 0
  # the initial MxB iterations
  cat("\nBurn-in phase:")
  # the kernel function returns item parameter estimates, person parameter estimates and sigma estimates
  x <- kernel(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,cor.matrix = TRUE,positive=positive,fix.sigma=fix.sigma)
  a <- x$a
  d <- x$d
  sigm <- x$sigm
  theta <- x$theta
  for(m in 1:M){ 
    
    cat("\n  # of batch = ",m)
    x <- kernel(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,cor.matrix = (m<M/2),positive=positive,fix.sigma=fix.sigma)
    plist[[m]] <- x$pmatrix
    a <- x$a
    d <- x$d
    sigm <- x$sigm
    theta <- x$theta
    # print(sigm)
    if(m>=burnin.maxitr) break
  }
  #terminates burn-in using Geweke statistic z? --- see Page 3 of the manual 
  mcmc.par <- coda::mcmc(t(Reduce(cbind,plist)))
  z <- coda::geweke.diag(mcmc.par,frac1 = frac1,frac2 = frac2)$z
  print(z)
  # m <- M
  while(sum(z^2)>=npar*eps1&&m<burnin.maxitr){ # if the chain is not stable....
    if(m==M){
      cat("  sum z^2 / npar = ",sum(z^2)/npar, " criterion = ",eps1)
    }else{
      cat("\n  # of batch = ",m," sum z^2 / npar = ",sum(z^2)/npar, " criterion = ",eps1)
    }
    
    x <- kernel(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,positive=positive,fix.sigma=fix.sigma)
    a <- x$a
    d <- x$d
    sigm <- x$sigm
    #print(sigm)
    theta <- x$theta
    plist <- plist[-1]
    plist[[M]] <- x$pmatrix
    
    #terminates burn-in using Geweke statistic z? --- see Page 3 of the manual
    mcmc.par <- coda::mcmc(t(Reduce(cbind,plist)))
    z <- coda::geweke.diag(mcmc.par,frac1 = frac1,frac2 = frac2)$z
    cat("  sum z^2 / npar = ",sum(z^2)/npar)
    burn.in.size <- burn.in.size + B
    m <- m+1
  }
  
  cat("\n  # of batch = ",m," sum z^2 / npar = ",sum(z^2)/npar, " criterion = ",eps1)
  total.number.of.batch <- m
  cat("\nBurn in batch = ",1 + burn.in.size/B," burn-in iterations = ",B + burn.in.size)
  
  ###########################
  #
  # determining chain length
  #
  ###########################  
  n <- M
  d.hat <- batch.var(plist,n=n)
    cat("\nAfter burn-in phase:")
  while(max(d.hat*N)>=eps2&&n<maxitr){
    cat("\n  # of valid batch = ",n," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)
    n <- n+1
    
    x <- kernel(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,positive=positive,fix.sigma=fix.sigma)
    a <- x$a
    d <- x$d
    sigm <- x$sigm
    theta <- x$theta
    plist[[n]] <- x$pmatrix
    
    d.hat <- batch.var(plist,n=n)
  }
    cat("\n  # of valid batch = ",n," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)
    
  cat("\nLenth of final MC chain = ",n*B)
  total.number.of.batch <- total.number.of.batch + n - M
  est <- param.vec.2.ads(rowMeans(Reduce(cbind,plist)),J,D,sigm=sigm,fix.sigma=fix.sigma)
  
  stopCluster(cl)
  s2 <- Sys.time()
  return(list(a=est$a,d=est$d,sigm=est$sigm,total.number.of.batch=total.number.of.batch,
              final.chain=n*B,burn.in.size=burn.in.size,plist=plist,
              timeused = s2 - s1, start.time = s1, end.time = s2))
}




