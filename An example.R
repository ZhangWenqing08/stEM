######### An example
source('stEM.R')
source('utils.R')
required.packages <- c("armspp", "doParallel","foreach","doRNG","coda","mvnfast","lvmcomp")
lapply(required.packages,require, character.only = TRUE)
D <- 6
nitem.per.dim <- 10
nblock <- D * nitem.per.dim / 3
set.seed(123456)

# item parameters
item.par <- data.frame(a=seq_len(D*nitem.per.dim))
item.par <- within(item.par,{
  a <- runif(D*nitem.per.dim,0.7,3)
  a[1:5] <- -1*a[1:5]
  a[8] <- -1*a[8]
  # a <- rlnorm(D*nitem.per.dim,sdlog = 0.5)
  b <- rnorm(D*nitem.per.dim)
  d <- a*b
})

#BID
BID <- data.frame(Block=rep(1:nblock,each=3),
                  Item=rep(1:3,nblock),
                  Dim=c(combn(D,3)[,sample(choose(D,3),nblock,replace = TRUE)]))

N <- 1000
v <- matrix(0.5,D,D)
diag(v) <- 1
eigen(v)$values
theta <- mvnfast::rmvn(N,seq(-1,1,length.out = D),sigma = v)

Y <- data.sim(item.par,theta,BID)

item.par$d <- c(t(aggregate(item.par$d,by=list(BID$Block),function(x)x-mean(x))[,-1]))



x <- stEM(Y,BID,maxitr = 100,positive = sign(item.par$a),sigma = v,fix.sigma = FALSE)
plot(c(x$a),item.par$a)
plot(c(x$d),item.par$d)
cor(c(x$a),item.par$a)
cor(c(x$d),item.par$d)
x$sigm
cat("\nEstimating person abilities...")

theta.est <- matrix(0,N,D)
for (i in 1:nrow(Y)){
  
  theta.est[i,] <- optim(par=theta[i,],fn = logLi,method = "L-BFGS-B",lower = rep(-5,D), upper = rep(5,D),
                         a=x$a, d=x$d,yi=Y[i,], BID=BID,prior=TRUE,sigma=x$sigm,control=list(fnscale=-1))$par
}
diag(cor(theta,theta.est))
cor(theta.est)
cor(theta)