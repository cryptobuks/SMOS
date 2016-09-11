#Gibbs IG Prior

#IG Gamma Prior for V
a1<-1
b1<-10
#IG Gamma Prior for W
a2<-1
b2<-10

#Starting values
psi1<-1
psi2<-1

mcmc.ig<-function(i,a1,b1,a2,b2,psi1,psi2,days,reps,thin,burn){
  set.seed(20)
  y=pix.long[[i]]
  y.meanadj<-y$tau-mean(pix.long[[i]]$tau,na.rm=TRUE)

  mod_level<-dlmModTrig(s=days,q=2,dV=1/psi1,dW=1/psi2)
  every<-thin+1
  mc<-reps*every
  gibbsV<-numeric(reps)
  gibbsW<-numeric(reps)
  n<-length(y.meanadj)
  gibbsTheta<-array(0,dim=c(n+1,4,reps))

  #find non missing values of y_t
   y.not.na <- which(!is.na(y.meanadj))
   t.star<-length(y.not.na)
   sh1<- a1+t.star/2
   sh2<-a2+2*n
  it.save<-0
  for (it in 1:mc){
  #draw the states: FFBS
  filt<-dlmFilter(y.meanadj,mod_level)
  level<-dlmBSample(filt)
  
  #drop missing values of y_t
  state<-tcrossprod(level[-1, , drop = FALSE], mod_level$FF)
  state.nomis<-state[y.not.na]
  
  #draw obs precision; only sum over obs y's
  y.center<-y.meanadj[y.not.na] - state.nomis
  SSEy <- drop(crossprod(y.center))
  rate1<-b1+.5*SSEy
  psi1<-rgamma(1,shape=sh1,rate=rate1)
  
  #draw state precision
  theta.center<-level[-1,]-(level[-(n+1),])%*%t(mod_level$GG)
  SSE.theta<-(diag(crossprod(theta.center)))
  rate2<-b2 + 0.5*(sum(SSE.theta))
  psi2<- rgamma(1, shape = sh2, rate = rate2)
  
  #update 
  V(mod_level)<-1/psi1
  diag(W(mod_level))<-1/psi2
  #save samples
  if (!(it%%every)) {
    it.save <- it.save + 1
    gibbsV[it.save]<-1/psi1
    gibbsW[it.save]<-1/psi2
    gibbsTheta[,,it.save]<-level
  }
  }
  out=list(gibbsV[-(1:burn)],gibbsW[-(1:burn)],gibbsTheta[-1,,-c(1:burn)])
  return(out)
}


#Cpp Version
library(Rcpp)
library(dlm)
library(coda)
sourceCpp("gibbscpp.cpp")

#gamma (1,1) prior
prior=list(a1=1,a2=1,b1=1,b2=1)
#starting values
initial=list(psi1=1,psi2=1)
model_2harm<-dlmModTrig(s=276,q=2,dV=1/initial$psi1,dW=1/initial$psi2)
y=pix.long[[30]]
dat<-y$tau-mean(pix.long[[1]]$tau,na.rm=TRUE)
dat<-matrix(dat,1656)
l=mcmc_local_seas(2500,10,dat,initial,prior,model_2harm)

#Run 4 Chains with Different Starting Values
ChainV=matrix(NA,nrow=2500,ncol=4)
ChainW=matrix(NA,nrow=2500,ncol=4)
Theta<-list()

start<-list(c(psi1=1,psi2=1),c(psi1=.1,psi2=.1),c(psi1=10,psi2=10),c(psi1=10,psi2=20))
for (i in 1:4){
  initial<-start[[i]]
  l=mcmc_local_seas(2500,10,dat,initial,prior,model_2harm)
  ChainV[,i]<-l$sigma2_e
  ChainW[,i]<-l$sigma2_w
  Theta[i]<-l$theta
}

#Rhat for V 
mc1=as.mcmc(ChainV[-c(1:500),1])
mc2=as.mcmc(ChainV[-c(1:500),2])
mc3=as.mcmc(ChainV[-c(1:500),3])
mc4=as.mcmc(ChainV[-c(1:500),4])
mclist = mcmc.list(mc1,mc2,mc3,mc4)
gelman.diag(mclist, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)

#Rhat for W 
mc1w=as.mcmc(ChainW[-c(1:500),1])
mc2w=as.mcmc(ChainW[-c(1:500),2])
mc3w=as.mcmc(ChainW[-c(1:500),3])
mc4w=as.mcmc(ChainW[-c(1:500),4])
mclistw = mcmc.list(mc1w,mc2w,mc3w,mc4w)
gelman.diag(mclistw, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)



#remove burn in of 500
W=l[[3]][-c(1:500)]
W=l[[3]]
V=l[[2]][-c(1:500)]
V=l[[2]]

#Plot Results
burn=500
use<-2500-burn
from<-.05*use

par(mfrow=c(1,1))
plot(ergMean(V[-(1:burn)],from),type="l",xaxt="n",ylab="V",xlab="iter")
at<-pretty(c(0,use),n=3)
at<-at[at>=from]
axis(1,at=at-from,labels=format(at))

#Theta Values
l.sub=lapply(l[[1]], function(x) x[-1,])
l.sub2<-l.sub[-c(1:burn)]
l.sub3<-lapply(l.sub2,function(x) x[,1]+x[,3])
df <- data.frame(matrix(unlist(l.sub3), nrow=1656, byrow=F))
m=apply(df, 1, mean)
l=apply(df, 1, quantile, probs = .025,  na.rm = TRUE)
u=apply(df, 1, quantile, probs = .975,  na.rm = TRUE)
plot(m,type="l",ylab="T",xlab="Year",xaxt="n",ylim=c(-.35,.35))
lines(l,lty=3);lines(u,lty=3)
points(dat, col = "seagreen",ylab=expression(tau),main="Mean Adjusted Pix 194406",xaxt="n")
abline(v=seq(0,1656,by=276))
color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
rect(276,-1,552,1,col=color)
rect(828,-1,1104,1,col=color)
rect(1380,-1,1656,1,col=color)
axis(1, at=c(0,276,552,828,1104,1380,1656), labels=c("2010", "2011", "2012","2013","2014","2015","2016"))


