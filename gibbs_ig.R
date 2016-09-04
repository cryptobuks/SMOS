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
sourceCpp("gibbscpp.cpp")

prior=list(a1=1,a2=10,b1=1,b2=10)
initial=list(psi1=1,psi2=1)
model_2harm<-dlmModTrig(s=276,q=2,dV=1/initial$psi1,dW=1/initial$psi2)
y=pix.long[[1]]
dat<-y$tau-mean(pix.long[[1]]$tau,na.rm=TRUE)
dat<-matrix(dat,1656)
l=mcmc_local_seas(2500,dat,initial,prior,model_2harm)

#remove burn in of 100
W=l[[3]][-c(1:100)]
W=l[[3]]
V=l[[2]][-c(1:100)]
V=l[[2]]

#Plot Results
burn=100
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
plot(m,type="l")

