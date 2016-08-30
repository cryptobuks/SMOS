#First Load Data from Fouriereg.R file
library(dlm)
#Look at 5 Pixels (from Hornbuckle): 1, 10* (275 days),20,27,30
horn=c(1,20,27,30)

#276 day pix
pixall=seq(1,30,1)
xday<-c(5,7,8,10,12,17,18,23,25)
dayplus<-pixall[-xday]

#build model with 3 harmonics; use fixed variance for all W;
buildFun<-function(parm){
  dlmModTrig(s=276,q=3,dV=exp(parm[1]),dW=(exp(parm[2])))
}

#build model with 2 harmonics 
buildFun2<-function(parm){
  dlmModTrig(s=276,q=2,dV=exp(parm[1]),dW=c(exp(parm[2])))
}

#build model with 1 harmonics
buildFun3<-function(parm){
  dlmModTrig(s=276,q=1,dV=exp(parm[1]),dW=exp(parm[2]))
}

#build model with 2 harmonics and random walk mean
dlm_nifa <- dlmModPoly(1) + dlmModTrig(s=276,q=2)
buildFun4 <- function(x) {
  diag(W(dlm_nifa))[2:5] <- exp(x[1])
  diag(W(dlm_nifa))[1]<-exp(x[2])
  V(dlm_nifa) <- exp(x[3])
  return(dlm_nifa)
}

fit4 <- dlmMLE(y.meanadj, parm = rep(0, 3), build = buildFun4,hessian=TRUE)
dlm_out <- buildFun4(fit4$par)
drop(V(dlm_out))
diag(W(dlm_out))[1:2]

#Smooth
nifSmooth <- dlmSmooth(y.meanadj, mod = dlm_out)
plot(y.meanadj, col = "seagreen",ylab=expression(tau),main="Mean Adjusted Pix 194406: One Step Ahead Prediction (2 Harmonics & Random Walk)",xaxt="n")
abline(v=seq(0,1656,by=276))
color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
rect(276,-1,552,1,col=color)
rect(828,-1,1104,1,col=color)
rect(1380,-1,1656,1,col=color)
axis(1, at=c(0,276,552,828,1104,1380,1656), labels=c("2010", "2011", "2012","2013","2014","2015","2016"))
lines(dropFirst(nifSmooth$s[,c(3)]), type = 'o', pch = 20, col = "red")

#Filter
dlmFilt4 <- dlmFilter(y.meanadj, dlm_out)
lines(dropFirst(dlmFilt4$f), pch = 20, col = "orange")


V=matrix(NA,nrow=30,ncol=4)
W=matrix(NA,nrow=30,ncol=5)
MSE=matrix(NA,nrow=30,ncol=4)

for (i in dayplus){
  y=pix.long[[i]]
  y.meanadj<-y$tau-mean(pix.long[[i]]$tau,na.rm=TRUE)

fit<-dlmMLE(y.meanadj,rep(1,2),build=buildFun,hessian=TRUE)
dlm1 <- buildFun(fit$par)
V[i,1]=V(dlm1)
W[i,1]=W(dlm1)[1]

fit2<-dlmMLE(y.meanadj,rep(0,2),build=buildFun2,hessian=TRUE)
dlm2 <- buildFun2(fit2$par)
V[i,2]=V(dlm2)
W[i,2]=W(dlm2)[1]

fit3<-dlmMLE(y.meanadj,rep(0,2),build=buildFun3,hessian=TRUE)
dlm3 <- buildFun3(fit3$par)
V[i,3]=V(dlm3)
W[i,3]=W(dlm3)[1]

fit4 <- dlmMLE(y.meanadj, parm = rep(0, 3), build = buildFun4,hessian=TRUE)
dlm4 <- buildFun4(fit4$par)
V[i,4]=(V(dlm_out))
W[i,4]=diag(W(dlm4))[1]
W[i,5]=diag(W(dlm4))[2]

#Filtering, Kalman Filter

dlm1Filt <- dlmFilter(y.meanadj, dlm1)
dlm1Filt2 <- dlmFilter(y.meanadj, dlm2)
dlm1Filt3 <- dlmFilter(y.meanadj, dlm3)
dlm1Filt4 <- dlmFilter(y.meanadj, dlm4)

#Compare Models
resid1=(residuals(dlm1Filt,type="raw",sd=FALSE))
resid2=(residuals(dlm1Filt2,type="raw",sd=FALSE))
resid3=(residuals(dlm1Filt3,type="raw",sd=FALSE))
resid4=(residuals(dlm1Filt4,type="raw",sd=FALSE))
MSE[i,1]=mean(abs(resid1),na.rm = TRUE)
MSE[i,2]=mean(abs(resid2),na.rm = TRUE)
MSE[i,3]=mean(abs(resid3),na.rm = TRUE)
MSE[i,4]=mean(abs(resid4),na.rm = TRUE)
}

par(mfrow=c(1,4))
hist(MSE[,1],main="3 Harm Model",xlab="MAPE")
hist(MSE[,2],main="2 Harm Model",xlab="MAPE")
hist(MSE[,3],main="1 Harm Model",xlab="MAPE")
hist(MSE[,4],main="2 Harm and RW",xlab="MAPE")



par(mfrow=c(1,4))
hist(V[,1],main="3 Harm",xlab="V")
hist(V[,2],main="2 Harm",xlab="V")
hist(V[,3],main="1 Harm",xlab="V")
hist(V[,4],main="2 Harm and RW",xlab="V")

par(mfrow=c(1,5))
hist(W[,1],main="3 Harm",xlab="W")
hist(W[,2],main="2 Harm",xlab="W")
hist(W[,3],main="1 Harm",xlab="W")
hist(W[,4],main="2 Harm and RW Var",xlab="W")
hist(W[,5],main="2 Harm Var and RW",xlab="W")

#Filtering, Kalman Filter

dlm1Filt <- dlmFilter(y.meanadj, dlm1)
dlm1Filt2 <- dlmFilter(y.meanadj, dlm2)
dlm1Filt3 <- dlmFilter(y.meanadj, dlm3)

#Compare Models
resid1=(residuals(dlm1Filt,type="raw",sd=FALSE))
resid2=(residuals(dlm1Filt2,type="raw",sd=FALSE))
resid3=(residuals(dlm1Filt3,type="raw",sd=FALSE))
mean(abs(resid1),na.rm = TRUE)
mean(abs(resid2),na.rm = TRUE)
mean(abs(resid3),na.rm = TRUE)


plot(y.meanadj, col = "seagreen",ylab=expression(tau),main="Mean Adjusted Pix 194406: One Step Ahead Prediction (1,2,3 Harmonics)",xaxt="n")
abline(v=seq(0,1656,by=276))
color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
rect(276,-1,552,1,col=color)
rect(828,-1,1104,1,col=color)
rect(1380,-1,1656,1,col=color)
axis(1, at=c(0,276,552,828,1104,1380,1656), labels=c("2010", "2011", "2012","2013","2014","2015","2016"))

#f is one step ahead predictive distribution of y_t
lines(dropFirst(dlm1Filt3$f), pch = 20, col = "orange")
lines(dropFirst(dlm1Filt2$f), pch = 20, col = "red")
lines(dropFirst(dlm1Filt$f), pch = 20, col = "blue")

#m is the filtered distributions of state vector Theta_t given y1:t
#a is the one step ahead predictive distribution of state vector

lines(dropFirst(dlm1Filt$m[,1]), type = 'o', pch = 20, col = "red")
lines(dropFirst(dlm1Filt$m[,2]), type = 'o', pch = 20, col = "blue")
lines(dropFirst(dlm1Filt$m[,3]), type = 'o', pch = 20, col = "pink")
lines(dropFirst(dlm1Filt$m[,4]), type = 'o', pch = 20, col = "yellow")
lines(dropFirst(dlm1Filt$m[,5]), type = 'o', pch = 20, col = "orange")
lines(dropFirst(dlm1Filt$m[,6]), type = 'o', pch = 20, col = "green")


#Dlm Smooth (Only Interested in the Odd States)
dlm1Smooth <- dlmSmooth(y$meanadj, dlm1)
dlm1Smooth2 <- dlmSmooth(y$meanadj, dlm2)

plot(y$meanadj, col = "seagreen",ylab=expression(tau),main="Mean Adjusted Pix 194406: Smoothed State Space (3 Harmonics)")
lines(dropFirst(dlm1Smooth$s[,1]), type = 'o', pch = 20, col = "red")
#lines(dropFirst(dlm1Smooth$s[,2]), type = 'o', pch = 20, col = "blue")
lines(dropFirst(dlm1Smooth$s[,3]), type = 'o', pch = 20, col = "pink")
#lines(dropFirst(dlm1Smooth$s[,4]), type = 'o', pch = 20, col = "yellow")
lines(dropFirst(dlm1Smooth$s[,5]), type = 'o', pch = 20, col = "orange")
#lines(dropFirst(dlm1Smooth$s[,6]), type = 'o', pch = 20, col = "green")

#residual checks
qqnorm(residuals(dlm1Filt,sd=FALSE))
qqline(residuals(dlm1Filt,sd=FALSE))
tsdiag(dlm1Filt)

#Bayes Approach using pg. 166; Dropping Missing y's in sigma_v calculation

#Select Pixel and Center It
i=1
y=pix.long[[i]]
y.meanadj<-y$tau-mean(pix.long[[i]]$tau,na.rm=TRUE)
y.meanadjfix<-y.meanadj

#IG Gamma Prior for V
a1<-1
b1<-10
#IG Gamma Prior for W
a2<-1
b2<-10

#Starting values
psi1<-1
psi2<-1

#Assume Pix has 276 days per cycle
mod_level<-dlmModTrig(s=276,q=2,dV=1/psi1,dW=1/psi2)
n.sample<-100
thin=10
burn<-10
every<-thin+1
mc<-n.sample*every
gibbsV<-numeric(n.sample)
gibbsW<-numeric(n.sample)
gibbsTheta<-array(0,dim=c(n+1,4,n.sample))

#find non missing values of y_t
y.not.na <- which(!is.na(y.meanadj))
n<-length(y.meanadj)
t.star<-length(y.not.na)
sh1<- a1+t.star/2
sh2<-a2+2*n


set.seed(10)
it.save<-0
for (it in 1:mc){
  #draw the states: FFBS
  filt<-dlmFilter(y.meanadj,mod_level)
  level<-dlmBSample(filt)
  
  #drop missing values of y_t
  state<-tcrossprod(level[-1, , drop = FALSE], mod_level$FF)
  state.nomis<-state[y.not.na]
  #impute=rnorm(length(y.na),state[y.na],sd=sqrt(1/psi1))
  #y.meanadj[y.na]=impute
  
  #draw obs precision; only sum over obs y's
  y.center<-y.meanadj[y.not.na] - state.nomis
  SSEy <- drop(crossprod(y.center))
  #rate<- b1+ crossprod(y.meanadj-drop((mod_level$FF%*%t(level[-1,]))))*.5
  rate1<-b1+.5*SSEy
  psi1<-rgamma(1,shape=sh1,rate=rate1)
  
  #draw state precision
  theta.center<-level[-1,]-(level[-(n+1),])%*%t(mod_level$GG)
  SSE.theta<-(diag(crossprod(theta.center)))
  #tt.theta.center<-t(level[-1,]-(level[-(n+1),])%*%t(mod_level$GG))
  #theta.center<-((level[-1,]-(level[-(n+1),])%*%t(mod_level$GG)))%*%(t(level[-1,]-(level[-(n+1),])%*%t(mod_level$GG)))
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

#Plot Results
use<-n.sample-burn
from<-.05*use

par(mfrow=c(1,1))
plot(ergMean(gibbsW[-(1:burn)],from),type="l",xaxt="n",ylab="W",xlab="iter")
at<-pretty(c(0,use),n=3)
at<-at[at>=from]
axis(1,at=at-from,labels=format(at))

#Look at ACF
acf(sqrt(gibbsV[-(1:burn)]),main="ACF for V")
acf(sqrt(gibbsW[-(1:burn)]),main="ACF for W")

mcmcMean(gibbsV[-(1:burn)])
mcmcMean(gibbsW[-(1:burn)])

#Plot Means
gibbsTheta2<-gibbsTheta[-1,,-c(1:10)]
gibbsTheta3<-gibbsTheta2[,1,]+gibbsTheta2[,3,]
thetaMean<-ts(apply(gibbsTheta3,1,mean),start=1, end=1656,frequency = 1)
LprobLim<-ts(apply(gibbsTheta3,1,quantile,probs=.025),start=1,end=1656,frequency=1)
UprobLim<-ts(apply(gibbsTheta3,1,quantile,probs=.975),start=1,end=1656,frequency=1)
plot(thetaMean,xlab="",ylab=expression(tau),ylim=c(-.75,.75),xaxt="n")
lines(LprobLim,lty=3);lines(UprobLim,lty=3)


#plot results
points(y.meanadj, col = "seagreen",ylab=expression(tau),main="Mean Adjusted Pix 194406",xaxt="n")
abline(v=seq(0,1656,by=276))
color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
rect(276,-1,552,1,col=color)
rect(828,-1,1104,1,col=color)
rect(1380,-1,1656,1,col=color)
axis(1, at=c(0,276,552,828,1104,1380,1656), labels=c("2010", "2011", "2012","2013","2014","2015","2016"))


#Try dlmGibbsDIG
gibbsOut<-dlmGibbsDIG(y.meanadj,mod=dlmModTrig(s=276,q=2),a.y = 1,b.y=10,a.theta=1,b.theta = 10,n.sample=1000)
