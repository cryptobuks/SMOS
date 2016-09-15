#Simulate Data with Different Harmonic Variances
require(dlm, quietly = TRUE)

#Try with 2 and 3 Harmonics


#3 harm initial values
inits3<-list(c(0,0,0,0,0,0),c(1,1,1,1,1,1),c(1,0,1,0,1,0),c(0,1,0,1,0,1),c(10,0,1,0,1,0),c(1,0,10,0,1,0),c(1,0,1,0,10,0),c(1,10,1,10,1,10))
#2 harm initial values
inits2<-list(c(0,0,0,0),c(1,1,1,1),c(1,0,1,0),c(0,1,0,1),c(10,0,1,0),c(1,0,10,0),c(1,10,1,0),c(1,0,1,10))

par(mfrow = c(4, 2))

#seed 353

sim<-function(q,start,seed){
set.seed(seed)
  for (i in 1:8){
    if (q==2){
        m=dlmModTrig(s=276,q=2,dV=1,dW=inits2[[i]],m0=c(0,0,0,0),C0=diag(c(1,1,1,1)))
  }
    else if (q==3){
        m=dlmModTrig(s=276,q=3,dV=1,dW=inits3[[i]],m0=c(0,0,0,0,0,0),C0=diag(c(1,1,1,1,1,1)))
    }
  o1 = dlmForecast(m,1656, sampleNew=1)
  dat = data.frame(y = as.numeric(o1$newObs[[1]]))
  
  
  plot(dat$y,type="l",ylab="y",main=paste("W=",start[i]))
  abline(v=seq(0,1656,by=276))
  color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
  rect(276,-500,552,500,col=color)
  rect(828,-500,1104,500,col=color)
  rect(1380,-500,1656,500,col=color)
  }
}

#joint = dlmGibbsDIG(o1$newObs[[1]], m1, a.y = 1, b.y = 1, a.theta = 1, b.theta = 1, n.sample = 1000, thin = 1)
