library(dlm)

#5,7,8,10,12,17,18,23,25 have 275 days per year

mod2h<-function(pix, days){
  i<-pix
  j<-days
  y=pix.long[[i]]
  y.meanadj<-y$tau-mean(pix.long[[i]]$tau,na.rm=TRUE)
  
  buildFun2<-function(parm){
  dlmModTrig(s=j,q=2,dV=exp(parm[1]),dW=c(exp(parm[2])))
}

  fit<-dlmMLE(y.meanadj,rep(0,2),build=buildFun2,hessian=TRUE)
  dlm2 <- buildFun2(fit$par)
  dlm1Smooth2 <- dlmSmooth(y.meanadj, dlm2)
  smooth.sum<-dropFirst(dlm1Smooth2$s[,1])+dropFirst(dlm1Smooth2$s[,3])

  #Make Plot
  plot(y.meanadj, col = "black",ylab="T",xlab="Year",main=paste("Mean Adjusted Pix",i,"Smoothed"),xaxt="n")
  abline(v=seq(0,j*6,by=j))
  color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
  rect(j,-1,2*j,1,col=color)
  rect(j*3,-1,j*4,1,col=color)
  rect(j*5,-1,j*6,1,col=color)
  axis(1, at=c(0,j,j*2,j*3,j*4,j*5,j*6), labels=c("2010", "2011", "2012","2013","2014","2015","2016"))
  lines(smooth.sum,col = "red",lwd=2)
  return(cbind(dlm2$V,dlm2$W[1]))
}