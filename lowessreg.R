#Let's bring files back in and add variable names
folder <- "C:/Users/Colin/Documents/Caragea/TSDat/"      # path to folder that holds multiple .csv files
files.pix <- list.files(path=folder, pattern="*.txt") # create list of all .txt files in folder
files.pix2<-gsub("\\.txt.*","",files.pix)
# read in each .csv file in file_list and create a data frame with the same name as the .csv file
for (i in 1:length(files.pix)){
  assign(files.pix2[i], 
         read.table(paste(folder, files.pix[i], sep=''),header=TRUE)
  )}

pix.dat = lapply(files.pix2, get)

#make column names the same and make new dataset
for (i in seq_along(pix.dat)) {
  names(pix.dat[[i]]) <- c("index","2010","2011", "2012","2013","2014","2015")
}



#switch data from wide to stacked format.  We now have 30 t.s going over all years
short<-c(5,7,8,10,12,17,18,23,25)
long=c(1,2)

#Plot Regression with Loess Curve
setwd("C:/Users/Colin/Documents/Caragea/Plots")
pdf("plot_lowess_long.pdf")
for (i in long){
mod.3<-lm(tau~sin((2*pi*day)/276)+cos((2*pi*day)/276)+sin((4*pi*day)/276)+cos((4*pi*day)/276)+sin((6*pi*day)/276)+cos((6*pi*day)/276),data=pix.long[[i]])
plot(tau~day,data=pix.long[[i]],type="n")
lines((pix.long[[i]]$day[!is.na(pix.long[[i]]$tau)]),mod.3$fitted,col=2,lwd=2,lty=4)
title(files.pix2[i])

#Loess (.05 bw)

pix.lo<-loess(tau~day,data=pix.long[[i]],degree=2,span=.05)
allday <- seq(from=1, to=1656, by=1)
hat1 <- predict(pix.lo,allday)
lines(allday,hat1,col="blue",type="l",lwd=2)
#Add Shading by Years
abline(v=seq(0,1656,by=276))
color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
rect(276,-1,552,1,col=color)
rect(828,-1,1104,1,col=color)
rect(1380,-1,1656,1,col=color)
}
dev.off()

#short 275 days
pdf("plot_lowess_short.pdf")
for (i in short){
  mod.3<-lm(tau~sin((2*pi*day)/275)+cos((2*pi*day)/275)+sin((4*pi*day)/275)+cos((4*pi*day)/275)+sin((6*pi*day)/275)+cos((6*pi*day)/275),data=pix.long[[i]])
  plot(tau~day,data=pix.long[[i]],type="n")
  lines((pix.long[[i]]$day[!is.na(pix.long[[i]]$tau)]),mod.3$fitted,col=2,lwd=2)
  title(files.pix2[i])
  
  #Loess (.05 bw)
  
  pix.lo<-loess(tau~day,data=pix.long[[i]],degree=2,span=.05)
  allday <- seq(from=1, to=1656, by=1)
  hat1 <- predict(pix.lo,allday)
  lines(allday,hat1,lty=4,col="blue",type="l",lwd=2)
  #Add Shading by Years
  abline(v=seq(0,1650,by=275))
  color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
  rect(275,-1,550,1,col=color)
  rect(825,-1,1100,1,col=color)
  rect(1375,-1,1650,1,col=color)
}
dev.off()

loessGCV <- function (x) {
  ## Modified from code by Michael Friendly
  ## http://tolstoy.newcastle.edu.au/R/help/05/11/15899.html
  if (!(inherits(x,"loess"))) stop("Error: argument must be a loess object")
  ## extract values from loess object
  span <- x$pars$span
  n <- x$n
  traceL <- x$trace.hat
  sigma2 <- sum(resid(x)^2) / (n-1)
  gcv  <- n*sigma2 / (n-traceL)^2
  result <- list(span=span, gcv=gcv)
  result
}
bestLoess <- function(model, spans = c(.05, .95)) {
  f <- function(span) {
    mod <- update(model, span = span)
    loessGCV(mod)[["gcv"]]
  }
  result <- optimize(f, spans)
  result
}


