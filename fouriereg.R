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
pix.long <- vector("list", length(pix.dat))
for (i in 1:30){
  pix.long[[i]] <- reshape(pix.dat[[i]], 
                           varying = c("2010", "2011", "2012", "2013", "2014","2015"), 
                           v.names = "tau",
                           timevar = "order", 
                           idvar = "index",
                           direction = "long")
}

#Add Day Index Column
for (i in 1:30){
  pix.long[[i]]$day<-seq.int(nrow(pix.long[[i]]))
}

#Fourier Regression
#Note: files.pix2 gives you pixel numbers for the list
#Fit Regression
mod.2<-lm(tau~sin((2*pi*day)/276)+cos((2*pi*day)/276)+sin((4*pi*day)/276)+cos((4*pi*day)/276),data=pix.long[[1]])
summary(mod.2)

setwd("C:/Users/Colin/Documents/Caragea/Plots")
pdf("plots.pdf")
#275days
short=c(5,7,8,10,12,17,18,23,25)
long=c(1,2,3,4,6,9,11,13,14,15,16,19,20,21,22,24,26,27,28,29,30)
for (i in long){
mod.3<-lm(tau~sin((2*pi*day)/276)+cos((2*pi*day)/276)+sin((4*pi*day)/276)+cos((4*pi*day)/276)+sin((6*pi*day)/276)+cos((6*pi*day)/276),data=pix.long[[i]])
plot(tau~day,data=pix.long[[i]],type="l")
lines((pix.long[[i]]$day[!is.na(pix.long[[i]]$tau)]),mod.3$fitted,col=2,lwd=3)
title(files.pix2[i])
#Add Shading by Years
abline(v=seq(0,1656,by=276))
color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
rect(276,-1,552,1,col=color)
rect(828,-1,1104,1,col=color)
rect(1380,-1,1656,1,col=color)
}
dev.off()

#Plots for Cells with 275 days per year
pdf("plots_short.pdf")
#275days
short=c(5,7,8,10,12,17,18,23,25)
for (i in short){
  mod.3<-lm(tau~sin((2*pi*day)/275)+cos((2*pi*day)/275)+sin((4*pi*day)/275)+cos((4*pi*day)/275)+sin((6*pi*day)/275)+cos((6*pi*day)/276),data=pix.long[[i]])
  plot(tau~day,data=pix.long[[i]],type="l")
  lines((pix.long[[i]]$day[!is.na(pix.long[[i]]$tau)]),mod.3$fitted,col=2,lwd=3)
  title(files.pix2[i])
  #Add Shading by Years
  abline(v=seq(0,1650,by=275))
  color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
  rect(275,-1,550,1,col=color)
  rect(825,-1,1100,1,col=color)
  rect(1375,-1,1650,1,col=color)
}
dev.off()

