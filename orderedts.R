library(plyr)
library(dplyr)
library(zoo)
setwd("C:/Users/Colin/Documents/Caragea/TSDat")

#Load Data (Note, total of 30 pixels, 180 files)
folder <- "C:/Users/Colin/Documents/Caragea/Data/"      # path to folder that holds multiple .csv files
file_list <- list.files(path=folder, pattern="*.csv") # create list of all .csv files in folder
file_list2<-gsub("\\.csv.*","",file_list)
# read in each .csv file in file_list and create a data frame with the same name as the .csv file
for (i in 1:length(file_list)){
  assign(file_list2[i], 
         read.csv(paste(folder, file_list[i], sep=''),header=TRUE)
  )}

#Get a list of order of pix
seqpix<-seq(1,180,by=6)
namespix=file_list2[seqpix]
spl=strsplit(namespix, "_",fixed=TRUE)
namesfin=unlist(lapply(spl,FUN=function(x){paste(x[1],x[2],sep="_")}))

#Make the list access the data
list_df = lapply(file_list2, get)


#make a list of all 30 pixels (note: can make as .as.matrix too)
series.days<-list()
for (i in 0:29){
  lower=1+i*6
  upper=6+i*6
  series.days[[i+1]]<-lapply(list_df[lower:upper], as.data.frame)
}

#add integer dates for all 30 pixels
out <- vector("list", length(series.days))
for (j in 1:30){
  out[[j]]=lapply(series.days[[j]],function(x) transform(x, doyint = as.integer(x[,1])))
}


#Drop NA on dates
days.nomiss <- vector("list", length(out))
for (j in 1:30){
  days.nomiss[[j]]=lapply(out[[j]],na.omit)
}


#Now return the list back to 180 files
avgdup<-unlist(days.nomiss, recursive = FALSE)

#Drop Obs with Tau=0
taupositive <- vector("list", length(avgdup))
for(i in 1:180){
  taupositive[[i]]=avgdup[[i]][avgdup[[i]]$tau>0,]
}



#how many duplicates per year?  Assume 276 days per year even though some are 275
dup.total<-sapply(taupositive, function(x) sum(duplicated(x[,7]))/276) 
hist(dup.total,main="Percent duplicate per year for 180 series")

#Are there more than 3 duplicate?
tridup<-NULL
for (i in 1:180){
  x=table(taupositive[[i]]$doyint)
  tridup[i]=max(x)
}
morethreedup=which(tridup %in% c(3))

#Look at Difference between times
sat.time <- vector("list", length(taupositive))
for (j in 1:180){
test=taupositive[[j]]
index=test[duplicated(test$doyint),7]
test3=test[test$doyint %in% index,]
l=diff(test3$doy)
sat.time[[j]]=l[l<.95 & l!=0]
}

plot(sat.time[[1]],type="l",xlim=c(0,50),ylim=c(.45,1),main="Time gap between days with double measurements",ylab="time (fractional day)")
for (i in 2:180){
lines(sat.time[[i]],type="l",xlim=c(0,50),ylim=c(.4,1))
}
abline(h=.46,lty=3)
abline(h=.58,lty=3)

#Average over duplicate days
no.dup<-lapply(taupositive, function(x) aggregate(tau~doyint,data=x,FUN = mean))

#What's maximum number of days with data?
maxnull<-NULL
for (i in 1:180){
  maxnull[i]=max(no.dup[[i]]$doyint)
}

#Turn all 180 files to zoo objects
pix.zoo <- lapply(no.dup, read.zoo, index=1)
print(pix.zoo[[1]],style = "vertical")

#Number of non-missing and no duplicate obs per pixel per year
ct=sapply(pix.zoo,length)
ct

#split up lists to sublists by pixel, merge using zoo, export file
for (i in 0:29){
  lower=1+i*6
  upper=6+i*6
  l1=pix.zoo[lower:upper]
  unionl1 <- Reduce(merge,l1)
  g=zoo(,seq(start(unionl1),end(unionl1)))
  z <- merge(unionl1,g)
  write.zoo(z,file = paste(namesfin[i+1], ".txt"))
}

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

t.length=sapply(pix.long, function (x) length(x$tau))
table(t.length)

#how many NA over 2010-2015 for each pixel
par(mfrow=c(1,1))
miss<-sapply(pix.long, function (x) sum(is.na(x$tau))/length(x$tau))
hist(miss,main="Percent Missing for Each Pixel")
table(miss)


#Add Day Index Column
for (i in 1:30){
pix.long[[i]]$day<-seq.int(nrow(pix.long[[i]]))
}

#Plots for Pix 23,26,28,30; add cool shading too
par(mfrow=c(4,1))
par(mar=c(3,4,1,1))
par(oma=c(0,0,2,0))
plot(pix.long[[23]]$tau,type="l",col="blue",ylab=expression(tau))
plot(pix.long[[26]]$tau,type="l",col="blue",ylab=expression(tau))
plot(pix.long[[28]]$tau,type="l",col="blue",ylab=expression(tau))
plot(pix.long[[30]]$tau,type="l",col="blue",ylab=expression(tau))
title("Pixels 200560, 201589, 202095, 203632", outer=T)

abline(v=seq(0,1656,by=276))
color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
rect(276,-1,552,1,col=color)
rect(828,-1,1104,1,col=color)
rect(1380,-1,1656,1,col=color)

#Plots for Pix 29,27,25,24, not pix 25 has 275 days
par(mfrow=c(4,1))
par(mar=c(3,4,1,1))
plot(pix.long[[29]]$tau,type="l",col="blue",ylab=expression(tau))
plot(pix.long[[27]]$tau,type="l",col="blue",ylab=expression(tau))
plot(pix.long[[25]]$tau,type="l",col="blue",ylab=expression(tau))
plot(pix.long[[24]]$tau,type="l",col="blue",ylab=expression(tau))
title("Pixels 202112, 201598, 201083, 201081", outer=TRUE)

abline(v=seq(0,1656,by=276))
color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)
rect(276,-1,552,1,col=color)
rect(828,-1,1104,1,col=color)
rect(1380,-1,1656,1,col=color)

