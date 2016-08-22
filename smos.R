library(plyr)
library(dplyr)
library(stats)
library(zoo)
#Load Data (Note, total of 30 pixels)
folder <- "C:/Users/Colin/Documents/Caragea/"      # path to folder that holds multiple .csv files
file_list <- list.files(path=folder, pattern="*.csv") # create list of all .csv files in folder
file_list2<-gsub("\\.csv.*","",file_list)
# read in each .csv file in file_list and create a data frame with the same name as the .csv file
for (i in 1:length(file_list)){
  assign(file_list2[i], 
         read.csv(paste(folder, file_list[i], sep=''),header=TRUE)
  )}

#Make the list access the data
list_df = lapply(file_list2, get)

#Recode All Missing Values to NA
for (i in 1:length(list_df)){
   list_df[[i]][[6]][list_df[[i]][[6]]==-999]<-NA
   list_df[[i]][[1]][list_df[[i]][[1]]<0]<-NA
}

#make column names the same and make new dataset
for (i in seq_along(list_df)) {
  names(list_df[[i]]) <- c("doy","dgg","lat", "lon","sm","tau")
}

setwd("C:/Users/Colin/Documents/Caragea/Data/")
n <- 1:length(file_list2)
lapply(n, function(ni) {
  write.table(file = paste(file_list2[ni], ".csv", sep = ""), 
              list_df[ni], sep = ',', row.names = F)
}
)

#How to load data NOW
library(plyr)
#Load Data (Note, total of 30 pixels, 180 files)
folder <- "C:/Users/Colin/Documents/Caragea/Data/"      # path to folder that holds multiple .csv files
file_list <- list.files(path=folder, pattern="*.csv") # create list of all .csv files in folder
file_list2<-gsub("\\.csv.*","",file_list)
# read in each .csv file in file_list and create a data frame with the same name as the .csv file
for (i in 1:length(file_list)){
  assign(file_list2[i], 
         read.csv(paste(folder, file_list[i], sep=''),header=TRUE)
  )}

#Make the list access the data
list_df = lapply(file_list2, get)



#Remove a pixel from list and turn it into a dataframe take 194406
xm1<-do.call("rbind", lapply(list_df[1:6],as.matrix))

#can also make this a data.frame, not sure what's better
xm2<-do.call("rbind", lapply(list_df[7:12], as.data.frame))

#turn this into loop
series<-list()
for (i in 0:29){
lower=1+i*6
upper=6+i*6
series[[i+1]]<-do.call("rbind", lapply(list_df[lower:upper], as.data.frame))
}

#Number of rows in each matrix with missing excluded
L <- lapply(series, na.exclude)
dim.series<-sapply(L,nrow)
dim.series
hist(dim.series,xlab="Numer of obs for 2010-2015",main="All 30 pixels Non-Missing Obs")

plot(xm2$tau,type = "l")
#linear interpolation between NA values
plot(na.approx(xm2$tau),type="l")

#Count Days Loop
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

test<-unlist(days.nomiss, recursive = FALSE)

no.dup <- vector("list", length(out))
no.dup2<-lapply(test, function(x) aggregate(tau~doyint,data=x,FUN = mean))


#Turn all 180 files to zoo objects
Lz <- lapply(no.dup2, read.zoo, index=1)
print(Lz[[1]],style = "vertical")

#Check order
ct=sapply(Lz,length)


#split up lists to sublists
list12=Lz[1:6]
df.list <- do.call(merge,list12,all=TRUE)
print(h,style="vertical")
m44=merge(Lz[[1]],Lz[[2]],all=TRUE)
#Merge List and add rows for missing days
jj=do.call("merge",list12)
z <- merge(zoo(,start(jj):end(jj)), jj)

#export
write.zoo(z,file="bday.txt")

x <- merge(jj, zoo(,seq(start(jj),end(jj),by="%d")), all=TRUE)

write.zoo(jj,file="testbday",index.name = "Index")
jj33=as.matrix(jj)
print(jj,style = "vertical")

match.all.hi=full_join(hi,hi2, by="doyint")

#Make a list with just 199032 and add integer days 
pixfull<-list_df[175:180]
pixfull2<-lapply(pixfull,function(x) transform(x, doyint = as.integer(x[,1])))
#Make Integer dates
lapply(series.days[[2]],function(x) transform(x, doyint = as.integer(x[,1])))

#Just Pull Day Integer Column
days.all <- vector("list", length(series.days))
for (j in 1:30){
  days.all[[j]]=lapply(out[[j]], "[", c(7))
}


#Drop NA Values
days.nomiss <- vector("list", length(series.days))
for (j in 1:30){
  days.nomiss[[j]]=lapply(days.all[[j]],na.omit)
}

#total days per pixel
test<-do.call("c", days.nomiss,recurs=F)
test<-unlist(days.nomiss, recursive = FALSE)


#join to find matching values
match.all <- vector("list", length(series.days))
for (i in 1:30){
match.all[[i]]=join_all(days.nomiss[[i]], by="doyint",type='inner',match="all")
}
#Get Unique Matches, i.e. drop duplicates
c <- vector(mode="numeric", length=30)
for (i in 1:30){
x<-length(unique(match.all[[i]]$doyint))
c[i]=x
}
table(c)
hist(c,xlab="number of matching days over all 5 years",main="By Pixel number of matching Days over all years")


#Pull days for one pixel
days<-lapply(pixfull2, "[", c(7))
#Drop NA
days2<-lapply(days, na.omit)
#just get unique values for each year, can also do this at the end
days3<-lapply(days2,unique)
obs.days<-sapply(days2, function(x) length(unlist(x[]))) 
days.list <- as.list(days)
sum(obs.days)
r<-merge(days2[1], days2[2], by="doyint")
r2<-merge(r[1], days2[3], by="doyint")
r3<-merge(r2[1], days2[4], by="doyint")
r4<-merge(r3[1], days2[5], by="doyint")
r5<-merge(r4[1], days2[6], by="doyint")

m=join_all(days2, by="doyint",type='inner',match="all")
m2=join_all(days3, by="doyint",type='inner',match="all")
length(unique(m$doyint))



#Count not NA, work with "out" list
na.all <- vector("list", length(out))
for (i in 1:30){
na.all[[i]]<-lapply(out[[i]],function(x) sum(!is.na(x[,7]))/length(x[,7]))
}
dsn.notmissing<-unlist(na.all)
hist(dsn.notmissing,xlab="Percent Non Missing",main="All Years and All Pixels (180 obs)")

#Length of each year (non missing)
na.countall <- vector("list", length(out))
for (i in 1:30){
  na.countall[[i]]<-lapply(out[[i]],function(x) sum(!is.na(x[,7])))
}
count.notmissing<-unlist(na.countall)
hist(count.notmissing,xlab="Number of non-missing obs by year",main="All Years and All Pixels (180 obs)")

#Count Duplicates within pixel
dup.all <- vector("list", length(out))
for (i in 1:30){
  dup.all[[i]]<-lapply(out[[i]],function(x) sum(duplicated(x[,7]))/length(x[,7]))
}
dsn.dup<-unlist(dup.all)
hist(dsn.dup,xlab="percent duplicates",main="All Years and All Pixels (180 obs)")


#Scaled Periodogram 
test<-na.omit(SM_199032_full)
tau<-test$tau-mean(test$tau)
I=abs(fft(tau)/sqrt(1678))^2
P = (4/1678)*I[1:840] # Only need the first (n/2)+1 values of the FFT result.
f = (0:839)/1678 # this creates harmonic frequencies from 0 to .5 in steps of 1/128.
plot(f, P, type="l") # This plots the periodogram;

#Fit Regression
mod1<-lm(tau~sin((2*pi*t)/279.67)+cos((2*pi*t)/279.67)+sin((4*pi*t)/279.67)+cos((4*pi*t)/279.67)+sin((6*pi*t)/279.67)+cos((6*pi*t)/279.67),data=r)
mod2<-lm(tau~sin((2*pi*t)/279.67)+cos((2*pi*t)/279.67)+sin((4*pi*t)/279.67)+cos((4*pi*t)/279.67),data=r)
summary(mod1)
plot(tau~t,data=r)
lines(r$t,mod1$fitted,col=2)
lines(r$t,mod2$fitted.values,col=3)



#Pickel 199032
par(mfrow=c(3,2))
plot(SM_199032._2015$doy,SM_199032._2015$tau,xaxt='n',xlab="Days",main="2015",xlim=c(59,335),ylab=expression(tau),col="red")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_199032_2014$doy,SM_199032_2014$tau,xaxt="n",xlab="Days",main="2014",xlim=c(59,335),ylab=expression(tau),col="black")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

plot(SM_199032_2013$doy,SM_199032_2013$tau,xaxt="n",xlab="Days",main="2013",xlim=c(59,335),ylab=expression(tau),col="blue")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_199032_2012$doy,SM_199032_2012$tau,xaxt="n",xlab="Days",main="2012",xlim=c(59,335),ylab=expression(tau),col="green")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

plot(SM_199032_2011$doy,SM_199032_2011$tau,xaxt="n",xlab="Days",main="2011",xlim=c(59,335),ylab=expression(tau),col="orange")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_199032_2010$doy,SM_199032_2010$tau,xaxt="n",xlab="Days",main="2010",xlim=c(59,335),ylab=expression(tau),col="brown")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

#Pixle 119538
SM_199538._2015$tau[SM_199538._2015$tau==-999] <- NA
SM_199538._2015$doy[SM_199538._2015$doy<0] <- NA
SM_199538_2014$tau[SM_199538_2014$tau==-999] <- NA
SM_199538_2014$doy[SM_199538_2014$doy<0] <- NA
SM_199538_2013$tau[SM_199538_2013$tau==-999] <- NA
SM_199538_2013$doy[SM_199538_2013$doy<0] <- NA
SM_199538_2012$tau[SM_199538_2012$tau==-999] <- NA
SM_199538_2012$doy[SM_199538_2012$doy<0] <- NA
SM_199538_2011$tau[SM_199538_2011$tau==-999] <- NA
SM_199538_2011$doy[SM_199538_2011$doy<0] <- NA
SM_199538_2010$tau[SM_199538_2010$tau==-999] <- NA
SM_199538_2010$doy[SM_199538_2010$doy<0] <- NA

#Pickel 199538
par(mfrow=c(3,2))
plot(SM_199538._2015$doy,SM_199538._2015$tau,xaxt='n',xlab="Days",main="2015",xlim=c(59,335),ylab=expression(tau),col="red")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_199538_2014$doy,SM_199538_2014$tau,xaxt="n",xlab="Days",main="2014",xlim=c(59,335),ylab=expression(tau),col="black")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

plot(SM_199538_2013$doy,SM_199538_2013$tau,xaxt="n",xlab="Days",main="2013",xlim=c(59,335),ylab=expression(tau),col="blue")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_199538_2012$doy,SM_199538_2012$tau,xaxt="n",xlab="Days",main="2012",xlim=c(59,335),ylab=expression(tau),col="green")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

plot(SM_199538_2011$doy,SM_199538_2011$tau,xaxt="n",xlab="Days",main="2011",xlim=c(59,335),ylab=expression(tau),col="orange")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_199538_2010$doy,SM_199538_2010$tau,xaxt="n",xlab="Days",main="2010",xlim=c(59,335),ylab=expression(tau),col="brown")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

#Pixle 119546
SM_199546._2015$tau[SM_199546._2015$tau==-999] <- NA
SM_199546._2015$doy[SM_199546._2015$doy<0] <- NA
SM_199546_2014$tau[SM_199546_2014$tau==-999] <- NA
SM_199546_2014$doy[SM_199546_2014$doy<0] <- NA
SM_199546_2013$tau[SM_199546_2013$tau==-999] <- NA
SM_199546_2013$doy[SM_199546_2013$doy<0] <- NA
SM_199546_2012$tau[SM_199546_2012$tau==-999] <- NA
SM_199546_2012$doy[SM_199546_2012$doy<0] <- NA
SM_199546_2011$tau[SM_199546_2011$tau==-999] <- NA
SM_199546_2011$doy[SM_199546_2011$doy<0] <- NA
SM_199546_2010$tau[SM_199546_2010$tau==-999] <- NA
SM_199546_2010$doy[SM_199546_2010$doy<0] <- NA

#Pickel 199546
par(mfrow=c(3,2))
plot(SM_199546._2015$doy,SM_199546._2015$tau,xaxt='n',xlab="Days",main="2015",xlim=c(59,335),ylab=expression(tau),col="red")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_199546_2014$doy,SM_199546_2014$tau,xaxt="n",xlab="Days",main="2014",xlim=c(59,335),ylab=expression(tau),col="black")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

plot(SM_199546_2013$doy,SM_199546_2013$tau,xaxt="n",xlab="Days",main="2013",xlim=c(59,335),ylab=expression(tau),col="blue")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_199546_2012$doy,SM_199546_2012$tau,xaxt="n",xlab="Days",main="2012",xlim=c(59,335),ylab=expression(tau),col="green")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

plot(SM_199546_2011$doy,SM_199546_2011$tau,xaxt="n",xlab="Days",main="2011",xlim=c(59,335),ylab=expression(tau),col="orange")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_199546_2010$doy,SM_199546_2010$tau,xaxt="n",xlab="Days",main="2010",xlim=c(59,335),ylab=expression(tau),col="brown")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

#Plot 2015 for all 4 pixles
SM_200052._2015$tau[SM_200052._2015$tau==-999] <- NA
SM_200052._2015$doy[SM_200052._2015$doy<0] <- NA

par(mfrow=c(2,2))

plot(SM_199546._2015$doy,SM_199546._2015$tau,xaxt='n',xlab="Days",main="Pixle 199546 (2015)",xlim=c(59,335),ylab=expression(tau),col="red")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_199032._2015$doy,SM_199032._2015$tau,xaxt="n",xlab="Days",main="Pixle 119032 (2015)",xlim=c(59,335),ylab=expression(tau),col="black")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

plot(SM_199538._2015$doy,SM_199538._2015$tau,xaxt="n",xlab="Days",main="Pixle 119538 (2015)",xlim=c(59,335),ylab=expression(tau),col="blue")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plot(SM_200052._2015$doy,SM_200052._2015$tau,xaxt="n",xlab="Days",main="Pixle 200052 (2015)",xlim=c(59,335),ylab=expression(tau),col="green")
axis(side = 1, at = seq(60, 330, by = 30), labels =c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))


