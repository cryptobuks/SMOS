#Merge Dates just for 199032 pixel
library(zoo)
head(SM_199032_2010)
full <- list(SM_199032_2010, SM_199032_2011, SM_199032_2012, SM_199032_2013, SM_199032_2014, SM_199032_2015)
fulldat<-lapply(full,function(x) transform(x, doyint = as.integer(x[,1])))
#Count NA
pctna<-lapply(fulldat,function(x) sum(is.na(x[,7]))/length(x[,7]))
pctna
#Count Duplicates within pixel
dup<-lapply(fulldat, function(x) sum(duplicated(x[,7]))) 


#play with 
test1<-list_df[1]
test2 <- lapply(test1, na.omit)
vdate=as.vector(sapply(test2[1], `[[`, 1))
vtau=as.vector(sapply(test2[1], `[[`, 6))
x.zo <- zoo(vtau, vdate)
test3<-list_df[2]
test4 <- lapply(test3, na.omit)
vdate2=as.vector(sapply(test4[1], `[[`, 1))
vtau2=as.vector(sapply(test4[1], `[[`, 6))
x.zo2 <- zoo(vtau2, vdate2)
print(x.zo2,style = "vertical")

#Then make this a dataframe
m=merge(x.zo,x.zo2,all=TRUE)

y1 <- zoo(matrix(1:10, ncol = 2), 1:5)
y2 <- zoo(matrix(rnorm(10), ncol = 2), 3:7)
m=merge(x, y1, y2, all = TRUE)


a <- list(1:3, 4:5, 6:9)
 b <- c(2, 3, 5, 8)
 g <- rep(seq_along(a), sapply(a, length))
 g[match(b, unlist(a))]
 
 
 #combine matrices with union
 matLis <- list(matrix(1:4, 2, 2), matrix(1:6, 3, 2), 
                matrix(2:1, 1, 2)) 
  n <- max(sapply(matLis, nrow)) 
 do.call(cbind, lapply(matLis, function (x) 
  rbind(x, matrix(, n-nrow(x), ncol(x))))) 

