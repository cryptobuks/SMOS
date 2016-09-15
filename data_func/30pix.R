#Bring data files back in and Add Variable names
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
