setwd("C:/Users/Colin/Documents/Caragea/")
#Pull Duplicate doy
#Pull All Duplicate Doy
dups <- vector("list", length(list_df))
for(i in 1:180){
  dups[[i]]=list_df[[i]] %>% group_by(doy) %>% filter(n() >= 2) %>% filter(!is.na(doy))
}

duplicates = do.call("rbind", dups)
dup2 <- duplicates[duplicates$tau>0,]
dup3<-na.omit(dup2)
dup.doy= dup3 %>% group_by(doy,dgg) %>% filter(n() >= 2) 

write.csv(dup.doy, file = "duplicatedoy.csv",row.names = TRUE)
