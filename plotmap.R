library(maps)

#grab first row from a list of data.frames

maps2=lapply(series,"[",1,2:4,drop=FALSE)
DF.maps <- as.data.frame(
  do.call(rbind.fill.matrix,
          lapply(maps2, function(l) {
            res <- unlist(l)
            t(res)
          })
  )
)


map('county', 'iowa',interior = T)
data(us.cities)
map.cities(us.cities, country="IA",label=TRUE)

points(DF.maps$lon,DF.maps$lat, col="blue")
text(DF.maps$lon, DF.maps$lat, DF.maps$dgg, col="red", cex=.5,pos=3)
