feduRF <- function(omics = list(), models=list(), k=NaN){

# concatenate models
lyves = vector("list", length(models))

for(xx in 1:length(models)){

lyves[[xx]] = predict(models[[xx]], omics[[1]], 
    type = "terminalNodes")$prediction

}

LYVES = Reduce("cbind", lyves)


# Get Similarity
SIMILARITY = matrix(0, dim(LYVES)[1], dim(LYVES)[1])

for (i in 1:(dim(SIMILARITY)[1]-1)){
 for (j in i:dim(SIMILARITY)[1]){
    res = LYVES[c(i,j),]
    hit = sum(apply(res,2,function(x){x[1]==x[2]}))
    SIMILARITY[i,j] = hit
    SIMILARITY[j,i] = hit
 }
}

DISTANCE_ALL = 1 - SIMILARITY/max(SIMILARITY)
DISTANCE_ALL = as.dist(DISTANCE_ALL)

if(is.na(k)){
    min.nc = 2
    max.nc = 10
    bestCl = NbClust(distance=NULL, diss=DISTANCE_ALL, 
    method="ward.D2", index="silhouette", min.nc=min.nc,
            max.nc=max.nc)$Best.nc[1]

    hc = hclust(DISTANCE_ALL, method="ward.D2")
    cl1_global = cutree(hc, bestCl)
}else{
    hc = hclust(DISTANCE_ALL, method="ward.D2")
    cl1_global = cutree(hc, k)
}

return(list(clusters=cl1_global, affinity=SIMILARITY))
}