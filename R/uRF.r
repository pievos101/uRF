# unsupervised random forest function
uRF <- function(omics = list(), k=NaN, mb=5, ntrees=100){

ntrees = ntrees
COUNTS = list()

for (z in 1:length(omics)){

train = omics[[z]]

###########
# Train a random forest classifier
model = uRF::ranger(x=train, #y=as.factor(target), 
            num.trees=ntrees, 
            #classification=NULL,
            clustering = TRUE, 
            probability=FALSE, 
            importance='impurity',
            #max.depth=maxd, 
            #min.node.size = 20, 
            mtry = 2, #ncol(train),# ceiling(ncol(train)/2),
            #replace=FALSE,
            min.bucket=mb, #30 in case of real-world cancer data
            oob.error = FALSE)

# Predictions on test data
#pred = predict(model, train)$predictions
#pred = apply(pred,1,function(x){which.max(x)-1})
pred.ranger <- predict(model, train, type = "terminalNodes")
lyves <- pred.ranger$prediction
 
# Get Similarity
SIMILARITY = matrix(0, dim(lyves)[1], dim(lyves)[1])

for (i in 1:(dim(SIMILARITY)[1]-1)){
 for (j in i:dim(SIMILARITY)[1]){
    res = lyves[c(i,j),]
    hit = sum(apply(res,2,function(x){x[1]==x[2]}))
    SIMILARITY[i,j] = hit
    SIMILARITY[j,i] = hit
    
 }
}

COUNTS[[z]] = SIMILARITY

}

SIMILARITY = Reduce("+", COUNTS)
DISTANCE_ALL = 1 - SIMILARITY/max(SIMILARITY)
DISTANCE_ALL = as.dist(DISTANCE_ALL)

if(is.na(k)){
   min.nc = 2
   max.nc = 10

   require(NbClust)
   bestCl = NbClust(distance=NULL, diss=DISTANCE_ALL, 
   method="ward.D2", index="silhouette", min.nc=min.nc,
         max.nc=max.nc)$Best.nc[1]

   hc = hclust(DISTANCE_ALL, method="ward.D2")
   cl = cutree(hc, bestCl)
}else{
   hc = hclust(DISTANCE_ALL, method="ward.D2")
   cl = cutree(hc, k)
}

#print(bestCl)
set.seed(Sys.time())
return(list(clusters=cl, affinity=SIMILARITY))
}


uRF2 <- function(omics = list()){

ntrees = 500
COUNTS = list()

model_all = list()

for (z in 1:length(omics)){

train = omics[[z]]

###########
# Train a random forest classifier
model = ranger(x=train, #y=as.factor(target), 
            num.trees=ntrees, 
            #classification=NULL,
            clustering = TRUE, 
            probability=FALSE, 
            importance='impurity',
            #max.depth=maxd, 
            #min.node.size = 20, 
            mtry = 2, #ncol(train),# ceiling(ncol(train)/2),
            #replace=FALSE,
            min.bucket=5, #30 in case of real-world cancer data
            oob.error = FALSE)

model_all[[z]] = model

# Predictions on test data
#pred = predict(model, train)$predictions
#pred = apply(pred,1,function(x){which.max(x)-1})
pred.ranger <- predict(model, train, type = "terminalNodes")
lyves <- pred.ranger$prediction
 
# Get Similarity
SIMILARITY = matrix(0, dim(lyves)[1], dim(lyves)[1])

for (i in 1:(dim(SIMILARITY)[1]-1)){
 for (j in i:dim(SIMILARITY)[1]){
    res = lyves[c(i,j),]
    hit = sum(apply(res,2,function(x){x[1]==x[2]}))
    SIMILARITY[i,j] = hit
    SIMILARITY[j,i] = hit
    
 }
}

COUNTS[[z]] = SIMILARITY

}

SIMILARITY = Reduce("+", COUNTS)
DISTANCE_ALL = 1 - SIMILARITY/max(SIMILARITY)
DISTANCE_ALL = as.dist(DISTANCE_ALL)

min.nc = 2
max.nc = 10

bestCl = NbClust(distance=NULL, diss=DISTANCE_ALL, 
  method="ward.D2", index="silhouette", min.nc=min.nc,
        max.nc=max.nc)$Best.nc[1]

hc = hclust(DISTANCE_ALL, method="ward.D2")
cl = cutree(hc, bestCl)

print(bestCl)
set.seed(Sys.time())
return(list(cl=cl, model_all=model_all))
}