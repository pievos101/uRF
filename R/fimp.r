fimp <- function(omic = NaN, clusters=NaN, mb=5, ntrees=100){

cl1 = clusters
n.cluster = length(unique(cl1))
cl = unique(cl1)

FIMP  = matrix(NaN, n.cluster, dim(omic)[2])
   

for (xx in 1:n.cluster){
   
   # one vs all
   resX = cl1
   ids = is.element(cl1, cl[-xx])
   resX[ids] = cl[xx] + 1 

   ###########
   # Train a random forest classifier
   model = ranger(x=omic, y=as.factor(resX), 
            num.trees=ntrees, 
            #classification=NULL,
            clustering = TRUE, 
            probability=FALSE, 
            importance='impurity',
            #max.depth=10, 
            #min.node.size = 3, 
            mtry = 2, #ncol(train),
            #replace=TRUE,
            min.bucket=mb,
            oob.error = FALSE)


FIMP[xx,] = model$variable.importance

}

rownames(FIMP) = paste("cluster",1:n.cluster)
colnames(FIMP) = colnames(omic)

return(FIMP)

}