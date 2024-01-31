# uRF
<p align="center">
<img src="https://github.com/pievos101/uRF/blob/main/uRF.jpg" width="400">
</p>


## Installation
The uRF R-package can be installed using devtools.

```{r}
install.packages("devtools")
library(devtools)
devtools::install_github("pievos101/uRF")
library(uRF)

```

## Usage

### Clustering

```{r}
library(uRF)

data(iris)
d = iris[,1:4]

# RF-based clustering
# Here, we just use one view, whereas it can be a list of views
res = uRF::uRF(list(d), k=3)

res$affinity
res$clusters

# Performance of the clustering
library(aricode)
ARI(res$clusters, iris[,5])
NMI(res$clusters, iris[,5])
```

### Feature Importance per Cluster

```{r}
# Now, that we have the clusters computed, 
# feature importance per cluster and feature can be determined
fimp = uRF::fimp(d, clusters=res$clusters)
fimp

```

### Federated Computation
```{r}

library(uRF)
data(iris)
d = iris[,1:4]

# Distribute data over three clients
n.clients = 3
samp = as.numeric(rownames(d))
df = sample(samp, length(samp))
spl = split(df, rep(1:n.clients, length.out = length(df), each = ceiling(length(df)/3)))

d1 = d[spl[[1]],]
d2 = d[spl[[2]],]
d3 = d[spl[[3]],]

res_1 = uRF::uRF(list(d1), k=3)
res_2 = uRF::uRF(list(d2), k=3)
res_3 = uRF::uRF(list(d3), k=3)

res_1_global = uRF::feduRF(list(d1), 
        models=list(res_1$model,res_2$model,res_3$model), k=3)
res_2_global = uRF::feduRF(list(d2), 
        models=list(res_1$model,res_2$model,res_3$model), k=3)
res_3_global = uRF::feduRF(list(d3), 
        models=list(res_1$model,res_2$model,res_3$model), k=3)


library(aricode)
print("Local performance (non-federated model)")
ARI(res_1$clusters, iris[spl[[1]],5])
ARI(res_2$clusters, iris[spl[[2]],5])
ARI(res_3$clusters, iris[spl[[3]],5])

print("Local performance (federated model)")
ARI(res_1_global$clusters, iris[spl[[1]],5])
ARI(res_2_global$clusters, iris[spl[[2]],5])
ARI(res_3_global$clusters, iris[spl[[3]],5])

```


## Citation
If you find the uRF package useful please cite

```
@misc{pfeifer2024federated,
      title={Federated unsupervised random forest for privacy-preserving patient stratification}, 
      author={Bastian Pfeifer and Christel Sirocchi and Marcus D. Bloice and Markus Kreuzthaler and Martin Urschler},
      year={2024},
      eprint={2401.16094},
      archivePrefix={arXiv},
      primaryClass={cs.LG}
}
```

## Acknowledgement 
The uRF package is based on the ranger package (https://github.com/imbs-hl/ranger).

