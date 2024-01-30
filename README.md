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

```{r}
library(uRF)

data(iris)
d = iris[,1:4]

# RF-based clustering
# Here, we just use one view, whereas it can be a list of views
res = uRF::uRF(list(d), k=3)

res$affinity
res$clusters
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

