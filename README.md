# nSBM-breast-cancer
Repository of the paper in preparation (Malagoli et al. 2022), preprint available here: 
https://www.biorxiv.org/content/10.1101/2022.09.13.507779v1

To reproduce our experiments you need to create a Anaconda environment following the below instructions. In this GitHub you will find the notebooks and the necessary helper files. All the data and the results can be downloaded from https://drive.google.com/drive/folders/1Tdwt9Z85ru6r1yrDnA7WyhBg8e3WNHdW?usp=sharing. You need to download all the material (from GitHub and the link) locally on your machine.

# Presentation
We propose a topic modeling approach to investigate the transcriptional heterogeneity of luminal and triple negative breast cancer cells using patient-derived xenograft models of acquired resistance to chemotherapy and targeted therapy. 
We show that using an integrative clustering that combines the information coming from mRNAs and lncRNAs treated as disjoint omic layers greatly improves the accuracy of cell classification. Topics associated with specific breast cancer subpopulations show a clear enrichment for pathways involved in subtyping and progression of breast cancer and to sets of lncRNA encoded in the open chromatin regions of breast cancer cell lines.



## HOWTO create the appropriate environment to reproduce our experiments

### 1) Install [graph-tool](https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions) as new conda environment

```
conda create --name malagolietal2022 -c conda-forge graph-tool
```


### 2) Activate the env

```
conda activate malagolietal2022
```


### 3) Install additional packages

```
conda install -n malagolietal2022 -c conda-forge ipython jupyter
```


### 4) Install [nSBM](https://github.com/BioPhys-Turin/nsbm)

```
conda install -c conda-forge nsbm
```


### 5) Install the requirements file

```
pip install -r requirements.txt
```


--------------------------------------------------


### You can simply create an environment using th .yml file you can find in this Github using the command line

```
conda env create -f /PATH_TO_THE_FILE/malagoli_et_al_2022.yml -n malagolietal2022
```

### The activate the environment 

```
conda activate malagolietal2022
```

# References

Our work is mainly based on technical findings of

Gerlach, M., Peixoto, T. P., & Altmann, E. G. (2018). A network approach to topic models. Science Advances, 4(7), eaaq1360. https://doi.org/10.1126/sciadv.aaq1360

C. C. Hyland, Y. Tao, L. Azizi, M. Gerlach, T. P. Peixoto, and E. G. Altmann Multilayer Networks for Text Analysis with Multiple Data Types , EPJ Data Science 10, 33 (2021). https://epjdatascience.springeropen.com/articles/10.1140/epjds/s13688-021-00288-5
Setup

https://github.com/martingerlach/hSBM_Topicmodel

and 

Valle, F.; Osella, M.; Caselle, M. Multiomics Topic Modeling for Breast Cancer Classification. Cancers 2022, 14, 1150. https://doi.org/10.3390/cancers14051150

https://github.com/BioPhys-Turin/nsbm

and on the experimental data of

Grosselin, K., Durand, A., Marsolier, J. et al. High-throughput single-cell ChIP-seq identifies heterogeneity of chromatin states in breast cancer. Nat Genet 51, 1060–1066 (2019). https://doi.org/10.1038/s41588-019-0424-9


# License
See [LICENSE](https://github.com/gmalagol10/nSBM-breast-cancer/blob/main/LICENSE)
