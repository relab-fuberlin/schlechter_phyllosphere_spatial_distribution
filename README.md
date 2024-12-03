# Differential responses of _Methylobacterium_ and _Sphingomonas_ species to multispecies interactions in the phyllosphere
Spatial distribution analysis of synthetic bacterial communities (SynCom) on the _Arabidopsis thaliana_ leaf surface.
This repository contains the scripts used to analyse the data published in Environmental Microbiology (accepted). Raw data is available in [Zenodo](https://zenodo.org/doi/10.5281/zenodo.100361160):
File | Data type | Content
:---: | :---: | :---:
`cfu.csv` | CSV table | CFU data of each bacterial population in arabidopsis
`comm_id.csv` | CSV table | Composition of each SynCom
`coordinates_S2.csv` | CSV table | Coordinates of cells within S2
`coordinates_S3.csv` | CSV table | Coordinates of cells within S3
`metadata.csv` | CSV table | Data of each community including biological replicates and images taken
`bacimg.tar.gz` | Images | Images taken for cell populations in C, S2 and S3
`Supplemental Material.pdf` | PDF | Supplemental Material
## Setting up 
First, you need to clone this repository
```
git clone https://github.com/roschlec/schlechter_phyllosphere_spatial_distribution.git
```
Now you have all the codes to run!

For reproducibility in data analysis, a renv.lock file is available to install and load all R packages (and corresponding versions) used in this paper.

In R:
```
install.packages("renv")
library(renv)
renv::restore()
```
## Download datasets
To download the datasets associated to the manuscript, run the following command in your terminal.
First, you need to install zenodo-get to download the files. Documentation about zenodo-get [here](https://gitlab.com/dvolgyes/zenodo_get).
You don't have to install this package if you use the conda environment provided.
```
pip install zenodo-get
```
Now you should be able to run the following code from the root directory
```
code/bash_download.sh
```
You can also download the datasets manually in [Zenodo](https://zenodo.org/doi/10.5281/zenodo.100361160).
## Changes in taxon-specific population density correlate with community complexity
Scripts used to generate `docs/01_CFU.Rmd`
```
code/bash_script1.sh data
```
This script will analyse the CFU data and create a `results` directory to store the processed data.
Additionally, the Rmarkdown containing the data analysis and plots will be rendered.
## Spatial distribution of individual strains depends on their community context
Scripts used to generate `docs/02_celldensity.Rmd`
```
code/bash_script2.sh
```
This script will analyse the single-cell data. Additionally, the Rmarkdown containing the data analysis and plots will be rendered.
## Effect of community complexity on intraspecific spatial relations
Scripts used to generate `docs/03_Intraspecific_spatial.Rmd`
```
code/bash_script3.sh
```
Using coordinate data, this script will perform the spatial analysis (K-estimates) within each population. Additionally, the Rmarkdown containing the data analysis and plots will be rendered.
## Effect of community complexity on interspecific spatial correlations
Scripts used to generate `docs/04_Interspecific_spatial.Rmd`
```
code/bash_script4.sh
```
Using coordinate data, this script will perform the spatial analysis (pair cross-correlation) between populations. Additionally, the Rmarkdown containing the data analysis and plots will be rendered.
## Resubmission
```
code/bash_resubmission.sh
```
