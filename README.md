# Bacterial community complexity in the phyllosphere penalises specialists over generalists

Spatial distribution analysis of syntehtic bacterial communities (SynCom) on the _Arabidopsis thaliana_ leaf surface.

This repository contains the scripts used to analyse the data published in: xxxxxx 

Raw data is storaged in [Zenodo](https://zenodo.org/doi/10.5281/zenodo.100361160):
File | Data type | Content
:---: | :---: | :---:
`cfu.csv` | CSV table | CFU data of each bacterial population in arabidopsis
`comm_id.csv` | CSV table | Composition of each SynCom
`coordinates_S2.csv` | CSV table | Coordinates of cells within S2
`coordinates_S3.csv` | CSV table | Coordinates of cells within S3
`metadata.csv` | CSV table | Data of each community including biological replicates and images taken
`bacimg.tar.gz` | Images | Images taken for cell populations in C, S2 and S3

## Setting up 
First, you need to clone this repository
```
git clone https://github.com/roschlec/schlechter_phyllosphere_spatial_distribution.git
```
Now you have all the codes to run!

It is advisable to run these codes in a conda environment. The file `environment.yml` will help you create an environment will the necessary packages. We recommend using [mamba](https://mamba.readthedocs.io/en/latest/index.html).

```
mamba env create -f environment.yml
conda activate spatial
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
Scripts used to generate docs/results1_bacdensity_communitycomplexity.Rmd
```
code/bash_script1.sh data
```
This script will analyse the CFU data and create a `results` directory to store the processed data.
Additionally, the Rmarkdown containing the data analysis and plots will be rendered.

## Spatial distribution of individual strains depends on their community context
Scripts used to generate docs/results2_celldensity_communitycomplexity.Rmd
```
code/bash_script2.sh
```
This script will analyse the single-cell data. Additionally, the Rmarkdown containing the data analysis and plots will be rendered.

## Effect of community complexity on intraspecific spatial relations
Scripts used to generate docs/results3_Kest_communitycomplexity.Rmd
```
code/bash_script3.sh
```
Using coordinate data, this script will perform the spatial analysis (K-estimates) within each population. Additionally, the Rmarkdown containing the data analysis and plots will be rendered.

## Effect of community complexity on interspecific spatial correlations
Scripts used to generate docs/results4_PCF_communitycomplexity.Rmd
```
code/bash_script4.sh
```
Using coordinate data, this script will perform the spatial analysis (pair cross-correlation) between populations. Additionally, the Rmarkdown containing the data analysis and plots will be rendered.
