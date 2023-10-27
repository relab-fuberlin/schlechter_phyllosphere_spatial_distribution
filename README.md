# schlechter_spatial_xxxxx_2023

UNDER CONSTRUCTION

Spatial distribution analysis of syntehtic bacterial communities (SynCom) on the _Arabidopsis thaliana_ leaf surface.

This repository contains the scripts used to analyse the data published in: xxxxxx 

Raw data is storaged in zenodo:
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
git clone https://github.com/roschlec/schlechter_spatial_xxxxx_2023.git
```
Now you have all the codes to run!

## Download datasets
To download the datasets associated to the manuscript, run the following command in your terminal. You need to be in the root directory of the repository.

```
code/bash_download.sh
```

You can also download the datasets manually in Zenodo (add link)

## Changes in taxon-specific population density correlate with community complexity
Scripts used to generate docs/results1_bacdensity_communitycomplexity.Rmd
```
code/bash_script1.sh
```
This script will analyse the CFU data and create a `results` directory to store the processed data.
Additionally, the Rmarkdown containing the data analysis and plots will be rendered.

## Spatial distribution of individual strains depends on their community context
Scripts used to generate docs/results2_celldensity_communitycomplexity.Rmd
```
code/bash_script2.sh
```

## Effect of community complexity on intraspecific spatial relations
Scripts used to generate docs/results3_Kest_communitycomplexity.Rmd
```
code/bash_script3.sh
```

## Effect of community complexity on interspecific spatial correlations
Scripts used to generate docs/results4_PCF_communitycomplexity.Rmd
```
code/bash_script4.sh
```
