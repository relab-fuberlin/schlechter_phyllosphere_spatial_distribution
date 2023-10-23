# schlechter_spatial_xxxxx_2023
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
`s2.tar.gz` | Images | Images taken for cell populations within S2
`s3.tar.gz` | Images | Images taken for cell populations within S3

## Changes in taxon-specific population density correlate with community complexity
Scripts used to generate docs/results1_bacdensity_communitycomplexity.Rmd

## Spatial distribution of individual strains depends on their community context
Scripts used to generate docs/results2_celldensity_communitycomplexity.Rmd

## Effect of community complexity on intraspecific spatial relations
Scripts used to generate docs/results3_Kest_communitycomplexity.Rmd

## Effect of community complexity on interspecific spatial correlations
Scripts used to generate docs/results4_PCF_communitycomplexity.Rmd
