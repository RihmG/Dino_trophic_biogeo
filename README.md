# Dinoflagellates trophic traits biogeography

Scripts and data used for obtaining results from Rihm et al., 2025 (submittted)

The data folder contains the observations drawn from metaPR2 (Vaulot et al., 2022), 
the trophic annotation data from Schneider et al. (2020) and Mitra et al. (2023), as well as the 
species list used for the species distribution modeling workflow (see folder 'sdms')
                                                              
The 'scripts' folder contains 2 types of scripts: R scripts and bash script. 

Scripts starting with 'Z' refer to the data tyding and multivariate analyses conducted prior to the sdms analysis. 
Scripts starting with 'ZM' (included in the 'sdms' folder refer to the species distribution modeling workflow. 

See the following details: 
  
  # Evenness/Multivariate analyses
  
    - Z1 script performs the trophic annotation of dinnoflagellate species.
    - Z2 (bash) and Z3 script merges the observations into unique XY coordinates. 
    - Z4 script tidies the environmental data prior to analysis. 
    - Z5 script builds the ASVs summary plot found in the supplementary materials.
    - Z6 script performs the evenness analysis and builds the Figure 1. 
    - Z7 script performs the multivariate/ordinations analyses and builds the Figure 2. 
  
  # SDMs
  
    - ZM1 creates monthly resolution species occurences from abundance matrix.
    - ZM2 creates rasterized occurence datasets.
    - ZM3 builds calibration datasets indluding rasterized occurence as well as environmental data fitted to the occurences.
    - ZM4 performs species-specific environmental predictors selection.
    - ZM5 performs individual models calibration.
    - ZM6 computes Jaccard index for individual models and filter the models based on the results.
    - ZM7 projects spatially the selected individual models and build ensemble models maps. 
    - ZM8 builds Figures 3,4 and 5 based on the models outputs. 
