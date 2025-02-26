
# ------ Calibration/Validation datasets building ---------


# The goal is to create dataframes with coordinates data, month information, and environmental variables
# fitted for each site

# In the end we obtain one dataframe per species, that gather all env. information for max 12 months. 

library(tidyverse)
library(tidyr)  
library(terra)


  for(m in 1:12)
  {
    # We load the monthly abundance matrix
    mothermat <- readRDS(paste0("month_occ/matrix/mothermat_", m, ".RDS")) 
    # Extract the rowxnames to get the species present in this month
    row.names(mothermat)[row.names(mothermat) == "Heterocapsa_nei/rotundata"] <- "Heterocapsa_nei_rotundata"
  
    sp_list <- rownames(mothermat)# list of species
    
    base <- readRDS(paste0("data_out/Env_data/monthly_base/base_",m, ".RDS"))# environnmental stack 
    
    for(sp in sp_list)
    {
      
      if(!dir.exists(paste0("month_occ/species_df/",sp, "/")))
      {
        dir.create(paste0("month_occ/species_df/", sp, "/"), recursive = T)
      }
      
        sp_occ <- readRDS(paste0("month_occ/merged_presabs/", m, "/", sp, m, ".RDS"))# load occurence data
        sp_occ$month <- m # add month information to the df
        
        env_val <- extract(base, sp_occ[, c("x", "y")]) # extraction of env. values for each locations of occurences
        sp_occ <- cbind(sp_occ, env_val) # merging with occurences dataframe 
        colnames(sp_occ) <- c("x","y", "ap", "month", "ID", "SST", "SSS", "MLD", "PAR", 
                              "O2", "Sistar", "Nstar", "no3", "po4", "si", "chla", "pic")
        saveRDS(sp_occ, paste0("month_occ/species_df/", sp, "/", sp, "_", m, ".RDS"))
      
    }
  }

  # =================================================================================================







