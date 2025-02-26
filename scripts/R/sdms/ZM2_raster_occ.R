
# ------ SPECIES OCCURRENCE RASTERIZATION ------


# This script takes as input the abundance matrix generated in script Z3 and creates specific occurence data
# for each modeled species


library(tidyverse)
library(tidyr)  
library(terra)
library(Rarity)
library(viridis)
library(reshape2)
library(plyr)


# ----- 1. DATA -------

    mothermat <- readRDS("data_out/mothermat.RDS") # abundance matrix 
    
    row.names(mothermat)[row.names(mothermat) == "Heterocapsa_nei/rotundata"] <- "Heterocapsa_nei_rotundata"
    
    test_rasta <- rast(rast("data_out/Env_data/baseline_logmean.tif"))
    
    
# ----- 2. MONTH ABUNDANCE RASTERIZATION --------
   
 for (m in 1:12)
 {

    #test_rasta = template raster
    
    mothermat <- readRDS(paste0("month_occ/matrix/mothermat_", m, ".RDS")) 
    
    row.names(mothermat)[row.names(mothermat) == "Heterocapsa_nei/rotundata"] <- "Heterocapsa_nei_rotundata"
    
    test_rasta <- rast("data_out/Env_data/baseline_logmean.tif")
    
    
    if(!dir.exists(paste0("month_occ/sp_rast/", m, "/")))
    {
      dir.create(paste0("month_occ/sp_rast/", m, "/"), recursive = T)
    }
    
    if(!dir.exists(paste0("month_occ/sp_rast_maps/", m, "/")))
    {
      dir.create(paste0("month_occ/sp_rast_maps/", m, "/"), recursive = T)
    }
    
    if(!dir.exists(paste0("month_occ/sp_rast_abs/", m, "/")))
    {
      dir.create(paste0("month_occ/sp_rast_abs/", m, "/"), recursive = T)
    }
    
    if(!dir.exists(paste0("month_occ/sp_rast_maps_abs/", m, "/")))
    {
      dir.create(paste0("month_occ/sp_rast_maps_abs/", m, "/"), recursive = T)
    }
    
    for(i in 1:length(rownames(mothermat)))
    {
      
      sp <- rownames(mothermat)[i]
      
      sp_grd <- readRDS(paste0("month_occ/sp_sites/",m,"/pres_", sp, ".RDS"))
      sp_grd <- sp_grd[,c(1,2)]
      sp_grd$pres <- 1
      sp_grd$pres <- as.integer(sp_grd$pres)
      
      sp_grd_vect <- vect(sp_grd, 
                          geom = c("longitude", "latitude"), 
                          crs = "EPSG:4326")
      
      sp_grd_rast <- terra::rasterize(x = sp_grd_vect,
        #x = as.matrix(sp_grd[,1:2][, c("longitude", "latitude")]), 
                               y = test_rasta,
                               fun = sum)
      
      sp_grd_rast[sp_grd_rast > 1] <- 1
      
      sp_grd.df <- as.data.frame(sp_grd_rast, xy = T)
      
      #writeRaster(sp_grd_rast, paste0("data_out/sp_sites_rast_tif/", sp, "_pres.tif"))
      saveRDS(sp_grd.df, paste0("month_occ/sp_rast/", m,"/pres_", sp, m, ".RDS"))
      
      p <- ggplot() +
        geom_sf(data = worldmap, color = "grey", fill = "grey") +
        theme(panel.background = element_rect(fill = "aliceblue")) +
        geom_point(sp_grd.df, mapping = aes(x = x, y = y), size = 1) + 
        ggtitle(sp)
      
      cat("-----PRES-----", 
          "\n", sp, " ", m)
      
      ggsave(paste0("month_occ/sp_rast_maps/", m , "/",sp, m,".pdf"))
      
    }
    
    # We also want to create a df with all the site locations, but after rasterization
    
    
    # We do the rasterization with Absence Data (from abs_ file in /data_out)
    
    for(i in 1:length(rownames(mothermat)))
    {
      sp <- rownames(mothermat)[i]
      
      sp_grd <- readRDS(paste0("month_occ/sp_sites_abs/",m ,"/abs_", sp, ".RDS"))  
      sp_grd <- sp_grd[,c(1,2)] 
      
      sp_grd_abs_vect <- vect(sp_grd, 
                              geom = c("longitude", "latitude"), 
                              crs = "EPSG:4326")
      
      sp_grd_abs_rast <- rasterize(x = sp_grd_abs_vect, 
                                   y = test_rasta,
                                   fun = sum)
      
      sp_grd_abs_rast[sp_grd_abs_rast >= 1] <- 0
      names(sp_grd_abs_rast) <- paste0("abs_",sp)
      
      sp_grd_abs.df <- as.data.frame(sp_grd_abs_rast, xy = T)
      
      #writeRaster(sp_grd_rast, paste0("data_out/sp_sites_rast_tif/", sp, "_pres.tif"))
      saveRDS(sp_grd_abs.df, paste0("month_occ/sp_rast_abs/", m, "/abs_", sp, m, ".RDS"))
      
      p <- ggplot() +
        geom_sf(data = worldmap, color = "grey", fill = "grey") +
        theme(panel.background = element_rect(fill = "aliceblue")) +
        geom_point(sp_grd_abs.df, mapping = aes(x = x, y = y), size = 1) + 
        ggtitle(sp)
      
      cat("-----ABS-----", 
          "\n", sp, " ", m)
      
      ggsave((paste0("month_occ/sp_rast_maps_abs/", m, "/abs_", sp, m, ".pdf")))

    }
  } 



# Need to check for double Presence/Absence due to the previous rasterization process
# If there are pixels that are presence and absences at the same time, we consider it is presence
for (m in 1:12)
  {
  
  if(!dir.exists(paste0("month_occ/merged_presabs/", m, "/")))
  {
    dir.create(paste0("month_occ/merged_presabs/", m, "/"), recursive = T)
  }
  
  mothermat <- readRDS(paste0("month_occ/matrix/mothermat_", m, ".RDS")) 
  
  row.names(mothermat)[row.names(mothermat) == "Heterocapsa_nei/rotundata"] <- "Heterocapsa_nei_rotundata"
  
  sp_list <- rownames(mothermat)
  
    for (sp in sp_list) 
    {
      
      sp_pres <- readRDS(paste0("month_occ/sp_rast/",m,"/pres_" ,sp,m,".RDS"))
      sp_abs <- readRDS(paste0("month_occ/sp_rast_abs/",m,"/abs_" ,sp,m,".RDS"))
      
      colnames(sp_abs) <- c("x", "y", "ap")
      colnames(sp_pres) <- c("x", "y", "ap")
      
      # We merge the PRESENCE and ABSENCE datasets
      
      sp.f <- rbind(sp_pres, sp_abs)
      
      # we test if there are double coordinates for absence and presence 
      
      wlf <- sp.f[,c(1,2)] # take only spatial coordinates
      dup_rows <- duplicated(wlf) | duplicated(wlf, fromLast = TRUE) # extract duplicated rows 
      df_dup <- wlf[dup_rows,]
      df_nodup <- wlf[!dup_rows,] # extract single rows 
      
      # We create 2 df that we will recombine after removing duplicated rows 
      sp.f1 <- sp.f[rownames(df_nodup),]
      sp.f2 <- sp.f[rownames(df_dup),]
      
      # We remove duplicated rows keeping only pixels with the presence (not the absence)
      sp.f2 %>%
        filter(ap == 1) %>%
        distinct(x, y, .keep_all = TRUE) -> sp.f2
      
      # We merge the 2 df 
      sp.ff <- rbind(sp.f1, sp.f2) 
      
      # final presence/absence data
      sp.ff <- sp.ff[order(sp.ff$ap),]
      
      saveRDS(sp.ff, paste0("month_occ/merged_presabs/",m,"/" ,sp, m, ".RDS"))
      
    }
}


# ============================================================================================





