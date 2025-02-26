
#------ Environmental Data --------

# This script creates raster stacks (annual and month stacks) usable for 
# following multivariate analysis and SDMs

  library(raster)
  library(terra)
  library(sf)


# ----- 1. Data -----

    rast_names = list.files("Env/raster_files/", pattern = "*.grd")  
    
    for (i in 1:length(rast_names)){
      assign(paste0("stack_", str_sub(rast_names[i], 17,30)), stack(paste0("Env/raster_files/",rast_names[i])))
    }
    
    rm(raster_mon_clim_bbp443_22_03_22.grd,raster_mon_clim_CHLA_22_03_22.grd,raster_mon_clim_DIC_22_03_22.grd,raster_mon_clim_Dist2coast_22_03_22.grd,raster_mon_clim_EddyKineticEnergy_22_03_22.grd,raster_mon_clim_Kd490_22_03_22.grd,raster_mon_clim_MLD_22_03_22.grd,raster_mon_clim_Nitrates_22_03_22.grd,raster_mon_clim_Nitrates_mld_22_03_22.grd,raster_mon_clim_Nitrates_mld.grad_22_03_22.grd,raster_mon_clim_O2_200m_22_03_22.grd,raster_mon_clim_O2_22_03_22.grd,raster_mon_clim_O2_mld_22_03_22.grd,raster_mon_clim_Omega_Aragonite_22_03_22.grd,raster_mon_clim_Omega_Calcite_22_03_22.grd,raster_mon_clim_PAR_22_03_22.grd,raster_mon_clim_pCO2_22_03_22.grd,raster_mon_clim_Phosphates_22_03_22.grd,raster_mon_clim_Phosphates_mld_22_03_22.grd,raster_mon_clim_Phosphates_mld.grad_22_03_22.grd,raster_mon_clim_PIC_22_03_22.grd,raster_mon_clim_Revelle_22_03_22.grd,raster_mon_clim_Sal_mld_22_03_22.grd,raster_mon_clim_Silicates_22_03_22.grd,raster_mon_clim_Silicates_mld_22_03_22.grd,raster_mon_clim_Silicates_mld.grad_22_03_22.grd,raster_mon_clim_SSS_22_03_22.grd,raster_mon_clim_SST_22_03_22.grd,raster_mon_clim_TAlk_22_03_22.grd,raster_mon_clim_Temp_mld_22_03_22.grd,raster_mon_clim_Temp_mld.grad_22_03_22.grd,raster_mon_clim_Wind.speed_22_03_22.grd,raster_mon_clim_Zeu_22_03_22.grd
    )
    
    paste0("stack_", str_sub(rast_names, 17,30) ,",")
    
    
    
    longhurst <- st_read("Env/Longhurst_world_v4_2010/Longhurst_world_v4_2010.shp")
    

    #we rename rasterlayers
    stack_CHLA <- stack_CHLA_22_03_22.
    stack_MLD <- stack_MLD_22_03_22.g
    stack_NO3 <- stack_Nitrates_22_03
    stack_O2 <- stack_O2_22_03_22.gr
    stack_PO4 <- stack_Phosphates_22_
    stack_Si <- stack_Silicates_22_0
    stack_SSS <- stack_SSS_22_03_22.g
    stack_SST <- stack_SST_22_03_22.g
    stack_PAR <- stack_PAR_22_03_22.g

    # We rename Layers inside the stack
    names(stack_SST) <- paste0("SST_", 01:12)
    names(stack_SSS) <- paste0("SSS_", 01:12)
    names(stack_PAR) <- paste0("PAR_", 01:12)
    names(stack_CHLA) <- paste0("CHLA_", 01:12)
    names(stack_MLD) <- paste0("MLD_", 01:12)
    names(stack_NO3) <- paste0("NO3_", 01:12)
    names(stack_O2) <- paste0("O2_", 01:12)
    names(stack_PO4) <- paste0("PO4_", 01:12)
    names(stack_Si) <- paste0("Si_", 01:12)


# ----- 2. Annually averaged Rasters -----
    
      # for computing analyses before modeling, it would be interesting to have annual averaged 
      #raster layers
      mean_sst <- mean(stack_SST)
      mean_sss <- mean(stack_SSS)
      mean_mld <- mean(stack_MLD)
      mean_po4 <- mean(stack_PO4)
      mean_no3 <- mean(stack_NO3)
      mean_si <- mean(stack_Si)
      mean_chla <- mean(stack_CHLA)
      mean_par <- mean(stack_PAR)
      mean_O2 <- mean(stack_O2)
  
  
      mean_par2 <- mean(stack_PAR, na.rm = T) # keep high latitude values 
      
      
      # update 07/01, addition of Particulate inorganic carbon
      stack_pic <- stack_PIC_22_03_22.g
      
      names(stack_pic) <- paste0("PIC_", 01:12)
      mean_pic <- mean(stack_pic)
      
      # and also calculation of nutrients indices 
      mean_Nstar <- mean_no3 - 16*mean_po4
      
      mean_Sistar <- mean_si - mean_no3
      
      stack_Nstar = stack_NO3 - 16*stack_NO3
      stack_Sistar = stack_Si - stack_NO3
      
      names(stack_Nstar) <- paste0("Nstar_", 01:12)
      names(stack_Sistar) <- paste0("Sistar_", 01:12)
      
  
    # We create a final stack 
    baseline <-  c(as(stack_SST, "SpatRaster"), 
                   as(stack_SSS, "SpatRaster"),
                   as(stack_MLD, "SpatRaster"),
                   as(stack_CHLA, "SpatRaster"),
                   as(stack_NO3, "SpatRaster"),
                   as(stack_PO4, "SpatRaster"),
                   as(stack_Si, "SpatRaster"),
                   as(stack_PAR, "SpatRaster"),
                   as(stack_pic, "SpatRaster"),
                   as(stack_O2, "SpatRaster"),
                   as(stack_Sistar, "SpatRaster"),
                   as(stack_Nstar, "SpatRaster"))
  
  
    writeRaster(baseline, "data_out/Env_data/baseline.tif", overwrite = T)
  
    
    baseline_mean <- stack(mean_sst,
                           mean_sss,
                           mean_mld,
                           mean_chla,
                           mean_no3,
                           mean_po4,
                           mean_si,
                           mean_par,
                           mean_pic,
                           mean_O2,
                           mean_Sistar,
                           mean_Nstar)
    
  names(baseline_mean) <- c("sst", "sss", "mld", "chla", "no3", "po4", "si", "par", "pic", "O2", "Sistar", "Nstar")
  
  raster::writeRaster(baseline_mean, "data_out/Env_data/baseline_mean.tif", overwrite = T)
  


  # Log Transformation 
  
    # We change the 0 values into extremely low values in order not to 
    # obtain -Inf values in log-transformed rasters
    
    stack_NO3[stack_NO3 == 0] <- 1e-10
    stack_PO4[stack_PO4 == 0] <- 1e-10
    stack_Si[stack_Si == 0] <- 1e-10
    stack_CHLA[stack_CHLA == 0] <- 1e-10
    
    # As nutrients rasters are unbalanced towards lower values, we use a logarithmic transformation to scale it 
    
    stack_logno3 <- log(stack_NO3)
    stack_logpo4 <- log(stack_PO4)
    stack_logsi <- log(stack_Si)
    stack_logchla <- log(stack_CHLA)
    stack_logpic <- log(stack_pic)
    
    names(stack_logno3) <- paste0("no3_", 01:12)
    names(stack_logpo4) <- paste0("po4_", 01:12)
    names(stack_logsi) <- paste0("si_", 01:12)
    names(stack_logchla) <- paste0("chla_", 01:12)
    names(stack_logpic) <- paste0("pic_", 01:12)
    
    
    
    mean_logno3 <- mean(stack_logno3)
    mean_logpo4 <- mean(stack_logpo4)
    mean_logsi <- mean(stack_logsi)
    mean_logchla <- mean(stack_logchla)
    mean_logpic <- mean(stack_logpic)
    
    # With na.rm
    
    mean_logchla2 <- mean(stack_logchla, na.rm = T)
    mean_logpic2 <- mean(stack_logpic, na.rm = T)
    
    baseline_logmean <- stack(mean_sst,
                              mean_sss,
                              mean_mld,
                              mean_par,
                              mean_O2,
                              mean_Sistar,
                              mean_Nstar, 
                              mean_logno3,
                              mean_logpo4,
                              mean_logsi,
                              mean_logchla)
    
    names(baseline_logmean) <- c("sst", "sss", "mld", "par", "O2", "Sistar", "Nstar", "logno3", "logpo4", "logsi", "logchla")
    
    raster::writeRaster(baseline_logmean, "data_out/Env_data/baseline_logmean.tif", overwrite = T)
    saveRDS(baseline_logmean,  "data_out/Env_data/baseline_logmean.RDS")
    
    # With na.rm
    
    baseline_logmean.NA <- stack(mean_sst,
                              mean_sss,
                              mean_mld,
                              mean_par2,
                              mean_O2,
                              mean_Sistar,
                              mean_Nstar, 
                              mean_logno3,
                              mean_logpo4,
                              mean_logsi,
                              mean_logchla2)
    
    names(baseline_logmean.NA) <- c("sst", "sss", "mld", "par", "O2", "Sistar", "Nstar", "logno3", "logpo4", "logsi", "logchla")
    
    raster::writeRaster(baseline_logmean.NA, "data_out/Env_data/baseline_logmean.NA.tif", overwrite = T)
    saveRDS(baseline_logmean.NA,  "data_out/Env_data/baseline_logmean.NA.RDS")
    
  
  # We create a monthly scaled log baseline
    
      log_stacks <- c(as(stack_logno3, "SpatRaster"),
                      as(stack_logpo4, "SpatRaster"), 
                      as(stack_logsi, "SpatRaster"),
                      as(stack_logchla, "SpatRaster"), 
                      as(stack_logpic, "SpatRaster")
      )
      
      baseline2 <- c(as(stack_SST, "SpatRaster"), 
                     as(stack_SSS, "SpatRaster"),
                     as(stack_MLD, "SpatRaster"),
                     as(stack_PAR, "SpatRaster"),
                     as(stack_O2, "SpatRaster"),
                     as(stack_Sistar, "SpatRaster"),
                     as(stack_Nstar, "SpatRaster"))
      
      
      baseline_log <- c(baseline2, log_stacks)
      
      saveRDS(baseline_log, "data_out/Env_data/baseline_log.RDS")


  

# ----- 3. Monthly averaged Rasters -----
   
    names_baseline_log <- names(baseline_log)
    
    
    jan_names <- grep("_1$", names_baseline_log, value = T)
    base_1 <- terra::subset(baseline_log, which(names(baseline_log) %in% jan_names))
    
    feb_names <- grep("_2$", names_baseline_log, value = T)
    base_2 <- terra::subset(baseline_log, which(names(baseline_log) %in% feb_names))
    
    mar_names <- grep("_3$", names_baseline_log, value = T)
    base_3 <- terra::subset(baseline_log, which(names(baseline_log) %in% mar_names))
    
    apr_names <- grep("_4$", names_baseline_log, value = T)
    base_4 <- terra::subset(baseline_log, which(names(baseline_log) %in% apr_names))
    
    may_names <- grep("_5$", names_baseline_log, value = T)
    base_5 <- terra::subset(baseline_log, which(names(baseline_log) %in% may_names))
    
    jun_names <- grep("_6$", names_baseline_log, value = T)
    base_6 <- terra::subset(baseline_log, which(names(baseline_log) %in% jun_names))
    
    jul_names <- grep("_7$", names_baseline_log, value = T)
    base_7 <- terra::subset(baseline_log, which(names(baseline_log) %in% jul_names))
    
    aug_names <- grep("_8$", names_baseline_log, value = T)
    base_8 <- terra::subset(baseline_log, which(names(baseline_log) %in% aug_names))
    
    sep_names <- grep("_9$", names_baseline_log, value = T)
    base_9 <- terra::subset(baseline_log, which(names(baseline_log) %in% sep_names))
    
    oct_names <- grep("_10$", names_baseline_log, value = T)
    base_10 <- terra::subset(baseline_log, which(names(baseline_log) %in% oct_names))
    
    nov_names <- grep("_11$", names_baseline_log, value = T)
    base_11 <- terra::subset(baseline_log, which(names(baseline_log) %in% nov_names))
    
    dec_names <- grep("_12$", names_baseline_log, value = T)
    base_12 <- terra::subset(baseline_log, which(names(baseline_log) %in% dec_names))
    
    
    # save it 
    
    if(!dir.exists(paste0("data_out/Env_data/monthly_base")))
    {
      dir.create(paste0("data_out/Env_data/monthly_base"), recursive = T)
    }
    
    #base_names <- grep("^base_", ls(), value = T)
    base_names <- c("base_1", "base_2", "base_3", "base_4", "base_5", "base_6", "base_7", 
                    "base_8", "base_9", "base_10", "base_11", "base_12")
    
    
    for(i in base_names)
    {
      saveRDS(get(i), paste0("data_out/Env_data/monthly_base/", i,".RDS"))
    }
    

    
  
  # ===============================================================================
    
    
    
    
    
    
    
    
    
    
    
