
# ----- Monthly Species occurrences Maps ------

  # This script takes as input monthly-separated abundance matrix as well as sites datasets built in script Z3 for 
  # generating monthly-separated occurences (presence AND absence) datasets for each species before to rasterize it
  
  library(tidyverse)

  # ---- 1. Presence ------

    if(!dir.exists(paste0("month_occ/sp_sites/")))
    {
      dir.create(paste0("month_occ/sp_sites/"), recursive = T)
    }
    
    if(!dir.exists(paste0("month_occ/sp_maps/")) )
    {
      dir.create(paste0("month_occ/sp_maps/"), recursive = T)
    }
    
    
    for (m in 1:12)
    {
      cat(paste("||||||||||", Sys.time(), m, "||||||||||\n", sep = " "))
      
      if(!dir.exists(paste0("month_occ/sp_sites/", m, "/")))
      {
        dir.create(paste0("month_occ/sp_sites/", m, "/"), recursive = T)
      }
      
      if(!dir.exists(paste0("month_occ/sp_maps/", m, "/")))
      {
        dir.create(paste0("month_occ/sp_maps/", m, "/"), recursive = T)
      }
    
        mothermat <- readRDS(paste0("month_occ/matrix/mothermat_", m, ".RDS")) # monthly resolved abundance matrix
        site.i <-  readRDS(paste0("month_occ/sites/month_occ_", m, ".RDS")) # monthly resolved sites metadata
        
        row.names(mothermat)[row.names(mothermat) == "Heterocapsa_nei/rotundata"] <- "Heterocapsa_nei_rotundata"
        
        for(i in 1:length(rownames(mothermat)))
        {
          
          sp <- rownames(mothermat)[i]
          cat(paste(".........", Sys.time(), sp, "........\n", sep = " "))
          
          #row <- as.data.frame(t(mothermat[sp,][1:nrow(mothermat)] > 0))
          row <- as.data.frame(t(mothermat[sp,] > 0))
          colnames(row) <- "bool"
          row$site <- rownames(row)
          which_row <- which(row$bool == TRUE)
          which_site <-  row[which_row,]$site
            
          #coord_sp <- site.i[which(which_site %in% site.i$site_ID),]
          coord_sp <- site.i[which(site.i$site_ID %in% which_site),]
          coord_sp$ap <- 1
          
         
          saveRDS(coord_sp, paste0("month_occ/sp_sites/", m, "/pres_", sp, ".RDS"))
            
          
          
        }# end for loop
      }
        
     
    # ----- 2. Absence -------
        
    
    for (m in 1:12)
    {
      
        if(!dir.exists(paste0("month_occ/sp_sites_abs/", m, "/")))
        {
          dir.create(paste0("month_occ/sp_sites_abs/", m, "/"), recursive = T)
        }
        
        mothermat <- readRDS(paste0("month_occ/matrix/mothermat_", m, ".RDS")) 
        site.i <-  readRDS(paste0("month_occ/sites/month_occ_", m, ".RDS")) 
        

        
          for(i in 1:length(rownames(mothermat)))
          {
            sp <- rownames(mothermat)[i]
            
            
            #row <- as.data.frame(t(mothermat[sp,][1:nrow(mothermat)] == 0))
            row <- as.data.frame(t(mothermat[sp,] == 0))
            colnames(row) <- "bool"
            row$site <- rownames(row)
            not_which_row <- which(row$bool == TRUE)
            not_which_site <-  row[not_which_row,]$site
            
            #coord_sp_abs <- site.i[which(not_which_site %in% site.i$site_ID),]
            coord_sp_abs <- site.i[which(site.i$site_ID %in% not_which_site),]
            coord_sp_abs$ap <- 0
            
            saveRDS(coord_sp_abs, paste0("month_occ/sp_sites_abs/", m, "/abs_", sp, ".RDS"))
        }
    }
    

  # ================================================================================================









    






