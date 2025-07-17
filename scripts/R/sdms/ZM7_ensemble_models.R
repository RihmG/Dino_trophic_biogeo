
# ====== Ensemble models - projection ======== 

# This script takes as input data the individual models and the jaccard matrix associated and build global 
# projections, as well as ensemble models based on the average and the standard-deviation of each individual maps
  
    library(tidyverse)
    library(tidyr)  
    library(terra)
    library(biomod2)
    library(viridis)
    library(corrplot)
    library(reshape2)
    library(plyr)
    library(gridExtra)

  # ----- 0. Data ------

    base_names <- c("base_1", "base_2", "base_3", "base_4", "base_5", "base_6", "base_7", 
                    "base_8", "base_9", "base_10", "base_11", "base_12")
    
    sp_tab <- read.csv2("sp_list.csv")
    
    sp_list <- sp_tab$sp
    # There is unresolved problem with Azadinium_concinuum so we create a sp_list without it 
    sp_list_aza <- sort(subset(sp_list, sp_list !="Azadinium_concinnum"))
    
    # We create sp lists depending on trophic type
    # Alphabetically ranked in order to follow the calibration 
    sp_list_phago <- sort(sp_tab[sp_tab$trophy == "phagotroph ",]$sp) # careful to the space after phagotroph... 
    sp_list_photo <- sort(sp_tab[sp_tab$trophy == "phototroph",]$sp)
    sp_list_photo_aza <- subset(sp_list_photo, sp_list_photo !="Azadinium_concinnum")
    sp_list_mixo <- sort(sp_tab[sp_tab$trophy == "mixotroph",]$sp)


    # Create an appropriate directory 
    
    if(!dir.exists(paste0("projection/")))
    {
      dir.create(paste0("projection/"), recursive = T) 
    }
    
    # Jaccard values
    ggjaccard.troph <- readRDS("evaluation/ggjaccard.troph.RDS")


  # ---- 1. Projection of all the individual models-----

    for (sp in sp_list_aza) # sp_list_photo, or sp_list_mixo 
    {
      cat(paste("IIIIIIII", Sys.time(), sp, " Projection initialised IIIIIIII\n", sep = " "))
      
      for (base in base_names)
      {
        cat(paste(".................", Sys.time(), base, ".................\n", sep = " "))
        
        iter <- gsub("base_", "", base) 
        
        
        # get the data 
        sel_vars <- readRDS(paste0("var_selection/sel_vars/sel_vars_f", sp, ".RDS"))
        
        mod2proj <- readRDS(paste0("models/", sp, "/model_runs_mod_", sp, ".RDS"))
        
        env_data <- readRDS(paste0("data_out/Env_data/monthly_base/", base, ".RDS"))
        
        #mod2proj@expl.var.names <- sel_vars[[1]]
        
        # transform the data 
        namenv <- names(env_data)
        namenv <- gsub(paste0("_", iter, "$"), "", namenv)
        names(env_data) <- namenv
        
        calib_env_data <- env_data[[sel_vars[[1]]]] # We fit the selected variables with baseline raster 
        
        if (length(sel_vars[[1]]) == 1)
        {
          names(calib_env_data) <- "expl.var"
        }
        
        projection_runs <- BIOMOD_Projection(bm.mod = mod2proj,
                                             proj.name = iter,
                                             new.env = calib_env_data, 
                                             models.chosen = "all", 
                                             binary.meth = NULL,
                                             #metric.binary = "TSS", 
                                             metric.filter = NULL, 
                                             build.clamping.mask = TRUE, 
                                             nb.cpu = 4)
        
        if(!dir.exists(paste0("projection/", sp)))
        {
          dir.create(paste0("projection/", sp), recursive = T) 
        }
        
        saveRDS(projection_runs, file = paste0("projection/", sp, "/proj_", iter, ".RDS"))
        
        
      }
    }





  # ---- 2. Selection of best individual models ------


    # After all the individual models have been projected, we select them depending on their Jaccard performances
    # For each species we build a stack of projections, after Jaccard validation
    # We build a df that will store the number of models that are kept after Jaccard selection (for performance statistics)

    # Removing GBM models 
      ggjaccard.troph <- subset(ggjaccard.troph, model != "GBM" )
      
      
      models_kept <- data.frame(species = character(0),
                                models_built = numeric(0),
                                models_kept = numeric(0),
                                ratio = numeric(0))
      
      sp_rank <- 1 # incrementation 
      sp_list_jaccard <- c() # List that store species with individual models selected 

      
      for (sp in sp_list_aza) # sp_list_photo, or sp_list_mixo 
      {
        cat(paste("IIIIIIII", Sys.time(), sp, " Selection initialised IIIIIIII\n", sep = " "))
        
        mod2proj <- readRDS(paste0("models/", sp, "/model_runs_mod_", sp, ".RDS"))
        
        # we look for jaccard max scores per indiv.models 
        jacc_df <- ggjaccard.troph[ggjaccard.troph$species == sp,]
        jacc_df<- na.omit(jacc_df[,1:5])
        threshold <- 0.3
        jacc_fil <- jacc_df[jacc_df$value > threshold,]
        
        if (nrow(jacc_fil) > 0)
        {
          
          # Extract the names of individual models 
          indiv_list <- c()
          for (k in 1:(nrow(jacc_fil)))
          {
            patt <- paste0(jacc_fil$cv.run[k], "_", jacc_fil$model[k])
            models2extract <- get_built_models(mod2proj)[grep(patt, get_built_models(mod2proj))]
            indiv_list <- c(indiv_list, models2extract)
            
          }
          
          # Projection with monthly env. data 
          for (base in base_names)
          {
            iter <- gsub("base_", "", base)
            
            ProjMods <- readRDS(paste0("projection/", sp, "/proj_", iter, ".RDS"))
            
            # extract raster layers
            projstack <- unwrap(ProjMods@proj.out@val)
            
            #Stack with filtered models 
            ProjMods_fil <- projstack[[indiv_list]]
            
            # Computation 
            # Mean
            mean_projmod <- mean(ProjMods_fil)
            
            # standard deviation of the individual models 
            sd_projmod <- app(ProjMods_fil, sd)
            
            # Now we penalise the mean with the SD
            mean_pen_mod <- mean_projmod - sd_projmod # Penalised probability of presence
            mean_pen_mod[mean_pen_mod < 0] <- 0
            
            if(!dir.exists(paste0("projection/", sp, "/base_", iter, "/ensemble_fil")))
            {
              dir.create(paste0("projection/", sp, "/base_", iter, "/ensemble_fil"), recursive = T) 
            }
            
            
            saveRDS(ProjMods_fil, file = paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/stack_fil", iter, ".RDS"))
            saveRDS(mean_projmod, file = paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/mean_em_fil_", iter, ".RDS"))
            saveRDS(sd_projmod, file = paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/sd_em_fil_", iter, ".RDS"))
            saveRDS(mean_pen_mod, file = paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/mean_em_pen_fil_", iter, ".RDS"))
            
            # Fill the performance stats array 
            
            
          }
          
          sp_list_jaccard <- c(sp_list_jaccard, sp) 
          
        } else { 
          cat("\n", "for ", sp, ", no individual models were retained", "\n")
        }
        
        models_kept[sp_rank,"species"] <- sp
        models_kept[sp_rank,"models_built"] <- nrow(jacc_df)
        models_kept[sp_rank, "models_kept"] <- nrow(jacc_fil)
        models_kept[sp_rank, "ratio"] <- nrow(jacc_fil) / nrow(jacc_df)
        
        
        sp_rank <- sp_rank + 1
        
      }


      saveRDS(sp_list_jaccard, "sp_list_jaccard.RDS")





  # ----- 3. Binary Projections ---------

    # Building of Binary Projection based on Jaccard cutoff 
    
    # We first create a list that takes only species with ensemble models built (i.e kept
    # after Jaccard selection)
    

      jaccard_cutoff <- readRDS(file = "evaluation/jaccard_cutoffs.RDS")
      
      for (sp in sp_list_jaccard) 
      {
        cat(paste("IIIIIIII", Sys.time(), sp, " Selection initialised IIIIIIII\n", sep = " "))
        
        for (base in base_names)
        {
          cat(paste("..........", Sys.time(), sp, base,  " initialised ..........\n", sep = " "))
          
          iter <- gsub("base_", "", base) 
          
          # stack of individual models kept
          projmod_fil <- readRDS(paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/stack_fil", iter, ".RDS"))
          #plot(projmod_fil)
          #stacks of ensemble models made with individual models kept 
          mean_projmod <- readRDS(file = paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/mean_em_fil_", iter, ".RDS"))
          sd_projmod <- readRDS(file = paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/sd_em_fil_", iter, ".RDS"))
          mean_pen_mod <- readRDS(file = paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/mean_em_pen_fil_", iter, ".RDS"))
          
          #biomod object
          mod2proj <- readRDS(paste0("models/", sp, "/model_runs_mod_", sp, ".RDS"))
          
          #Input data 
          input_data <- get(load(mod2proj@formated.input.data@link))
          sp_coords <- input_data@coord
          
          #cutoff data 
          tmp <- reshape2::melt(jaccard_cutoff[sp,,])
          cur_cutoffs <- tmp$value
          sp.name <- gsub("_", ".", sp)
          names(cur_cutoffs) <- paste(sp.name, "allData", tmp$cv.run, tmp$model, sep = "_")
          
          # Building of binary stacks 
          
          # Absence/presence by jaccard cutoff 
          binary_stack <- bm_BinaryTransformation(projmod_fil, cur_cutoffs[names(projmod_fil)]) #on convertit les baselines en binaire 
          saveRDS(binary_stack, paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/stack_bin_", iter, ".RDS"))
          
          # create a plot with all the individual models binary projections 
          png(paste0("projection/", sp, "/base_", iter, "/indiv_bin", iter, ".png"), width = 2000, height = 1500)
          plot(binary_stack, ncol = 5, nrow = 3, cex.main = 2.5)
          dev.off()
          
          
          obs_data <- input_data@data.species #1 and 0
          
          pred_data <- terra::extract(mean_projmod, sp_coords) #on extrait la valeur du modèle d'ensemble aux coordonées 
          jaccard_test <- NULL
          
          for(cutoff in seq(0, 1000, by = 1))
          {
            pred.pa <- pred_data$mean
            pred.pa[pred.pa < cutoff] <- 0
            pred.pa[pred.pa >= cutoff] <- 1
            TP <- length(which(obs_data == 1 & pred.pa == 1))
            FN <- length(which(obs_data == 1 & pred.pa == 0))
            FP <- length(which(obs_data == 0 & pred.pa == 1))
            jaccard <- TP / (TP + FP + FN)
            jaccard_test <- rbind.data.frame(jaccard_test,
                                             data.frame(cutoff = cutoff,
                                                        TP = TP,
                                                        FN = FN,
                                                        FP = FP,
                                                        jaccard = jaccard))
          }
          em_cutoff <- mean(jaccard_test$cutoff[which(jaccard_test$jaccard == max(jaccard_test$jaccard))])
          
          # Ensemble model transformé en présence-absence
          mean_em_binary <- bm_BinaryTransformation(mean_projmod, em_cutoff)
          # Transformation of probabilities into presence-absence, penalized by uncertainty
          mean_em_binary_penalised <- bm_BinaryTransformation(mean_pen_mod, em_cutoff)
          
          writeRaster(mean_em_binary, paste0("projection/", sp, "/base_", iter, "/mean_bin", sp, "_", iter, ".tif"), overwrite = T)
          writeRaster(mean_em_binary_penalised, paste0("projection/", sp, "/base_", iter, "/meanpen_bin", sp,"_", iter, ".tif"), overwrite = T)
          
          
        } 
        
      }

      # We check how much individual model we lose depending on the chosen Jaccard cutoff
      
      hist(ggjaccard.troph$value)
      
      ggjaccard.troph.na <- ggjaccard.troph[which(!is.na(ggjaccard.troph$value)),]
      
      
      evol_cut <- data.frame(
        cutoff = rep(NA, 100),
        n_models = rep(NA, 100)
      )
      iter <- 1
      for (cut in seq(0, 1, by = 0.01))
      {
        n_mod <- nrow(subset(ggjaccard.troph.na, value > cut))
        evol_cut[iter, "cutoff"] <- cut
        evol_cut[iter, "n_models"] <- n_mod
        iter <- iter + 1
      }
      
      ggplot(evol_cut, aes(x = cutoff, y = n_models)) + 
        geom_point() +
        theme_bw()
      
     

  #  ----- 4. Algorithms-based projections --------
      
      # Projection based on the modeling method (supplementary figures)

      jaccard_cutoff <- readRDS(file = "evaluation/jaccard_cutoffs.RDS")
      
      for (sp in sp_list_jaccard) 
      {
        cat(paste("IIIIIIII", Sys.time(), sp, " Selection initialised IIIIIIII\n", sep = " "))
        
        for (base in base_names)
        {
          cat(paste("..........", Sys.time(), sp, base,  " initialised ..........\n", sep = " "))
          
          iter <- gsub("base_", "", base) 
          
          # stack of individual models kept
          projmod_fil <- readRDS(paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/stack_fil", iter, ".RDS"))
          #plot(projmod_fil)
          #stacks of ensemble models made with individual models kept 
          mean_projmod <- readRDS(file = paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/mean_em_fil_", iter, ".RDS"))
          sd_projmod <- readRDS(file = paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/sd_em_fil_", iter, ".RDS"))
          mean_pen_mod <- readRDS(file = paste0("projection/", sp, "/base_", iter, "/ensemble_fil", "/mean_em_pen_fil_", iter, ".RDS"))
          
          #biomod object
          mod2proj <- readRDS(paste0("models/", sp, "/model_runs_mod_", sp, ".RDS"))
          
          #Input data 
          input_data <- get(load(mod2proj@formated.input.data@link))
          sp_coords <- input_data@coord
          
          #cutoff data 
          tmp <- reshape2::melt(jaccard_cutoff[sp,,])
          cur_cutoffs <- tmp$value
          sp.name <- gsub("_", ".", sp)
          names(cur_cutoffs) <- paste(sp.name, "allData", tmp$cv.run, tmp$model, sep = "_")
          
          # Building of binary stacks 
          
          # Absence/presence by jaccard cutoff 
          binary_stack <- bm_BinaryTransformation(projmod_fil, cur_cutoffs[names(projmod_fil)]) #on convertit les baselines en binaire 
          
          # Now we extract the layers depending on the algorithm used
          if(!dir.exists(paste0("projection/", sp, "/base_", iter, "/algo_proj/")))
          {
            dir.create(paste0("projection/", sp, "/base_", iter, "/algo_proj/"), recursive = T) 
          }
          
          
          GLM_layers <- grep("_GLM$", names(binary_stack), value = TRUE)
          GAM_layers <- grep("_GAM$", names(binary_stack), value = TRUE)
          ANN_layers <- grep("_ANN$", names(binary_stack), value = TRUE)
          MARS_layers <- grep("_MARS$", names(binary_stack), value = TRUE)
          
          # Extract the layers only if there is corrsponding algorithm within rasterstack
          
          
          if(length(GLM_layers) > 0)
          {
            binary_GLM <- binary_stack[[GLM_layers]]
            mean_GLM <- mean(binary_GLM)
            saveRDS(mean_GLM, paste0("projection/", sp, "/base_", iter, "/algo_proj/GLM_mean_algox", iter, ".RDS"))
          } else {
          }
          
          if(length(GAM_layers) > 0)
          {
            binary_GAM <- binary_stack[[GAM_layers]]
            mean_GAM <- mean(binary_GAM)
            saveRDS(mean_GAM, paste0("projection/", sp, "/base_", iter, "/algo_proj/GAM_mean_algox", iter, ".RDS"))
          } else {
          }
          
          if(length(ANN_layers) > 0)
          {
            binary_ANN <- binary_stack[[ANN_layers]]
            mean_ANN <- mean(binary_ANN)
            saveRDS(mean_ANN, paste0("projection/", sp, "/base_", iter, "/algo_proj/ANN_mean_algox", iter, ".RDS"))
          } else {
          }
          
          if(length(MARS_layers) > 0)
          {
            binary_MARS <- binary_stack[[MARS_layers]]
            mean_MARS <- mean(binary_MARS)
            saveRDS(mean_MARS, paste0("projection/", sp, "/base_", iter, "/algo_proj/MARS_mean_algox", iter, ".RDS"))
          } else {
          }
          
        }
        
      }



# ======================================================================================================




