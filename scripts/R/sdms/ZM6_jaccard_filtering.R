

# ----- JACCARD INDEX FILTERING -----

    # This script takes as input data species list and individual models from script ZM5 and compute a jaccard 
    # score for each models, allowing then a model selection based on this metrics. Supplementary plots are 
    # built 

    # Definition
      # Jaccard = TP / (TP + FP + FN)
      
      #                 True Pres  |  True Abs  
      #      ---------+-----------+---------+
      #   Pred Pres  |      TP   |    FP   |
      #     ---------+----------+---------+
      #   Pred Abs  |     FN   |    TN   |
      #   ----------+---------+---------+


    # species list
    sp_tab <- read.csv2("sp_list2.csv")


    # We use the data structure of one arbitrary species to create storage objects 
    sp <- "Heterocapsa_nei_rotundata"
    iter <- 1
    model_runs <- readRDS(paste0("models/", sp, "/model_runs_mod_", sp, ".RDS"))


    # We build an array where all the Jaccard indices will be stored
    
    allsp_jaccard <- my_array <- array(
      dim = c(
        nrow(sp_tab),
        length(unique(model_runs@models.evaluation@val$algo)),
        length(unique(model_runs@models.evaluation@val$run))
      ),
      dimnames = list(
        species = sp_tab$sp,
        model = unique(model_runs@models.evaluation@val$algo),
        cv.run = unique(model_runs@models.evaluation@val$run)
      )
    ) 


    # we build a similar array where the maximum jaccard threshlods are stored in order to build the binary 
    # ensemble maps
    
    allsp_cutoffs <- allsp_jaccard
    
    
    jaccard_list <- list()
    
    
    sp_list_aza <- sort(subset(sp_list, sp_list !="Azadinium_concinnum"))
    
    for (sp in sp_list_aza)
    {
      cat(paste("IIIIIIII", Sys.time(), sp, " evaluation initialised IIIIIIII\n", sep = " "))
      
      model_runs_mod <- readRDS(paste0("models/", sp, "/model_runs_mod_", sp, ".RDS"))
      
      
      input_data <- get(load(model_runs_mod@formated.input.data@link))
      calib_lines <- get(load(model_runs_mod@calib.lines@link)) # array
      model_preds <- get(load(model_runs_mod@models.prediction@link))
      
      for (cv in unique(model_runs_mod@models.evaluation@val$run))
      {
        obs_data <- input_data@data.species[which(!calib_lines[,
                                                               paste0("_", cv), 
                                                               "_allData"])] # we keep calib_lines == F
                                                                             # (Means only validation rows) 
        
        cur_eval <- which(!calib_lines[, 
                                       paste0("_", cv), 
                                       "_allData"])
        
        
        for (md in unique(model_runs_mod@models.evaluation@val$algo))
        {
          # Only the rows that were not used for calibration are evaluated
          cur_preds <- model_preds[which(model_preds$points %in% cur_eval & # Evaluation lines
                                           model_preds$algo == md & # Algorithm
                                           model_preds$run == cv),]
          
          
          
          if(nrow(cur_preds) != 0) #  If there are NAs, the model failed, so we don't compute the Jaccard index
          {
            # Compute the Jaccard index for all thresholds between 0 and 1
            jaccard_test <- NULL
            
            for(cutoff in seq(0, 1000, by = 1)) # Thresholds are created (1000 times here)
            {
              pred_pa <- cur_preds$pred
              pred_pa[pred_pa < cutoff] <- 0 #On the vector of predicted presences/absences, values below the threshold are assigned 0
              pred_pa[pred_pa >= cutoff] <- 1 #Assign 1 to values above the threshold
              TP <- length(which(obs_data == 1 & pred_pa == 1)) # How many times did we get true positives?
              FN <- length(which(obs_data == 1 & pred_pa == 0)) 
              FP <- length(which(obs_data == 0 & pred_pa == 1))
              #TN <- length(which(obs_data == 0 & pred_pa == 0))
              jaccard <- TP / (TP + FP + FN)
              jaccard_test <- rbind.data.frame(jaccard_test,
                                               data.frame(cutoff = cutoff,
                                                          TP = TP,
                                                          FN = FN,
                                                          FP = FP,
                                                          jaccard = jaccard))
            }
            jaccard_list[[sp]][[paste0(cv, "_", md)]] <- jaccard_test
            
            #  The threshold with the highest Jaccard index is selected
            #  If the best threshold is found multiple times, take their average
            allsp_cutoffs[sp, md, cv] <- mean(jaccard_test$cutoff[
              which(jaccard_test$jaccard == max(jaccard_test$jaccard))])
            
            # The Jaccard index is extracted at the mean best threshold
            allsp_jaccard[sp, md, cv] <- jaccard_test$jaccard[
              which(jaccard_test$cutoff == round(allsp_cutoffs[sp, md, cv]))]
          } else 
          {   # -->  "If" statement above; below, the array is filled with NAs
            
            jaccard_list[[sp]][[paste0(cv, "_", md)]] <- NA
            
            allsp_cutoffs[sp, md, cv] <- NA 
            
            allsp_jaccard[sp, md, cv] <- NA
          }
        }
      }
    }



      if(!dir.exists(paste0("evaluation/")))
      {
        dir.create(paste0("evaluation/"), recursive = T) 
      }



      # Save the cutoffs for potential reuse in the ensemble model
      saveRDS(allsp_cutoffs, file = "evaluation/jaccard_cutoffs.RDS")
      # And we save the Jaccard evaluations
      saveRDS(allsp_jaccard, file = "evaluation/jaccard_evals.RDS")
      saveRDS(jaccard_list, file = "evaluation/jaccard_tests.RDS")



     #-------- Plots --------- 

      allsp_jaccard <- readRDS("evaluation/jaccard_evals.RDS")
      
      ggjaccard <- reshape2::melt(allsp_jaccard)
      
      saveRDS(ggjaccard, file = "evaluation/ggjaccard.RDS")


      # ----- 1. Jaccard selection outputs ------
      
        # in order to see if the threshold that we chose has an influence with the type of data we select 
        ggjaccard.troph <- left_join(ggjaccard, sp_tab[1:3], by = c("species" = "sp"))
        
        saveRDS(ggjaccard.troph, file = "evaluation/ggjaccard.troph.RDS")
        
        ggplot(subset(ggjaccard, ggjaccard$species == "Azadinium_caudatum"), aes(x = model, y = value)) +
          geom_boxplot() +
          theme_minimal() 
        # +
        #   xlim(0,1) + ylim(0,1) 
        # 
        # + facet_wrap (~ species)
        
        ggplot(ggjaccard, aes(x = model, y = value)) +
          geom_boxplot() +
          theme_minimal()  -> jaccard_all
        
        ggplot(ggjaccard, aes(x = value, y = species)) +
          geom_boxplot() +
          theme_minimal() 
        
        ggjaccard.troph$trophy <- factor(ggjaccard.troph$trophy, levels = c("mixotroph",
                                                                            "phototroph",
                                                                            "phagotroph"))

        ggplot(subset(ggjaccard.troph, ggjaccard$value > 0.3), aes(x = value, y = fct_reorder(species, value, .fun = mean), color = trophy)) +
          geom_boxplot() +
          theme_minimal() +
          ylab("Species") +
          xlab("Jaccard index value") -> jaccard_species

      # general plot 
      
      ggarrange(jaccard_all, jaccard_species, ncol = 2, widths= c(1, 3)) -> jaccard_suppmat
      ggsave("suppmat/jaccard.pdf", jaccard_suppmat, width = 15, height = 8)




      # ------ 2. Response curves per trophic type, Jaccard filtered --------


        if(!dir.exists(paste0("resp_curves")))
        {
          dir.create(paste0("resp_curves"), recursive = T) 
        }



        # empty dataframe to store response data 
        resp_full <- data.frame(
          Index = numeric(),
          Variable = factor(),
          Var.Value = numeric(), 
          Model = factor(),
          Response = numeric(), 
          Trophy = character()
        )


        for (sp in sp_list_aza) # sp_list_photo, or sp_list_mixo 
        {
          
          mod_data <- readRDS(paste0("models/", sp, "/model_runs_mod_", sp, ".RDS"))
          input_data <- readRDS(paste0("models/", sp, "/run_data", sp, ".RDS"))
          
          # Variables used for calibration 
          cur_vars <- mod_data@expl.var.names
          
          #  Response curves computing 
          resp <- bm_PlotResponseCurves(bm.out = mod_data,
                                        fixed.var = "mean",
                                        data_species = input_data@data.species,#argument caché
          )$tab
          
          
          colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")
          
          # we order variables 
          resp$Variable <- factor(resp$Variable, levels = cur_vars)
          resp$Trophy <- sp_tab[which(sp_tab$sp == sp),]$trophy
          
          resp_curve  <- ggplot(resp, aes(x = Var.value, y = Response))+
            geom_line(alpha = 0.2, aes(group = Model)) +
            stat_smooth() +
            facet_wrap(~Variable, scales = "free_x") +
            theme_bw() +
            ylim(0, 1.1) +
            xlab("Variable value") +
            ggtitle(paste0(sp))
          
          
          ggsave(paste0("resp_curves/", sp, ".jpg"), plot = resp_curve, dpi = 300)
          dev.off()
          
          resp_full <- rbind(resp_full, resp)
          
          
        }


         resp_full$algo <- sub(".+_(\\w+)$", "\\1", resp_full$Model)
         
         # We change expl.var for SST
         resp_full$Variable <- gsub("expl.var", "SST", resp_full$Variable)

       
      # --> Remove 0.3 Jscores from curves
        
        # We spot the model names that have > 0.3 Jscore
        
        ggjaccard.troph03 <- ggjaccard.troph
        
        ggjaccard.troph03$mod_name <- paste0(gsub("_", ".", ggjaccard.troph03$species), 
                                             "_allData_",
                                             ggjaccard.troph03$cv.run,
                                             "_", ggjaccard.troph03$model)
        
        ggjaccard.troph03 <- subset(ggjaccard.troph03, value > 0.3)
        
        resp_full3 <- resp_full2
        
        resp_full3 <- resp_full3[which(resp_full3$Model %in% ggjaccard.troph03$mod_name),]
        
        
        resp_mixo03 <- resp_full3[which(resp_full3$Trophy == "mixotroph"),]
        
        ggplot(resp_mixo03, aes(x = Var.value, y = Response, color = algo))+
          geom_line(alpha = 0.05, aes(group = Model)) +
          stat_smooth() +
          facet_wrap(~Variable, scales = "free_x") +
          theme_bw() +
          ylim(0, 1.1) +
          xlab("Variable value") +
          ggtitle("Mixotrophic species Response Curves") +
          theme(strip.background = element_rect(fill = "white", color = "black") , 
                strip.text = element_text(color = "firebrick1", face = "bold", size = 10)) -> mixo_curves03
        
        ggsave("suppmat/mixo_curve.pdf", mixo_curves03, width = 10, height = 7)
        
        resp_photo03 <- resp_full3[which(resp_full3$Trophy == "phototroph"),]
        
        ggplot(resp_photo03, aes(x = Var.value, y = Response, color = algo))+
          geom_line(alpha = 0.05, aes(group = Model)) +
          stat_smooth() +
          facet_wrap(~Variable, scales = "free_x") +
          theme_bw() +
          ylim(0, 1.1) +
          xlab("Variable value") +
          ggtitle("Phototrophic species Response Curves") +
          theme(strip.background = element_rect(fill = "white", color = "black") , 
                strip.text = element_text(color = "olivedrab3", face = "bold", size = 10)) -> photo_curves03
        
        ggsave("suppmat/photo_curve.pdf", photo_curves03, width = 10, height = 7)
        
        resp_phago03 <- resp_full3[which(resp_full3$Trophy == "phagotroph"),]
        
        ggplot(resp_phago03, aes(x = Var.value, y = Response, color = algo))+
          geom_line(alpha = 0.05, aes(group = Model)) +
          stat_smooth() +
          facet_wrap(~Variable, scales = "free_x") +
          theme_bw() +
          ylim(0, 1.1) +
          xlab("Variable value") +
          ggtitle("Phagotrophic species Response Curves") +
          theme(strip.background = element_rect(fill = "white", color = "black") , 
                strip.text = element_text(color = "dodgerblue", face = "bold", size = 10)) -> phago_curves03
        
        ggsave("suppmat/phago_curve.pdf", phago_curves03, width = 10, height = 7)


# --------------------------------------------------------------------------------











  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  