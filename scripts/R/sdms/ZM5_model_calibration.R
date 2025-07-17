

  # ----- MODELS CALIBRATION -------

  # This script takes as input the species presence/absence dataframes built in ZM3 and ZM4 and creates 
  # individual models 

        library(tidyverse)
        library(tidyr)  
        library(terra)
        library(biomod2)
        library(virtualspecies)
        library(Rarity)
        library(viridis)
        library(corrplot)
        library(reshape2)
        library(plyr)
        library(gridExtra)

  # ------- 1. Data importation ---------

      # SPECIES DATA 
      
        sp_tab <- read.csv2("sp_list2.csv")
        
        sp_list <- sp_tab$sp
        # There is unresolved problem with Azadinium_concinuum so we create a sp_list without it 
        sp_list_aza <- subset(sp_list, sp_list !="Azadinium_concinnum")
        
        # We create sp lists depending on trophic type
        # Alphabetically ranked in order to follow the calibration 
        sp_list_phago <- sort(sp_tab[sp_tab$trophy == "phagotroph",]$sp) # careful to the space after phagotroph... 
        sp_list_photo <- sort(sp_tab[sp_tab$trophy == "phototroph",]$sp)
        sp_list_photo_aza <- subset(sp_list_photo, sp_list_photo !="Azadinium_concinnum")
        sp_list_mixo <- sort(sp_tab[sp_tab$trophy == "mixotroph",]$sp)
      
      
      # ENV DATA
  
        # Monthly env. data
        baseline <- rast("data_out /Env_data/baseline.tif")
        baseline_log <- readRDS("data_out/Env_data/baseline_log.RDS") # with log nutrients
  
        base_files <- list.files("data_out/Env_data/monthly_base/")
  
        base_names <- c("base_1", "base_2", "base_3", "base_4", "base_5", "base_6", "base_7",
                        "base_8", "base_9", "base_10", "base_11", "base_12")
  
        for(i in base_files)
        {
          base_object <- readRDS(paste0("data_out/Env_data/monthly_base/",i))
          assign(gsub(".RDS", "", i), base_object)
        }
        


  # ------- 2. Models formationg and calibration  ---------- 

      algos = c('GLM',
                "GAM",
                "ANN",
                "MARS",
                #"GBM"
                #,"FDA",
                #"RF",
                #"MAXNET"
      )
  
      # Repositery creation 
      
        
      for (sp in sp_list_aza) # start of the loop
      {
        cat(paste("IIIIIIII", Sys.time(), sp, " evaluation initialised IIIIIIII\n", sep = " "))
        
        if(!dir.exists(paste0("models/", sp)))
        {
          dir.create(paste0("models/", sp), recursive = T) 
        }
    
        sp.ff <- readRDS(paste0("var_selection/sp_df/", sp, ".RDS"))
        
        sel_vars <- readRDS(paste0("var_selection/sel_vars/sel_vars_f", sp, ".RDS"))
        
    
     # Model Options 
        
        cat(paste("-------", Sys.time(), sp, base, " evaluation initialised -------\n", sep = " "))
        
        mod_options <- BIOMOD_ModelingOptions(
          GLM = list( type = 'quadratic',
                      
                      interaction.level = 0,
                      
                      myFormula = NULL,
                      
                      test = 'AIC',
                      
                      family = binomial("logit"),
                      
                      mustart = 0.5,
                      
                      control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),
          
          GAM = list(algo = 'GAM_mgcv',
                     type = 's_smoother',
                     k = 5,
                     interaction.level = 0,
                     myFormula = NULL,
                     family = binomial("logit"),
                     method = 'GCV.Cp',
                     optimizer = c('outer','newton'),
                     select = FALSE,
                     knots = NULL,
                     paraPen = NULL,
                     control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
                                    , maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
                                    , rank.tol = 1.49011611938477e-08
                                    , nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
                                    , optim = list(factr=1e+07)
                                    , newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
                                    , outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE) ), 
          
          GBM = list( distribution = 'bernoulli',
                      n.trees = 100,
                      
                      interaction.depth = 1,
                      
                      n.minobsinnode = 10,
                      
                      shrinkage = 0.001,
                      
                      bag.fraction = 0.5,
                      
                      train.fraction = 1,
                      
                      cv.folds = 2,
                      
                      keep.data = FALSE,
                      
                      verbose = FALSE,
                      
                      perf.method = "cv") 
        )
        
    
        # Calibration datasets 
        calib_env_data <- sp.ff[,sel_vars[[1]]]
        coorxy <- sp.ff[, c("x", "y")]
        presabs <- sp.ff[, "ap"]
        
        # Model formating 
        
        run_data_mod <- BIOMOD_FormatingData(resp.name = sp,
                                             resp.var = presabs, 
                                             expl.var = calib_env_data, 
                                             dir.name = "models",
                                             resp.xy = coorxy) 
        
        saveRDS(run_data_mod, file = paste0("models/", sp, "/run_data", sp, ".RDS"))
        
        
        # Model calibration 
        
        
        model_runs_mod <- BIOMOD_Modeling(bm.format = run_data_mod, 
                                          modeling.id = '1', 
                                          models = algos,
                                          bm.options = mod_options,
                                          nb.rep = 5, # CV runs 
                                          data.split.perc = 70, # % data for calibration 
                                          do.full.models = FALSE, 
                                          weights = NULL, 
                                          prevalence = 0.7, # prevalence, if > 0.5, presence have more weight than absences
                                          var.import = 4, # Randomization Runs for variable importance computing 
                                          metric.eval = c("TSS", "ROC"), # Eval metrics 
                                          nb.cpu = 4, # parallelization 
                                          do.progress = TRUE)
        
        saveRDS(model_runs_mod, file = paste0("models/", sp, "/model_runs_mod_", sp, ".RDS"))
    
        # get_variables_importance(model_runs_mod)
        # variable_importance <- get(load(model_runs_mod@variables.importance@link))
        # 
        
        
        # VARIABLE IMPORTANCE
          gg_varimp.m <- get_variables_importance(model_runs_mod)
          
          colnames(gg_varimp.m) <- c("id", 
                                     "PA.Run",
                                     "CV.Run", "Model", 
                                     "Variable", 
                                     "VI.run",
                                     "Variable.importance")
          
          
          gg_varimp.m$Variable <- reorder(gg_varimp.m$Variable,  
                                          gg_varimp.m$Variable.importance,
                                          median, 
                                          na.rm=TRUE)
          
          
          if(!dir.exists(paste0("models/", sp, "/var_import")))
          {
            dir.create(paste0("models/", sp, "/var_import"), recursive = T) 
          }
    
  
          # save plots
            varimp.m <- ggplot(gg_varimp.m, aes(y = Variable, x = Variable.importance)) +
            geom_boxplot() + geom_jitter(alpha = .5, aes(col = Model))  +  theme_classic() + ggtitle(paste0(sp))
          plot(varimp.m)
          ggsave(paste0("models/", sp, "/var_import/var_import_", sp, ".jpg"), plot = varimp.m, dpi = 300)
          dev.off()
      
          # synthetic dataframe
          gg_varimp.m$Variable <- factor(gg_varimp.m$Variable, levels = rev(levels(gg_varimp.m$Variable)))
          quantile.50.m <- aggregate(Variable.importance ~ Variable, data = gg_varimp.m, FUN = median)
          
          
          if(!dir.exists(paste0("models/", sp, "/var_import/var_import_tab")))
          {
            dir.create(paste0("models/", sp, "/var_import/var_import_tab"), recursive = T) 
          }
          
          saveRDS(quantile.50.m, file = paste0("models/", sp, "/var_import/var_import_tab/var_import_tab", sp, ".RDS"))
      
        
        # RESPONSE CURVES
        
          # formatted data and calibrated models loading 
          input_data <- get(load(model_runs_mod@formated.input.data@link))
          
          # Variables used for calibration 
          cur_vars <- model_runs_mod@expl.var.names
          
          #  Response curves computing 
          resp <- bm_PlotResponseCurves(bm.out = model_runs_mod,
                                        fixed.var = "mean",
                                        data_species = input_data@data.species,#argument caché
          )$tab
          
          
          colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")
          
          # we order variables 
          resp$Variable <- factor(resp$Variable, levels = cur_vars)
          
          if(!dir.exists(paste0("models/", sp, "/resp_curves")))
          {
            dir.create(paste0("models/", sp, "/resp_curves"), recursive = T) 
          }
          
          
          # plot
          resp_curve  <- ggplot(resp, aes(x = Var.value, y = Response))+
            geom_line(alpha = 0.2, aes(group = Model)) +
            stat_smooth() +
            facet_wrap(~Variable, scales = "free_x") +
            theme_bw() +
            ylim(0, 1.1) +
            xlab("Variable value") +
            ggtitle(paste0(sp))
          
      
      
            # save it 
            ggsave(paste0("models/", sp, "/resp_curves/", sp, ".jpg"), plot = resp_curve, dpi = 300)
            dev.off()
    
          
        # MODELS PERFORMANCES
    
          if(!dir.exists(paste0("models/", sp, "/models_performances")))
          {
            dir.create(paste0("models/", sp, "/models_performances"), recursive = T) 
          }
          
          eval_scores<-get_evaluations(model_runs_mod)
          
          bm_eval_plot <- bm_PlotEvalMean(bm.out =  model_runs_mod)
          
          eval_plot <- ggplot(data = bm_eval_plot$tab, aes(x = mean1, y = mean2, color = name)) +
            geom_point() +
            labs(x = "ROC", y = "TSS") +
            theme_minimal() +
            xlim(0,1) +
            ylim(0,1) +
            geom_linerange(data = eval_scores$tab, aes(ymin = mean2 - sd2, ymax = mean2 + sd2),
                           width = 0.2, position = position_dodge(0.5)) +
            geom_linerange(data = eval_scores$tab, aes(xmin = mean1 - sd1, xmax = mean1 + sd1),
                           width = 0.2, position = position_dodge(0.5)) 
          
          ggsave(paste0("models/", sp, "/models_performances/", sp, ".jpg"), plot = eval_plot, dpi = 300)
          dev.off()
          
          
        } # end of the loop
  
  
  
        
      
    # ------ 3. Supplementary validations -------
      
        # Synthetic plot of variable importance depending on trophic types 
  
            var_imp_files <- list.files(path = paste0("models/"), pattern = "var_import_tab", recursive = TRUE, full.names = TRUE)
            
            var_import_all <- data.frame(Variable = character(0),
                                         Variable.importance = numeric(0),
                                         species = character(0))
            
            for (path in var_imp_files)
            {
              cur_file <- readRDS(path)
              name_sp <- unlist(strsplit(path, "/")[1])[3] # get species name
              cur_file$species <- name_sp
              var_import_all <- rbind(var_import_all, cur_file)
            }
  
            # we have to replace the "expl.var" in the variable name with the actual name 
            
            explvar <- var_import_all[which(var_import_all$Variable == "expl.var"),]$species # list of species with unique variable 
            
            
            for (sp in explvar)
            {
              var.df <- readRDS(paste0("var_selection/sel_vars/sel_vars_f", sp, ".RDS"))[[1]]
              var_import_all[var_import_all$species == sp,]$Variable <- var.df
            }
            
            
            var_import_all <- left_join(var_import_all, sp_tab[,1:2], by = c("species" = "sp"))
            
            var_import_all$trophy <- factor(var_import_all$trophy, levels = c("mixotroph",
                                                                              "phototroph",
                                                                              "phagotroph"))
            
            
            
            ggplot(var_import_all, aes(x = fct_reorder(Variable, Variable.importance, .fun = mean, .desc = T), y = Variable.importance, fill = trophy)) +
              geom_boxplot() + 
              theme_bw() +
              labs(fill = "Trophic type") +
              xlab("Environmental predictor") -> var_imp.plot
            
            
            ggsave("suppmat/varimp.pdf", var_imp.plot, width = 10, height = 6)
  
  
  
  
        # Table of predictors per species (Supp. Mat.) 
  
  
            pred_tab <- data.frame(Species = character(0),
                                   Environmental_predictor =  character(0))
            
            i = 1
            for (sp in sp_list_aza) # sp_list_photo, or sp_list_mixo 
            {
              sel_vars <- readRDS(paste0("var_selection/sel_vars/sel_vars_f", sp, ".RDS"))
              pred_tab[i,1] <- names(sel_vars)
              pred_tab[i,2] <- paste(sel_vars[[1]], collapse = " ")
              i = i + 1
            }
            
             nrow(pred_tab)
            
  
            evenn_troph <- readRDS("evenn_troph.RDS")
            sp_success <- readRDS("sp_success.RDS")
            
            pred_tab.troph <- left_join(pred_tab, evenn_troph[,4:6], by = c("Species" = "species"))
            
            pred_tab.troph[is.na(pred_tab.troph)] <- "-"
            
            colnames(pred_tab.troph) <- c("Species", "Environmental predictor", "Trophic type", "Mixotrophic type")
            
            pred_tab.troph[pred_tab.troph$Species == "Heterocapsa_nei_rotundata",]$`Trophic type` <- "mixotroph"
            
            pred_tab.troph_succ <- left_join(pred_tab.troph, sp_success, by = "Species")
            
            stargazer(as.data.frame(pred_tab.troph_succ),
                      summary = F,
                      rownames = F,
                      title = "Environmental predictors used for model calibration for each species")
            


  #  ==============================================================================



















