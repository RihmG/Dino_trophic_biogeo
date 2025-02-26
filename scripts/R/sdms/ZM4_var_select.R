
# ---------- VARIABLE SELECTION -----------

# Variable Selection script 

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


# -------- 1. Data importation and cleaning ----------

    # SPECIES DATA 

    sp_tab <- read.csv2("sp_list.csv")
    
    
    if(!dir.exists(paste0("data_out/merged_presabs/")))
    {
      dir.create(paste0("data_out/merged_presabs/"), recursive = T)
    }

    # first we import species data for each sp, and we merge presence and absence layers
    
    sp_list <- sp_tab$sp
    # There is unresolved problem with Azadinium_concinuum so we create a sp_list without it 
    sp_list_aza <- subset(sp_list, sp_list !="Azadinium_concinnum")
    
    # We create sp lists depending on trophic type
    # Alphabetically ranked in order to follow the calibration 
    sp_list_phago <- sort(sp_tab[sp_tab$trophy == "phagotroph ",]$sp) # careful to the space after phagotroph... 
    sp_list_photo <- sort(sp_tab[sp_tab$trophy == "phototroph",]$sp)
    sp_list_mixo <- sort(sp_tab[sp_tab$trophy == "mixotroph",]$sp)
    
    
    # sp_crash <- sp_list[52:length(sp_list)]
    
        # ENV DATA
    
    # Monthly env. data 
    base_files <- list.files("data_out/Env_data/monthly_base/")
    
    for(i in base_files)
    {
      base_object <- readRDS(paste0("data_out/Env_data/monthly_base/",i))
      assign(gsub(".RDS", "", i), base_object)
    }
    

# -------- Variables selection ---------

    # The first step is the same for any species, it only depends on the environmental data 
    # The numbers correspond to the months 
    
    #base_lay <- gsub("\\.RDS$", "", base_files) # This is the list of monthly environmental stack to browse
    base_names <- c("base_1", "base_2", "base_3", "base_4", "base_5", "base_6", "base_7", 
                    "base_8", "base_9", "base_10", "base_11", "base_12")
    
    for (base in base_names)
      {
          cat(paste(base, "\n"))
      
          if(!dir.exists(paste0("var_selection/all_vars/interco_var")))
          {
            dir.create(paste0("var_selection/all_vars/interco_var"), recursive = T)
          }
      
      # Intercorrelated Variables 
    
      pdf(paste0("var_selection/all_vars/interco_var/", base, ".pdf"))
      groups_mean <- removeCollinearity(raster::stack(get(base)), plot = T,
                                        multicollinearity.cutoff = 0.5,
                                        method = "spearman")
      
      dev.off()
      saveRDS(groups_mean, paste0("var_selection/all_vars/", base, ".RDS"))
      # Correlation table 
          if(!dir.exists(paste0("var_selection/all_vars/corrplot")))
          {
            dir.create(paste0("var_selection/all_vars/corrplot"), recursive = T)
          }
      
      cor.tab <- na.omit(as.data.frame(get(base)))
      
      pdf((paste0("var_selection/all_vars/corrplot/", base, ".pdf")))
      corrplot(cor(cor.tab, method = "spearman"), addCoef.col = 'black')
      
      dev.off()
      
    }
    
    
    if(!dir.exists(paste0("var_selection/dendro/")))
    {
      dir.create(paste0("var_selection/dendro/"), recursive = T)
    }
    
    
    for (sp in sp_list_aza)
    {
      sp.df <- readRDS(paste0("month_occ/species_env_df/", sp, ".RDS"))
      cordata <- sp.df[,6:length(sp.df)]
      
      # cormat <- cor(cordata, method = "spearman", use = 'complete.obs')
      # bool.df <- as.data.frame(cormat > 0.7 & cormat < 1 | cormat < -0.7)
      # 
      # paires_col <- as.data.frame(which(bool.df == TRUE, arr.ind = TRUE))
      # 
      # paires_col$var1 <- sub("\\..*$", "", rownames(paires_col))
      # paires_col$var2 <- colnames(bool.df)
      
      cormat <- 1 - abs(cor(cordata, method = "spearman", use = 'complete.obs'))
      dist_mat <- stats::as.dist(cormat)
      clust <- stats::hclust(dist_mat, method = 'complete')
      
      pdf(paste0("var_selection/dendro/", sp, ".pdf"))
      plot(clust, main = 'Dendrogramm')
      rect.hclust(clust, h = 0.3)
      title(main = sp, line = 0)
      dev.off()
      
      # 
      # dist_mat <- as.dist(1 - cormat)
      # 
      # clust <- hclust(dist_mat, method = 'complete')
      # 
      # groups <- as.data.frame(as.matrix(stats::cutree(clust, h = 0.3)))
      # 
      # pdf(paste0("var_selection/dendro/", sp, ".pdf"))
      # plot(clust, main = 'Dendrogramm')
      # rect.hclust(clust, h = 0.3)
      # title(main = sp, line = 0)
      # dev.off()
    
    }
    
    # After inspection of all the dendrogramm, we keep a pool of variables that look not correlated
    # with a threshold of 0.7 (Spearman coeff.)
    
    if(!dir.exists(paste0("var_selection/sp_df/")))
    {
      dir.create(paste0("var_selection/sp_df/"), recursive = T)
    }
    
    for (sp in sp_list_aza)
    {
      sp.df <- readRDS(paste0("month_occ/species_env_df/", sp, ".RDS"))
      sp.df <- sp.df[,c(1:6, 8,9,11,13,14,15,16,17)]
      saveRDS(sp.df, paste0("var_selection/sp_df/", sp,".RDS"))
    }
    
    
    
    # SST, Chla, PO4, N03, Si, Sistar, MLD, PIC, PAR
    
    
    # We now run models in order to chose variables depending on the importance of each predictor for 
    # each species
    
    for (sp in sp_list_aza) 
    {
    
      
      #sp.ff <- readRDS(paste0("data_out/merged_presabs/merged_", sp,".RDS"))
      sp.ff <- readRDS(paste0("var_selection/sp_df/", sp, ".RDS"))
      #sp.ff <- sp.ff[,(1:3)]
        
        # for (base in base_names)
        #   {
          #iter <- gsub("base_", "", base) # We find the month corresponding to the variables dataframe 
          #vars_per_sp <- readRDS(paste0("var_selection/all_vars/sel_vars/sel_vars_", iter, ".RDS"))
          #vars_per_sp[1,1] <- sp
          
          # sp_env_stack <- get(base)[[
          #   vars_per_sp[vars_per_sp$species == sp, 
          #               2:ncol(vars_per_sp)]]] # We match the env.  data with the env. raster in order to get pre-selected env. variables
          # 
        
          coorxy <- sp.ff[, c("x", "y")]
          presabs <- sp.ff[, "ap"]
          env <- sp.ff[, 6:14 ]
          
          if(!dir.exists(paste0("var_selection/run_data/")))
          {
            dir.create(paste0("var_selection/run_data/"), recursive = T)
          }
    
          run_data <- BIOMOD_FormatingData(resp.name = sp,
                                           resp.var = presabs, 
                                           expl.var = env,
                                           dir.name = paste0("var_selection/run_data/"),
                                           resp.xy = coorxy, 
                                           PA.nb.rep = 0, # No need of background points since we have true absences 
                                           PA.nb.absences = 0) 
                                           #PA.strategy = 'random')
        
          saveRDS(run_data, file = paste0("var_selection/run_data/", sp, ".RDS"))
          
          models <- c('GLM', 'MARS', 'GAM', 'GBM', 'ANN')
          runs_CV <- 10
        
        
          model_runs <- BIOMOD_Modeling(bm.format = run_data, 
                                        modeling.id = "1", 
                                        models =  models,
                                        # bm.options = ,
                                        nb.rep = runs_CV, 
                                        data.split.perc = 80, 
                                        do.full.models = FALSE, 
                                        weights = NULL, 
                                        prevalence = 0.7,#  prevalence, if > 0.5, presence have more weight than absences
                                        var.import = 10, 
                                        nb.cpu = 4, 
                                        do.progress = TRUE)
        
          if(!dir.exists(paste0("var_selection/model_runs/")))
          {
            dir.create(paste0("var_selection/model_runs/"), recursive = T)
          }
          
          
          saveRDS(model_runs, file = paste0("var_selection/model_runs/model_runs_",sp,".RDS"))
        
          model_runs <- readRDS(paste0("var_selection/model_runs/model_runs_",sp,".RDS"))
        
          gg_varimp <- get_variables_importance(model_runs)
        
          colnames(gg_varimp) <- c("id", 
                                   "PA.Run",
                                   "CV.Run", "Model", 
                                   "Variable", 
                                   "VI.run",
                                   "Variable.importance")
        
        
          gg_varimp$Variable <- reorder(gg_varimp$Variable,  
                                        gg_varimp$Variable.importance,
                                        median, 
                                        na.rm=TRUE)
        
        
        # hist(extract(mean_sst, p_dongh1[,1:2]), range = 100 )
        
          # ggplot(gg_varimp, aes(y = Variable, x = Variable.importance)) +
          #   geom_boxplot() + theme_classic() +ggtitle(sp)
          # 
          # ggplot(gg_varimp, aes(y = Variable, x = Variable.importance, color = Model)) +
          #   geom_boxplot() + theme_classic() +ggtitle(sp)
        
          if(!dir.exists(paste0("var_selection/var_import/")))
          {
            dir.create(paste0("var_selection/var_import/"), recursive = T)
          }
          
          
          # png(paste0("var_selection/all_vars/", sp, "/var_import/var_import_", iter, ".png"))
          # ggplot(gg_varimp, aes(y = Variable, x = Variable.importance)) +
          #   geom_boxplot() + geom_jitter(alpha = .5, aes(col = Model))  +  theme_classic() + ggtitle(paste0(sp, "_", iter))
          # dev.off()
          
          
          varimp <- ggplot(gg_varimp, aes(y = Variable, x = Variable.importance)) +
             geom_boxplot() + geom_jitter(alpha = .5, aes(col = Model))  +  theme_classic() + ggtitle(paste0(sp))
          plot(varimp)
          ggsave(paste0("var_selection/var_import/var_import_", sp, ".jpg"), plot = varimp, dpi = 300)
          dev.off()
          
          
          gg_varimp$Variable <- factor(gg_varimp$Variable, levels = rev(levels(gg_varimp$Variable)))
          
          
          quantile.50 <- aggregate(Variable.importance ~ Variable, data = gg_varimp, FUN = median)
        
          
          
          if(!dir.exists(paste0("var_selection/sel_vars/")))
          {
            dir.create(paste0("var_selection/sel_vars/"), recursive = T)
          }
        
          sel_vars<-list()
          sel_vars[[sp]] <- as.character(quantile.50$Variable[which(quantile.50$Variable.importance >= 0.1)])
          sel_vars 
          
          # Finally we remove as much variables than needed according to the presence/var number ratio
          npres <- length(which(presabs == 1)) 
          nvar_kept <- floor(npres/10)
          sel_varsf <- sel_vars
          
          
          if(!dir.exists(paste0("var_selection/varimp_df/")))
          {
            dir.create(paste0("var_selection/varimp_df/"), recursive = T)
          }
          
    
          if(length(sel_vars[[sp]]) > nvar_kept)
          {
            sel_varsf[[sp]] <- head(sel_vars[[sp]], nvar_kept)
            gg_varimp.f <- gg_varimp %>% 
              filter(grepl(paste(sel_varsf[[sp]], collapse = "|"), Variable))
          } else {
            gg_varimp.f <- gg_varimp %>% 
              filter(grepl(paste(sel_vars[[sp]], collapse = "|"), Variable))
          }
    
         
          saveRDS(sel_varsf, paste0("var_selection/sel_vars/sel_vars_f", sp, ".RDS"))
          saveRDS(gg_varimp.f, paste0("var_selection/varimp_df/varimp_df_", sp, ".RDS"))
          
          flush.console() # it empties the console before to go to the next iteration of the loop
          # otherwise it could create infinite loop 
        
    }
    

      # Summary of variables number per species
      
      nvar.df <- data.frame(species = character(0), # we build a summary dataframe to check the number of variables 
                            nvar = numeric(0))
      
      iter = 1
      for (sp in sp_list_aza)
      {
        
        gg_varimp.f <- readRDS(paste0("var_selection/varimp_df/varimp_df_", sp, ".RDS"))
        sel_varsf <- readRDS(paste0("var_selection/sel_vars/sel_vars_f", sp, ".RDS"))
        nvar <- length(sel_varsf[[1]])
        nvar.df[iter,1] <- sp
        nvar.df[iter,2] <- nvar
        iter = iter + 1
        
      }

# =================================================================================================
      
      




