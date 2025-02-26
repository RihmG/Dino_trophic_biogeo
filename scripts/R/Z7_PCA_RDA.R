
# =========== Environmental PCA - RDA ==============

# This script takes as input the abundance matrix of dinoflagellates species across sites and sites metadata
# (environment) to compute multivariate analyses (PCA, RDA) and create figure 2 plot.

# Libraries 
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggsci)
library(funrar)
library(vegan)
library(ggrepel)
library(pastecs)
library(ggalt)
library(ggpubr)
library(rstatix)
library(stargazer)


# ------- 0. Function to remove empty rows/columns -------

    rm_rowzeros <- function(mat){
      mat$sum = rowSums(mat)
      mat %>%
        subset(sum > 0) -> mat
      mat[, -which(names(mat) %in% c("sum"))] -> mat
      return(mat)
    }
    
    rm_colzeros <- function(mat){
      mat = as.data.frame(t(mat))
      mat$sum = rowSums(mat)
      mat %>%
        subset(sum > 0) -> mat
      mat[, -which(names(mat) %in% c("sum"))] -> mat
      mat = as.data.frame(t(mat))
      return(mat)
    }
    

# ----- 1. DATA -------

    mothermat <- readRDS("data_out/mothermat.RDS")
    new_mixo <- readRDS("data_out/new_mixo.RDS")
    

# ----- 2. PCA ------- 

    # Climatologies 
    
      base_files <- list.files("data_out/Env_data/monthly_base/")
      
      for(i in base_files)
      {
        base_object <- readRDS(paste0("data_out/Env_data/monthly_base/",i))
        assign(gsub(".RDS", "", i), base_object)
      }
      
      baseline_logmean.NA <- readRDS("data_out/Env_data/baseline_logmean.NA.RDS")
      
    # Sampling locations 
      
      sites2<-read.table("Lat_Long_Proj_ProjD_NbS_List2_2023_2", sep = "\t", h=T)# from a .sh script (-->README)
      sites2<-sites2[-1,]
      rownames(sites2) <- 1:nrow(sites2)
      
      sites2_coord <- rev(sites2[,c(1:2)])
      sites2_coord$site <- paste0("site", rownames(sites2_coord))
      
    # Extract env. values for samples coordinates
      
      env.extract <- as.data.frame(extract(baseline_logmean.NA, sites2_coord[, c("longitude", "latitude")]))
      env.extract.lab <- env.extract
      env.extract.lab$site <- paste0("site", rownames(env.extract.lab))
      #env.extract[is.nan(env.extract)] <- NA # remove NaN
      
      env.extract <- as.data.frame(na.omit(env.extract)) # remove all rows containing NA
      env.extract.lab <- as.data.frame(na.omit(env.extract.lab)) 
      
      sites_extract <- left_join(env.extract.lab, sites2_coord, by = "site")
      
    # PCA 
      
      # standardization
      env.extract <- decostand(env.extract, method = 'standardize')
      
      # check parameters 
      round(apply(env.extract, 2, mean), 1) # mean = 0
      apply(env.extract, 2, sd) # standard deviation = 1
      
      pca.NA <- PCA(env.extract, ncp = 5)

      pca.NA_coord <- as.data.frame(pca.NA$ind$coord)
      
      ggplot(pca.NA_coord, aes(x = Dim.1, y = Dim.2)) + 
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_point()+
        theme_classic()
      
      # Clusterisation of points 
      
      hcpc.NA <- HCPC(pca.NA, graph = T)
      
      # CLuster graphs 
      fviz_cluster(hcpc.NA,
                   repel = TRUE,           
                   show.clust.cent = TRUE, 
                   palette = "jco",         
                   ggtheme = theme_minimal(),
                   main = "Factor map"
                   # ,
                   # labelsize = F
      )
      
      fviz_dend(hcpc.NA, show_labels = FALSE, k_colors = c("#D43F3AFF",
                                                           "#EEA236FF",
                                                           "#5CB85CFF",
                                                           "#46B8DAFF")) -> dendro_clust # Figure S5
      
      
      
      # We match the sites coordinates with the retained sites from the pca
      
      sites_clust.NA <- cbind(sites_extract[,c("site","longitude", "latitude")], hcpc.NA$data.clust$clust)
      
      colnames(sites_clust.NA)[4] <- 'clust'
      
      #we project the points on a map to see clusters
      
      
      # ====== Figure 2A =======
      
      ggplot() +
        geom_sf(data = worldmap, color = "grey", fill = "grey") +
        theme(panel.background = element_rect(fill = "aliceblue"),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        geom_point(sites_clust.NA, mapping = aes(x = longitude, y = latitude, color = clust)) +
        scale_color_locuszoom() + 
        coord_sf(expand = F) -> rda_map
      
      ggsave("last_figures/rda_map.pdf",rda_map , width = 9, height = 6)
      #...............................................................................................
      
      
      pca.NA_coords <- pca.NA$ind$coord 
      pca.NA_coords_clust <- cbind(pca.NA_coords, as.data.frame(hcpc.NA$data.clust$clust))
      
      colnames(pca.NA_coords_clust)[6] <- 'clust'
      
      
      # ====== Figure S5 ======
      ggplot(pca.NA_coords_clust, mapping = aes(x = Dim.1, y = Dim.2, color = clust)) + 
        geom_point()+
        scale_color_locuszoom() +
        # stat_ellipse() +
        xlab("PCA1 (58.80%)") +  ylab("PCA2 (14.80%)") +
        theme_classic() + 
        labs(color = "Cluster")  -> color_clust
  
      ggarrange(color_clust, dendro_clust, ncol = 2) -> dendro_S5
      
      ggsave("suppmat/dendro_S5.pdf", dendro_S5, width = 15, height = 7)
      #..............................................................................................
          
  # ----- 3. RDA -------
      
      # dataset standardization
      
      mothermat_hel <- decostand(t(mothermat), method = "hellinger")
      
      # --> log_pca_data is 744 sites long
      # --> mothermat is 895 sites long
      
      rda_mat.NA <- as.data.frame(mothermat_hel)
      
      # We remove "X" species 
      rda_mat.NA <- rda_mat.NA[,!grepl("X", colnames(rda_mat.NA))]
      
      rownames(rda_mat.NA) <- gsub(" ", "", rownames(rda_mat.NA))
      
      # So we select pixels of abundance matrix which are present in env. matrix 
      #rownames(log_pca_data) <- paste("site", rownames(log_pca_data))
      
      rda_matf.NA <- rda_mat[which(rownames(rda_mat.NA) %in% env.extract.lab$site),]
      
      # check for 0 in rows and cols
      rda_matf.NA <- rm_colzeros(rda_matf.NA)
      rda_matf.NA <- rm_rowzeros(rda_matf.NA) # 747 sites, 233 species after cleaning 
      
      # We can now make the rda 
      
      # We first select species with Escoufier's equivalent vector method 
      
      # Variables selection
      
      rda_matf.sel.NA <- escouf(rda_matf.NA)
      
      plot(rda_matf.sel.NA)
      rda_matf.sel.NA$level <- 0.90
      rda_matf.sel.NA.df <- pastecs::extract(rda_matf.sel.NA)

      # We remove the potential rows that are not present in the new selected matrix
      
      rownames(rda_matf.sel.NA.df) <- gsub(" ", "", rownames(rda_matf.sel.NA.df))
      
      which(env.extract.lab$site %in% rownames(rda_matf.sel.NA.df) == F) # 216, 245, 576
      rda.env.data <- env.extract.lab[-which(env.extract.lab$site %in% rownames(rda_matf.sel.NA.df) == FALSE),]
      
      rdaf.NA <- rda(rda_matf.sel.NA.df ~ ., rda.env.data[,-12]) # we dont add the index as an explanatory variable
      
      RsquareAdj(rdaf)
      summary(rdaf)
      
      # rdaf.step = ordistep(rda(rda_matf.sel ~ 1, data = log_pca_data), 
      #                          scope=formula(rdaf), 
      #                          direction="both", 
      #                          pstep=5000)
      
      rdaf.step.NA = ordistep(rda(rda_matf.sel.NA.df ~ 1, data = rda.env.data[,-12]), 
                              scope=formula(rdaf.NA), 
                              direction="both", 
                              pstep=5000)
      
      # rdaff <- rda(formula = formula(rdaf.step), data = log_pca_data)
      
      rdaff.NA <- rda(formula = formula(rdaf.step.NA), data = rda.env.data[,-12])
      plot(rdaff.NA)
      RsquareAdj(rdaff.NA)
      
      # Permutation test 
      anova.cca(rdaff.NA, step = 1000) # significative 
      rdaff.NA$CA$eig[rdaff.NA$CA$eig > mean(rdaff.NA$CA$eig)]
      
      # Percentages per RDA axes
      perc.NA <- round(100*(summary(rdaff.NA)$cont$importance[2, 1:2]), 2)
      
      
      # We add latitude/longitude as covariates to run partial RDA
      
      sites2_coord$sites <- paste("site", rownames(sites2_coord))
      
      pca.esc.part <- left_join(pca.esc, sites2_coord, by = "sites")
      pca.esc.part1 <- pca.esc.part[,1:9] # env. data 
      pca.esc.part2 <- pca.esc.part[,11:12] # coordinates data 
      
      rdap <- rda(rda_matf.sel, pca.esc.part1, pca.esc.part2)
      summary(rdap)
      
      anova.cca(rdap, step = 1000)
      
      RsquareAdj(rdap)
      round(100*(summary(rdap)$cont$importance[2, 1:2]), 2)
      
      
      # Variance paritionning 
      
      varpart <- varpart(rda_matf.sel, pca.esc.part1, pca.esc.part2)
      varpart$part
      
      plot(varpart,
           Xnames = c("Env.", "Space"), # name the partitions
           bg = c("green", "blue"), alpha = 80, # colour the circles
           digits = 2, # only show 2 digits
           cex = 1.5)
      
      anova.cca(rda(rda_matf.sel, pca.esc.part1))
      anova.cca(rda(rda_matf.sel, pca.esc.part2))
      anova.cca(rda(rda_matf.sel, pca.esc.part2))
      
      
      stations_rda <- as.data.frame(scores(rdaff, scaling=2)$sites)
      species_rda <- as.data.frame(scores(rdaff, scaling=2)$species)
      species_rda$species <- rownames(species_rda)
      env_rda <- scores(rdaff, choices = 1:2, display="bp", scaling=2) # we get rid of all qualitative variables
      
      
      
      # ===== RDA Plot ======
      
      
      plot(rdaff.NA, scaling = 2, main = "RDA - Scaling II")
      
      
      stations_rda.NA <- as.data.frame(scores(rdaff.NA, scaling=2)$sites)
      species_rda.NA <- as.data.frame(scores(rdaff.NA, scaling=2)$species)
      species_rda.NA$species <- rownames(species_rda.NA)
      env_rda.NA <- scores(rdaff.NA, choices = 1:2, display="bp", scaling=2) # we get rid of all qualitative variables
      
      sites_clust.rda.NA <- sites_clust.NA[which(sites_clust.NA$site %in% rda.env.data$site),]
      
      #we project the points on a map to see clusters
      
      ggplot() +
        geom_sf(data = worldmap, color = "grey", fill = "grey") +
        theme(panel.background = element_rect(fill = "aliceblue")) +
        geom_point(sites_clust.rda.NA, mapping = aes(x = longitude, y = latitude, color = clust)) +
        scale_color_locuszoom() 
      
      # OK
      
      # And now we join the cluster information on RDA scores dataset
      
      stations_rda.NA$site <- rownames(stations_rda.NA)
      
      stations_rda.NA2 <- left_join(stations_rda.NA, sites_clust.rda.NA, by = "site")
      
      g1.NA <- ggplot() + 
        geom_point(stations_rda.NA2, mapping = aes(x = RDA1, y = RDA2, color = clust), alpha = 0.5)+
        scale_color_locuszoom(name = "Environmental clusters") +
        xlab("RDA1 (26.08%)") +  ylab("RDA2 (3.07%)") +
        theme_classic()
      
      species_rda.NA2 <- left_join(species_rda, evenness_sites_troph_X24[,c(4:6)], by = c("species"))
      
      # We create a palette for trophic types
      color.troph <- as.character(species_rda.NA2$Trophy)
      color.troph[which(color.troph == "mixotroph")] = "firebrick1"
      color.troph[which(color.troph == "phagotroph")] = "dodgerblue2"
      color.troph[which(color.troph == "phototroph")] = "olivedrab3"
      color.troph[which(is.na(color.troph))] = "grey"
      
      
      # We select the species to label, according to their contribution to the RDA1 and RDA2 variance explanation
      
      ggplot(species_rda.NA2, aes(x = reorder(species, abs(RDA1)), y = abs(RDA1), label = species)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        labs(x = "species",
             y = "Absolute RDA1 values") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      
      ggplot(species_rda.NA2, aes(x = reorder(species, abs(RDA2)), y = abs(RDA2), label = species)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        labs(x = "species",
             y = "Absolute RDA2 values") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      sp_to_label2 <- c("Prorocentrum_sp.", "Gyrodinium_fusiforme", "Gyrodinium_dominans", "Torodinium_robustum",
                       "Warnowia_sp.", "Lepidodinium_chlorophorum", "Gyrodinium_sp.", "Islandinium_tricingulatum", 
                       "Heterocapsa_sp.", "Tripos_furca", "Karlodinium_veneficum", "Tripos_fusus")
      
      sp_scores_to_label2 <- species_rda.NA2[species_rda.NA2$species %in% sp_to_label2, ]
      
      sp_scores_to_label2$sp_name <- gsub("_", " ", sp_scores_to_label2$species)
      
      sp_scores_to_label3 <- sp_scores_to_label2[,-3]
      sp_scores_to_label3$RDA1 <- round(sp_scores_to_label3$RDA1, 2)
      sp_scores_to_label3$RDA2 <- round(sp_scores_to_label3$RDA2, 2)
      sp_scores_to_label3$Trophy <- as.character(sp_scores_to_label3$Trophy)
      sp_scores_to_label3 <- as.data.frame(lapply(sp_scores_to_label3, function(x) {replace(x, is.na(x), "-")}))
      sp_scores_to_label3 <- sp_scores_to_label3[,c("sp_name", "Trophy", "typeMX", "RDA1", "RDA2")]
      
      colnames(sp_scores_to_label3) <- c("Species", "Trophic type", "Mixotrophic type", "RDA1", "RDA2")
      
      
      # Table generation
      stargazer(sp_scores_to_label3,
                summary = F,
                rownames = F,
                title = "XXX")
      
   
      # We add species rows 
      g2.NA <- g1.NA  + geom_segment(aes(xend = species_rda.NA2$RDA1, yend = species_rda.NA2$RDA2),x=0,y=0,
                               arrow = arrow(length = unit(0.2,"cm")),size = 0.5, color=color.troph) +
        geom_text_repel(data = sp_scores_to_label2, 
                        aes(x = RDA1, y = RDA2, label = sp_name), 
                        size = 4,
                        color = c("firebrick1",
                                  "firebrick1",
                                  "grey",
                                  "grey",
                                  "firebrick1",
                                  "dodgerblue2",
                                  "dodgerblue2",
                                  "dodgerblue2",
                                  "dodgerblue2",
                                  "firebrick1",
                                  "firebrick1",
                                  "grey"),
                        fontface = "italic"
        ) 
      #geom_text(aes(env_rda[,1]*1.7, env_rda[,2]*1.7, label=rownames(env_rda)), color="steelblue4")
      
      
      
      
      # we ad the 0 centered  axes
      g3.NA <- g2.NA +  geom_hline(yintercept = 0, linetype='dotted') +
        geom_vline(xintercept = 0, linetype='dotted') +
        labs(title="RDA") +
        theme(plot.title=element_text(hjust=0.5))
      
      
      
      env_rda.NA.df <- as.data.frame(env_rda.NA)
      env_rda.NA.df$variables <- rownames(env_rda.NA.df)
      env_rda.NA.df[env_rda.NA.df$variable == "bio_chla"]
      
      #Change the names of predictors for plotting 
      env_rda.NA.df$var_names <- c("SST", "Chla", "SSS", "O2", "PO4", "Sistar", "Nstar","PAR", "MLD", "NO3", "Si")
      
      
      g4.NA <- g3.NA + geom_segment(data=as.data.frame(env_rda.NA.df),
                              aes(xend = RDA1*1.5, yend = RDA2*1.5),
                              linetype = "dashed",
                              x=0,y=0,size = 0.5,
                              color = 'steelblue4',
                              arrow = arrow(length = unit(0.2,"cm"), ends = "last", type = "closed"),
                              arrow.fill = 'steelblue4') +
        geom_text(aes(env_rda.NA.df[,1]*1.7, env_rda.NA.df[,2]*1.7, label= env_rda.NA.df$var_names), color="steelblue4"
                  # geom_text_repel(data=as.data.frame(env_rda),
                  #                 aes(x = RDA1, y = RDA2, label = rownames(env_rda)), 
                  #                 size = 4,
                  #                 color = "steelblue4"
        ) 
      
      g4.NA
      
      ggsave(g4.NA, filename = "graphiques/rda.NA.svg")
      ggsave(g4.NA, filename = "last_figures/rda.NA.pdf", width = 13.5, height = 9)
      
   
# ----------------------------------------------------------------------------
      
      
      
      











