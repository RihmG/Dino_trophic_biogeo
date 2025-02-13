
# =========== EXPLORATORY ANALYSIS ==============

# This script takes as input the abundance matrix of dinoflagellates species across sites and sites metadata 
# ("mothermat", created in script Z2) to compute multivariate analyses (diversity indexes, PCA, RDA) and create Figure 1. 

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
library(gridExtra)
library(grid)



# ----- DATA -------

  matfinalemoy2 <- readRDS("data_out/mothermat.RDS")


# ----- 1. Evenness ------

      #We first build an evenness graph 
      #to do so, some indices should be computed 
    
      #matfinalemoy2 <- as.data.frame(t(decostand(as.data.frame(t(matfinalemoy2)), method = 'total')))
      
      #Shannon index inverted
      Hinvsites<-diversity(matfinalemoy2, index = "shannon")
      
      #Inverted Species richness
      SRinvsites<-specnumber(matfinalemoy2)
      
      #Inverted Pielou's evenness 
      Jinvsites <- Hinvsites/log(SRinvsites)
    
      #Occupancy calculation
      matbin<-decostand(matfinalemoy2, method = 'pa')
      occupancy <- as.data.frame(rowSums(matbin))
      
      # occupancy %>%
      #   rename(occupancy = 'rowSums(matbin)') -> occupancy #occupancy is actually equal to inverted SR
      
      colnames(occupancy) <- 'occupancy'
      
      
      #Total abundance of species per sites
      abund.sites<-as.data.frame(rowSums(matfinalemoy2))
      
      # abund.sites %>%
      #   rename(Tot_abundance = 'rowSums(matfinalemoy2)') -> abund.sites
      
      colnames(abund.sites) <- 'Tot_abundance'
      
      
      #now we merge the data and add trophic information
      evenness_sites <-cbind(Jinvsites, occupancy, abund.sites)
      evenness_sites$species<-rownames(evenness_sites)
      
      evenness_sites_troph<-left_join(evenness_sites, trophic.sp, by = "species") 
      
      #we remove the "XXX" species, i.e those which are not taxonomically defined further than family
      evenness_sites_troph_X<-evenness_sites_troph[!grepl("X", evenness_sites_troph$species),]
      
      
      # Now we add the new trophic types from MDB Mitra 2023
      evenness_sites_troph_X23 <- evenness_sites_troph_X
      
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Dissodinium_pseudolunula"), "Trophy"] <- "mixoplankton"
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Dissodinium_pseudolunula"), "typeMX"] <- "CM"
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Gonyaulax_spinifera"), "Trophy"] <- "mixoplankton"
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Gonyaulax_spinifera"), "typeMX"] <- "CM"
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Kryptoperidinium_foliaceum"), "Trophy"] <- "mixoplankton"
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Kryptoperidinium_foliaceum"), "typeMX"] <- "CM" 
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Lepidodinium_chlorophorum"), "Trophy"] <- "mixoplankton"
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Lepidodinium_chlorophorum"), "typeMX"] <- "CM"
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Protoceratium_reticulatum"), "Trophy"] <- "mixoplankton"
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Protoceratium_reticulatum"), "typeMX"] <- "CM"  
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Torodinium_teredo"), "Trophy"] <- "mixoplankton"
      evenness_sites_troph_X23[which(evenness_sites_troph_X23$species == "Torodinium_teredo"), "typeMX"] <- "CM"  
      
      
      # Now we add the new_mixo dataframe 
      
      new_mixo <- readRDS("data_out/new_mixo.RDS")
      
      # for loop 
      
      evenness_sites_troph_X24 <- evenness_sites_troph_X23
      for(i in 1:(nrow(new_mixo))){
        
        sp = new_mixo[i,"Species.Name"]
        trophy = new_mixo[i, "Trophy"]
        typemx = new_mixo[i, "typeMX"]
        evenness_sites_troph_X24[which(evenness_sites_troph_X24$species == sp), "Trophy"] <- trophy
        evenness_sites_troph_X24[which(evenness_sites_troph_X24$species == sp), "typeMX"] <- typemx
        
      }
      
      # we remove empty spaces 
      evenness_sites_troph_X24$Trophy <- gsub(" ", "", evenness_sites_troph_X24$Trophy)
      
      
      # we change names of trophy 
      evenness_sites_troph_X24$Trophy <- as.factor(evenness_sites_troph_X24$Trophy)
      levels(evenness_sites_troph_X24$Trophy) <- c("mixotroph", "phototroph", "phagotroph")
      
      
      saveRDS(evenness_sites_troph_X24, "evenn_troph.RDS")
      
      evenness.plot2 <-ggplot(evenness_sites_troph_X24, aes(x = occupancy, 
                                                            y = Jinvsites, 
                                                            size = Tot_abundance, 
                                                            color = Trophy)) + 
        geom_point()+ 
        geom_text_repel(subset(evenness_sites_troph_X24, occupancy > 200 | species == "Dinophysis_acuminata"), 
                        mapping = aes(x = occupancy,
                                      y= Jinvsites,
                                      label = sp_name), 
                        size= 3, 
                        fontface = "italic",
                        point.padding = 4) +
        scale_color_manual(values = c("mixotroph" = "brown2",
                                      "phototroph"="chartreuse3",
                                      "phagotroph"="blue"),
                           name = "Trophic type") +
        scale_size(name = "Total reads abundance") +
        ggtitle("Species evenness among sampling locations", ) +
        xlab("Occupancy") + ylab("Pielou's evenness") + 
        theme_classic() + 
        theme(legend.position = c(0.80, 0.20), legend.box = "horizontal",
              legend.background = element_blank()) 
      
      # +
      #   scale_y_continuous(expand = c(0.005, 0.01)) 
      # 
      
      ggsave("evenn_plot.pdf", plot = evenness.plot2, width = 12, height = 8)
      
      # Let's check the distribution of the evenness (Y), and occupancy (X) 
      
      densXX <- ggplot(evenness_sites_troph_X24, aes(x=occupancy, color = Trophy)) +
        geom_density(alpha=0.5) +
        scale_color_manual(values = c("mixotroph" = "brown2",
                                      "phototroph"="chartreuse3",
                                      "phagotroph"="blue")) +
        theme_classic() +
        theme(axis.title.x=element_blank())+
        theme(axis.title.y=element_blank())+
        theme(axis.text.y=element_blank()) +
        theme(axis.ticks.y=element_blank()) +
        theme(legend.position = "none") +
        theme(plot.margin = margin(0, 0, 0, 0)) +
        scale_x_continuous(limits = range(evenness_sites_troph_X24$occupancy))
    
      
      densYY <- ggplot(evenness_sites_troph_X24, aes(x=Jinvsites, color = Trophy)) +
        geom_density(alpha=0.5) +
        scale_color_manual(values = c("mixotroph" = "brown2",
                                      "phototroph"="chartreuse3",
                                      "phagotroph"="blue")) +
        theme_classic() +
        coord_flip() +
        theme(axis.title.y=element_blank()) +
        theme(axis.title.x=element_blank()) +
        theme(axis.text.x=element_blank()) +
        theme(axis.ticks.x=element_blank()) +
        theme(legend.position = "none") +
        theme(plot.margin = unit(c(0,0.2,0.6,0), "cm")) +
        scale_x_continuous(limits = range(evenness_sites_troph_X24$Jinvsites))
      
      p_void <- ggplot() + theme_void()
      
      
      ggarrange(densXX, p_void, evenness.plot2, densYY, 
                ncol = 2, nrow = 2, 
                align = "hv", 
                widths = c(4.5, 1), heights = c(1, 2.5)) -> f_evenness_p
      
      ggsave("last_figures/evenn_f.pdf", f_evenness_p, width = 13.5, height = 11)
    
      
      # Tests 
        
        # Kruskal-Wallis test 
        
        # Occupancy
        occ.kruskal <- evenness_sites_troph_X24 %>% kruskal_test(occupancy ~ Trophy)
        # No difference 
        
        # Evenness
        eve.kruskal <- evenness_sites_troph_X24 %>% kruskal_test(Jinvsites ~ Trophy)
        # No difference
        
        # Abundance
        abu.kruskal <- evenness_sites_troph_X24 %>% kruskal_test(Tot_abundance ~ Trophy)
        
        
        # Wilcoxon test 
        
        # Occupancy
        occ.wilcox <- evenness_sites_troph_X24 %>% 
          wilcox_test(occupancy ~ Trophy)
        # No difference 
        
        # Evenness
        eve.wilcox <- evenness_sites_troph_X24 %>% 
          wilcox_test(Jinvsites ~ Trophy, p.adjust.method = "bonferroni")
        # No difference 
        
        # Abundance
        abu.wilcox <- evenness_sites_troph_X24 %>% 
          wilcox_test(Tot_abundance ~ Trophy, p.adjust.method = "bonferroni")
        
        ggdensity(evenness_sites_troph_X24, x = "Tot_abundance",
                  add = "mean", rug = TRUE,
                  color = "Trophy", fill = "Trophy") +
                  xlim(0,100000)
        
        
        evenness_sites_troph_X24$cat <- ifelse(evenness_sites_troph_X24$Jinvsites < 100, "A",
                               ifelse(evenness_sites_troph_X24$Jinvsites >= 100 & evenness_sites_troph_X24$Jinvsites <= 200, "B", "C"))
        
        evenness_sites_troph_X24 <- evenness_sites_troph_X24 %>%
          mutate(cat = case_when(
            occupancy < 100 ~ "A",
            occupancy >= 100 & occupancy <= 200 ~ "B",
            occupancy > 200 ~ "C"
          ))
        
        
        levene_test <- leveneTest(Jinvsites ~ cat, data = evenness_sites_troph_X24)
        # p value < 0.05, not homogenous variances
        
        kruskal_test <- kruskal_test(Jinvsites ~ cat, data = evenness_sites_troph_X24)
        
        wilcox_test <- wilcox_test(Jinvsites ~ cat, data = evenness_sites_troph_X24)
        
        JPCA <- PCA(evenness_sites_troph_X24[,(1:3)])
        hcpc_J <- HCPC(JPCA)
        Jclust <- hcpc_J$data.clust
        
        nrow(Jclust)
        length(unique(Jclust$Tot_abundance))
        
        Jclust2 <- left_join(Jclust, evenness_sites_troph_X24[,3:6])
        
        wilcox_test2 <- wilcox_test(Jinvsites ~ clust,  data = Jclust)
        
        
        # Summary statistics
        nrow(Jclust2[which(Jclust2$clust == 2),])
        summary(Jclust2[which(Jclust2$clust == 2),]$occupancy)
        sd(Jclust2[which(Jclust2$clust == 2),]$Jinvsites)
        sd(Jclust2[which(Jclust2$clust == 2),]$occupancy)
        
        # turn NaN into NA 
        evenness_sites_troph_X24[] <- lapply(evenness_sites_troph_X24, function(col) replace(col, is.nan(col), NA))

        
        # package for plotting
        
        library(ggpubr)
        
        compare_means(Jinvsites ~ clust,  data = Jclust)
        
        jplot <- ggboxplot(Jclust2, x = "clust", y = "Jinvsites",
                       color =  "clust", palette = "jco",
                       add = "jitter") + 
          stat_compare_means(comparisons = list( c("1", "2"), c("1", "3"), c("2", "3") )) +
          xlab("Group") + ylab("Pielou's Evenness value") + labs(color = "Cluster") + theme(legend.position = "none")
        
          
        
        # PLot with the clusters 
        
        library(paletteer)
        
        evenness.plot_clust <- ggplot(Jclust2, aes(x = occupancy, 
                                                              y = Jinvsites, 
                                                              size = Tot_abundance, 
                                                              color = clust)) + 
          geom_point() + 
          scale_color_jco(labels = c("1" = "Rare", "2" = "Intermediate", "3" = "Ubiquitous")) +
          # scale_color_paletteer_d("fishualize::Acanthurus_olivaceus", 
          #                         labels = c("1" = "Endemic",
          #                                   "2" = "Intermediate",
          #                                   "3" = "Ubiquitous"
          #                                   )) +
          #ggtitle("Species evenness among sampling locations, colored according their \n ubiquitousness profiles") +
          xlab("Occupancy") + ylab("Pielou's evenness") + labs(size = "Total reads abundance", color = "Clusters") +
          theme_classic() +
          xlim(0,800) +
          ylim(0,1)
        
        
        # ====== Figure S4 =======
        
        grid.arrange(jplot, evenness.plot_clust, ncol = 2) -> evenn_clust
        
        ggsave("suppmat/evenn_clust.pdf", evenn_clust, width = 12, height = 7)
        # ..........................................................................
        
        
        
        # Creating datatables
      
        stargazer(as.data.frame(occ.wilcox),
                  summary = F,
                  rownames = F,
                  title = "Wilcoxon tests for occupancy distribution between trophic types")
        
        stargazer(as.data.frame(eve.wilcox),
                  summary = F,
                  rownames = F,
                  title = "Wilcoxon tests for evenness distribution between trophic types")
        
        stargazer(as.data.frame(abu.wilcox),
                  summary = F,
                  rownames = F,
                  title = "Wilcoxon tests for abundance distribution between trophic types")
        
        
        # Sampling periods 
        
        date<-as.data.frame(data2511[,"date"])

        colnames(date) <- 'date'
        
        date$date <- as.numeric(substr(date$date, 1, 4))
        
        date_sample <- ggplot(date, aes(x = date)) + 
          geom_bar(width = 0.5, colour = 'black', fill = 'grey' ) +
          ggtitle("MetaPR2 extracted data sampling period") +
          xlab("Years") + ylab("Number of samples") + 
          theme_classic() 
        
        my_cols <- c("firebrick", "limegreen", "blue3")
      
        readRDS("numsamples.RDS") -> numsamples

        # ======= Figure S2 =====
        
        ggsave("suppmat/temp.pdf", date_sample, width = 8, height = 5)
        ggsave("suppmat/spat.pdf", numsamples, width = 10, height = 5)
        
        
# =========================================================================================
        
        
        
        
        
        
        
        
        





