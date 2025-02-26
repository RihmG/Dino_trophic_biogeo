
#--------- SUMMARY PLOT -----------

#Inspired from Malviya et al., 2016

# This script creates plots that are summaries of the abundance dataframe, with per species:
#   - abundance of reads 
#   - number of unique ASVs
#   - Depth category (Surface/Euphotic)
#   - Size Fraction
#   - Occupancy among sites

# --> Figure S3

# Libraries

library(tidyverse) # /!\ conflicts between dplyr and plyr (e.g. function count()...)
library(egg)
library(cowplot)



# ==== 1. ASV per species ====


#metaPR2 data, downloaded 25/11/2022

data2511<- read.table(file = 'metapr2_ASVs_selected_abundance_Dinophyceae_2022-11-25.tsv',
                      sep = '\t', header = TRUE)


  #we first need to count the number of different ASVs per species
  data2511 %>%
    group_by(species) %>%
    count(as.factor(asv_code)) -> count.asv
  
  count.asv %>%
    rename(asv_code = 'as.factor(asv_code)') -> count.asv
  
  
  count.asv %>%
    group_by(species) %>%
    summarise(n_asv=length(asv_code), .groups = 'drop') -> count.asv.sp
    

# ==== 2. Reads per species ====

  #then we count the number of reads per species
  data2511 %>%
    group_by(species) %>%
    count(n_reads) -> count.reads
  
  count.reads$total <- count.reads$n_reads * count.reads$n
  
  count.reads %>%
    group_by(species) %>%
    summarise(ntot=sum(total),
              .groups = 'drop') -> count.reads.sp


# ==== 3. Summary of reads per depth per species ====

  #then we summarise the depths
  data2511 %>%
    group_by(species) %>%
    count(depth_level) -> sum_depth
  
  list.data2511 <- split(data2511, data2511$species)
  
  
  #we create a function which counts the number of reads per depth_level (surface/euphotic), per species
  
  count.depths.reads <- function(i) {
    list.data2511[[i]] %>%
      group_by(depth_level) %>%
      summarise(tot_reads=sum(n_reads),
                .groups = 'drop') -> count
    return(count)
  
  }
  
  count.depths.reads.final<-lapply(as.list(1:264), count.depths.reads)
  
  reads_per_depth<-do.call(rbind, count.depths.reads.final)
  
  # now we have the nummber of observations AND the total number of reads per surface per depth_level
  sum_depth2<-cbind(sum_depth, reads_per_depth)
  
  sum_depth3<-sum_depth2[,c(1,2,3,5)]
  
  sum_depth3 %>%
    rename(
      depth_level = 'depth_level...2',
      n_obs = 'n'
    ) -> sum_depth3

# ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
  split.depth <- split(sum_depth3, sum_depth3$species)
  
  #we create a function which takes the sub-table of the list and return a dataframe which resume within
  #three columns the previous dataframe (i.e instead of having one line for each depth level, we 
  #create one column 'euphotic' and one 'surface' and store the information inside) We also create 
  #a total column to make percentage calculation easier
  
  sum_split_depth <- function(i){
    species <- c()
    surface <- c()
    euphotic <-c()
    tot <- c()
    #we need to be careful of the length of each sub-table : some species are only 'surface' or 'euphotic'
    #so we check if the sub-table has 2 or only one depth_level and store a 0 in case no individuals are 
    #found in the following level
    if (length(split.depth[[i]][split.depth[[i]]$depth_level == 'surface',]$tot_reads) == 0){
      species <-c(species, unique(split.depth[[i]]$species))
      surface <- c(surface, 0)
      euphotic <- c(euphotic, split.depth[[i]][split.depth[[i]]$depth_level == 'euphotic',]$tot_reads)
      tot <- c(tot, 0 + 
                 split.depth[[i]][split.depth[[i]]$depth_level == 'euphotic',]$tot_reads)
    }
    else if (length(split.depth[[i]][split.depth[[i]]$depth_level == 'euphotic',]$tot_reads) == 0){
      species <-c(species, unique(split.depth[[i]]$species))
      surface <- c(surface, split.depth[[i]][split.depth[[i]]$depth_level == 'surface',]$tot_reads)
      euphotic <- c(euphotic, 0)
      tot <- c(tot, split.depth[[i]][split.depth[[i]]$depth_level == 'surface',]$tot_reads + 
                 0)
    }
    else {
      species <-c(species, unique(split.depth[[i]]$species))
      surface <- c(surface, split.depth[[i]][split.depth[[i]]$depth_level == 'surface',]$tot_reads)
      euphotic <- c(euphotic, split.depth[[i]][split.depth[[i]]$depth_level == 'euphotic',]$tot_reads)
      tot <- c(tot, split.depth[[i]][split.depth[[i]]$depth_level == 'surface',]$tot_reads + 
                 split.depth[[i]][split.depth[[i]]$depth_level == 'euphotic',]$tot_reads )
    }
    return(data.frame(species, surface, euphotic, tot))
  }
  
  
  final_depth<-lapply(as.list(1:264), sum_split_depth)
  
  final_depth.df<-do.call(rbind, final_depth)
  
  #finally we compute percentages 
  
  final_depth.df %>%
    mutate(perc_surf = round((100*surface)/tot, digits = 2),
           perc_euph = round((100*euphotic)/tot, digits = 2)) -> final_depth.df
  

# ==== 4. Summary of reads per size per species ====

  # Then we summarise the samples size fractions
  
  data2511 %>%
    group_by(species) %>%
    count(fraction_name) -> sum_size
  
  
  count.size.reads <- function(i) {
    list.data2511[[i]] %>%
      group_by(fraction_name) %>%
      summarise(tot_reads=sum(n_reads),
                .groups = 'drop') -> count
    return(count)
    
  }
  
  count.size.reads.final<-lapply(as.list(1:264), count.size.reads)
  
  reads_per_size<-do.call(rbind, count.size.reads.final)
  
  # now we have the nummber of observations AND the total number of reads per surface per depth_level
  sum_size2<-cbind(sum_size, reads_per_size)
  
  sum_size3<-sum_size2[,c(1,2,3,5)]
  
  sum_size3 %>%
    rename(
      fraction_name = 'fraction_name...2',
      n_obs = 'n'
    ) -> sum_size3
  
  
  split.size <- split(sum_size3, sum_size3$species)
  
  
  #we create a function which takes the rank i of the species within the list and add the information of the list
  #inside a 264 x 7 matrix initially full of zeros
  
  sum_split_size <- function(i){
  
    df0<-as.data.frame(matrix(0, 264, 7))
    colnames(df0) <- c("species", "micro", "nano", "nano-micro", "pico", "pico-nano", "total")
    df0$species[i] <- unique(split.size[[i]]$species)
    
    #for each column of the matrix (i.e size fraction), it compares the content of the sub-tables of the list 
    #with the name of the matrix, if the function finds out that the name of the column is present in the 
    #sub-table, it fills the corresponding cell of the matrix with the good information
    
    if (isTRUE(colnames(df0)[2] %in% split.size[[i]]$fraction_name)){
      df0$micro[i] <- split.size[[i]][split.size[[i]]$fraction_name == 'micro',]$tot_reads
    }
  
    if (isTRUE(colnames(df0)[3] %in% split.size[[i]]$fraction_name)){
      df0$nano[i] <- split.size[[i]][split.size[[i]]$fraction_name == 'nano',]$tot_reads
    }
  
    if (isTRUE(colnames(df0)[4] %in% split.size[[i]]$fraction_name)){
      df0$`nano-micro`[i] <- split.size[[i]][split.size[[i]]$fraction_name == 'nano-micro',]$tot_reads
    }
  
    if (isTRUE(colnames(df0)[5] %in% split.size[[i]]$fraction_name)){
      df0$pico[i] <- split.size[[i]][split.size[[i]]$fraction_name == 'pico',]$tot_reads
    }
  
    if (isTRUE(colnames(df0)[6] %in% split.size[[i]]$fraction_name)){
      df0$`pico-nano`[i] <- split.size[[i]][split.size[[i]]$fraction_name == 'pico-nano',]$tot_reads
    }
    
    if (isTRUE(colnames(df0)[7] %in% split.size[[i]]$fraction_name)){
      df0$total[i] <- split.size[[i]][split.size[[i]]$fraction_name == 'total',]$tot_reads
    }
  
  return(df0)
    
  }
  
  final_size<-lapply(1:264, sum_split_size)
  
  final_size.df<-do.call(rbind, final_size)
  
  #Problem : the final df is 264 times too big, so we remove all the rows where the species names are not given
  #(means that the row is full of zeros)
  
  final_size.df %>%
    subset(species != 0) -> final_size.df
  
  # Finally we compute percentages 
  
  final_size.df$tot <- final_size.df$micro + final_size.df$nano + final_size.df$`nano-micro` + 
                          final_size.df$pico + final_size.df$`pico-nano` + final_size.df$total
  
  final_size.df %>%
    mutate(perc_micro = round((100*micro)/tot, digits = 2),
           perc_nano = round((100*nano)/tot, digits = 2),
           perc_nanomicro = round((100*`nano-micro`)/tot, digits = 2),
           perc_pico = round((100*pico)/tot, digits = 2),
           perc_piconano = round((100*`pico-nano`)/tot, digits = 2),
           perc_total = round((100*total)/tot, digits = 2)) -> final_size.df
  
  
  # We do the same but we split the percentages between total and non-total
  
  final_size.df2<-do.call(rbind, final_size)
  
  final_size.df2 %>%
    subset(species != 0) -> final_size.df2
  
  final_size.df2$tot <- final_size.df2$micro + final_size.df2$nano + final_size.df2$`nano-micro` + 
    final_size.df2$pico + final_size.df2$`pico-nano` + final_size.df2$total
  
  final_size.df2$non_tot <- final_size.df2$micro + final_size.df2$nano + final_size.df2$`nano-micro` + 
    final_size.df2$pico + final_size.df2$`pico-nano`
  
  
  # And now we compute the perc depending on total/non-total
  
  final_size.df2 %>%
    mutate(perc_micro = round((100*micro)/non_tot, digits = 2),
           perc_nano = round((100*nano)/non_tot, digits = 2),
           perc_nanomicro = round((100*`nano-micro`)/non_tot, digits = 2),
           perc_pico = round((100*pico)/non_tot, digits = 2),
           perc_piconano = round((100*`pico-nano`)/non_tot, digits = 2),
           perc_nontotal = round((100*non_tot)/tot, digits = 2),
           perc_total = round((100*total)/tot, digits = 2)) -> final_size.df2
  

# ==== 5. Merging of data ==== 

  summary_tab <- cbind(count.asv.sp, 
                       count.reads.sp$ntot, 
                       final_depth.df$perc_surf, final_depth.df$perc_euph,
                       final_size.df$perc_micro, final_size.df$perc_nano, final_size.df$perc_nanomicro,
                       final_size.df$perc_pico, final_size.df$perc_piconano, final_size.df$perc_total)
  
  #add occupancy
  summary_tab.f <- left_join(summary_tab, evenness_sites_troph[,c(2,4,5,6)], by = ("species" = "species"))
  
  #rename columns 
  summary_tab.f %>%
    rename(n_reads = "count.reads.sp$ntot",
           perc_surf = "final_depth.df$perc_surf",
           perc_euph = "final_depth.df$perc_euph",
           perc_micro = "final_size.df$perc_micro", 
           perc_nano = "final_size.df$perc_nano", 
           perc_nanomicro = "final_size.df$perc_nanomicro", 
           perc_pico = "final_size.df$perc_pico", 
           perc_piconano = "final_size.df$perc_piconano", 
           perc_total = "final_size.df$perc_total") -> summary_tab.f
  
  
  #We remove the X-species from the dataset
  
  summary_tab.f.X <- summary_tab.f.X[!grepl("X",summary_tab.f$species),]
  
  
  # Redo it with the second df 
  
  summary_tab2 <- cbind(count.asv.sp, 
                       count.reads.sp$ntot, 
                       final_depth.df$perc_surf, final_depth.df$perc_euph,
                       final_size.df2$perc_micro, final_size.df2$perc_nano, final_size.df2$perc_nanomicro,
                       final_size.df2$perc_pico, final_size.df2$perc_piconano, final_size.df2$perc_total, 
                       final_size.df2$perc_nontotal)
  
  #add occupancy
  summary_tab.f2 <- left_join(summary_tab2, evenness_sites_troph[,c(2,4,5,6)], by = ("species" = "species"))
  
  #rename columns 
  summary_tab.f2 %>%
    rename(n_reads = "count.reads.sp$ntot",
           perc_surf = "final_depth.df$perc_surf",
           perc_euph = "final_depth.df$perc_euph",
           perc_micro = "final_size.df2$perc_micro", 
           perc_nano = "final_size.df2$perc_nano", 
           perc_nanomicro = "final_size.df2$perc_nanomicro", 
           perc_pico = "final_size.df2$perc_pico", 
           perc_piconano = "final_size.df2$perc_piconano", 
           perc_total = "final_size.df2$perc_total",
           perc_nontotal = "final_size.df2$perc_nontotal") -> summary_tab.f2
  
  
  #We remove the X-species from the dataset
  
  summary_tab.f2.X <- summary_tab.f2[!grepl("X",summary_tab.f2$species),]
  



# ==== 6. Plots =====


    # ASV
    
      summary_tab.f2.X$alpha <- rank(tolower(summary_tab.f2.X$species))
      summary_tab.f2.X$species <-reorder(summary_tab.f2.X$species, -summary_tab.f2.X$alpha)
      
      
      asv_plot <- ggplot(summary_tab.f2.X, aes(x = n_asv, y = species)) +
                         geom_bar(width = 0.5, stat = "identity") + 
                         theme_classic() +
                         theme(text = element_text(size = 10))
      
      asv_plot <- ggplot(summary_tab.f2.X, aes(x = n_asv, y = species)) +
        # geom_col() + 
        geom_bar(width = 0.5, stat = "identity") + 
        theme_classic() + theme(text = element_text(size = 10)) +
        labs(x = "b. ASVs")
    
    
    
    #Reads 
    
      reads_plot <- ggplot(summary_tab.f2.X, aes(x = n_reads, y = species)) + 
                        geom_bar(width = 0.5, stat = "identity") + 
                        theme_classic() + theme(text = element_text(size = 10)) +
        labs(x = "a. Reads")
    
    #Depth
    
      #we first need to change the data shape
      
      depth_tall <- summary_tab.f2.X[,c(1,4,5)]
      
      depth_tall %>%
        gather(key = depth_level, value = perc, perc_surf:perc_euph) -> depth_tall
      
      depth_plot <- ggplot(depth_tall, aes(x = species, y = perc , fill = depth_level )) + 
        geom_bar(stat = "identity", width = 1) + 
        coord_flip() +
        theme_classic() + theme(text = element_text(size = 10), legend.position = "none") +
        scale_fill_brewer(palette = "Set1") +
        labs(y = "c. Depth")
      
        
    
    #Size 
    
      #Like for depth, we change table shape 
      
      size_tall <- summary_tab.f.X[,c(1,6,7,8,9,10,11)]
      
      size_tall %>%
        gather(key = fraction_name, value = perc, perc_micro:perc_total) -> size_tall
      
      size_plot <- ggplot(size_tall, aes(x = species, y = perc , fill = fraction_name )) + 
        geom_bar(stat = "identity") + 
        coord_flip() +
        theme_classic() + theme(text = element_text(size = 10)) 
      
      # we redo the size plots with the second df 
      
      size_fraction_non <- summary_tab.f2.X[,c(1,6,7,8,9,10)]
      size_fraction_total <- summary_tab.f2.X[,c(1,11,12)]
      
      size_fraction_non %>%
        gather(key = fraction_name, value = perc, perc_micro:perc_piconano) -> size_fraction_non
      
      size_fraction_non$fraction_name <- as.factor(size_fraction_non$fraction_name)
      levels(size_fraction_non$fraction_name) <- c("micro", "nano", "nano_micro", "pico", "pico_nano")
      
      size_plot_non <- ggplot(size_fraction_non, aes(x = species, y = perc , fill = fraction_name )) + 
        geom_bar(stat = "identity", width = 1) + 
        coord_flip() +
        theme_classic() +
        scale_fill_brewer(palette = "Dark2", limits = c("pico", "pico_nano", "nano", "nano_micro", "micro")) +
        theme(legend.position="none", text = element_text(size = 10)) +
        labs(y = "e. Size fraction")
      
      size_fraction_total %>%
        gather(key = fraction_name, value = perc, perc_total:perc_nontotal) -> size_fraction_total
      
      palette <- brewer.pal(n = 7, name = "Set1")
      
      size_plot_tot <- ggplot(size_fraction_total, aes(x = species, y = perc , fill = fraction_name)) + 
        geom_bar(stat = "identity", width = 1) + 
        coord_flip()  +
        theme_classic() +
        theme(legend.position="none", text = element_text(size = 10)) +
        scale_fill_manual(values = palette[c(3,5)]) +
        labs(y = "d. Size (Total/Non-Total)")
    
    #Occupancy
    
      occ_plot <- ggplot(summary_tab.f2.X, aes(x = occupancy, y = species)) + 
                    geom_bar(width = 0.5, stat = "identity") + 
                    theme_classic() +
                    theme(legend.position="none", text = element_text(size = 10)) +
                    labs(x = "f. Occupancy")
      
      
      # Now we combine them all
      
      ggarrange(reads_plot,
                asv_plot + theme(axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank() ),
                depth_plot + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank() ),
                size_plot_tot + theme(axis.text.y = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.title.y = element_blank() ),
                size_plot_non + theme(axis.text.y = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.title.y = element_blank() ),
                occ_plot + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank() ),
                nrow = 1,
                widths = c(1,1,1,1,1,1)) -> sum_plot
      
      
      ggsave("suppmat/sum_plot.pdf", sum_plot, width = 15, height = 25)
    
      # 2-parts plots 
      
        # split the dataset in two parts
        summary_tab.f2.X_1 <- summary_tab.f2.X[1:129,]
        summary_tab.f2.X_2 <- summary_tab.f2.X[130:258,]
        
          # first dataset
          col_labels_1 <- paste0(summary_tab.f2.X_1$Trophy, summary_tab.f2.X_1$typeMX)
          col_labels_1 <-  gsub("mixoplanktonCM", "firebrick1", col_labels_1)
          col_labels_1 <-  gsub("mixoplanktoneSNCM", "purple1", col_labels_1)
          col_labels_1 <-  gsub("mixoplanktonpSNCM", "sienna3", col_labels_1)
          col_labels_1 <-  gsub("phytoplanktonNA", "olivedrab3", col_labels_1)
          col_labels_1 <-  gsub("protozooplanktonNA", "dodgerblue2", col_labels_1)
          col_labels_1 <- ifelse(is.na(col_labels_1), "grey30", gsub("NANA", "grey30", col_labels_1))
          
          col_labels_1.rev <- rev(col_labels_1)
          
          # Second dataset
          col_labels_2 <- paste0(summary_tab.f2.X_2$Trophy, summary_tab.f2.X_2$typeMX)
          col_labels_2 <-  gsub("mixoplanktonCM", "firebrick1", col_labels_2)
          col_labels_2 <-  gsub("mixoplanktoneSNCM", "purple1", col_labels_2)
          col_labels_2 <-  gsub("mixoplanktonpSNCM", "sienna3", col_labels_2)
          col_labels_2 <-  gsub("phytoplanktonNA", "olivedrab3", col_labels_2)
          col_labels_2 <-  gsub("protozooplanktonNA", "dodgerblue2", col_labels_2)
          col_labels_2 <- ifelse(is.na(col_labels_2), "grey30", gsub("NANA", "grey30", col_labels_2))
          
          col_labels_2.rev <- rev(col_labels_2)
          
        
        # Reads
        reads_plot_1 <- ggplot(summary_tab.f2.X_1, aes(x = n_reads, y = species)) + 
          geom_bar(width = 0.5, stat = "identity") + 
          theme_classic() + theme(axis.title.x = element_text(size = 10),
                                  axis.text.y = element_text(size = 8)) +
          labs(x = "Reads") +
          theme(axis.text.y = element_text(colour = col_labels_1.rev),
                axis.text.x = element_text(size = 6))
        
        reads_plot_2 <- ggplot(summary_tab.f2.X_2, aes(x = n_reads, y = species)) + 
          geom_bar(width = 0.5, stat = "identity") + 
          theme_classic() + theme(text = element_text(size = 10),
                                  axis.text.y = element_text(size = 8)) +
          labs(x = "Reads") +
          theme(axis.text.y = element_text(colour = col_labels_2.rev),
                axis.text.x = element_text(size = 6))
        
        
        
        
      
        # ASVs 
        asv_plot_1 <- ggplot(summary_tab.f2.X_1, aes(x = n_asv, y = species)) +
          # geom_col() + 
          geom_bar(width = 0.5, stat = "identity") + 
          theme_classic() + theme(text = element_text(size = 10)) +
          labs(x = "ASVs") +
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 6))
        
        asv_plot_2 <- ggplot(summary_tab.f2.X_2, aes(x = n_asv, y = species)) +
          # geom_col() + 
          geom_bar(width = 0.5, stat = "identity") + 
          theme_classic() + theme(text = element_text(size = 10)) +
          labs(x = "ASVs") +
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 6))
    
            
          # Occupancy 
          
          occ_plot_1 <- ggplot(summary_tab.f2.X_1, aes(x = occupancy, y = species)) + 
            geom_bar(width = 0.5, stat = "identity") + 
            theme_classic() + 
            theme(text = element_text(size = 10)) +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = 6))
            
          
          occ_plot_2 <- ggplot(summary_tab.f2.X_2, aes(x = occupancy, y = species)) + 
            geom_bar(width = 0.5, stat = "identity") + 
            theme_classic() + 
            theme(text = element_text(size = 10)) +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = 6))
          
          # Depth 
          
          depth_tall_1 <- depth_tall[which(depth_tall$species %in% summary_tab.f2.X_1$species),]
          
          depth_plot_1 <- ggplot(depth_tall_1, aes(x = species, y = perc , fill = depth_level )) + 
            geom_bar(stat = "identity", width = 1) + 
            coord_flip() +
            theme_classic() + theme(text = element_text(size = 10), legend.position = "none") +
            scale_fill_brewer(palette = "Set1") +
            labs(y = "Depth") +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = 6))
          
          depth_tall_2 <- depth_tall[which(depth_tall$species %in% summary_tab.f2.X_2$species),]
          
          depth_plot_2 <- ggplot(depth_tall_2, aes(x = species, y = perc , fill = depth_level )) + 
            geom_bar(stat = "identity", width = 1) + 
            coord_flip() +
            theme_classic() + theme(text = element_text(size = 10), legend.position = "none") +
            scale_fill_brewer(palette = "Set1") +
            labs(y = "Depth") +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = 6))
          
          
          # Size fraction
          
            # Tot 
          
            size_fraction_total_1 <- size_fraction_total[which(size_fraction_total$species %in% summary_tab.f2.X_1$species),]
            
            size_plot_tot_1 <- ggplot(size_fraction_total_1, aes(x = species, y = perc , fill = fraction_name)) + 
              geom_bar(stat = "identity", width = 1) + 
              coord_flip()  +
              theme_classic() +
              theme(legend.position="none", text = element_text(size = 10)) +
              scale_fill_manual(values = palette[c(3,5)]) +
              labs(y = "Size (Total/Non-Total)") +
              theme(axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(size = 6))
            
            size_fraction_total_2 <- size_fraction_total[which(size_fraction_total$species %in% summary_tab.f2.X_2$species),]
            
            size_plot_tot_2 <- ggplot(size_fraction_total_2, aes(x = species, y = perc , fill = fraction_name)) + 
              geom_bar(stat = "identity", width = 1) + 
              coord_flip()  +
              theme_classic() +
              theme(legend.position="none", text = element_text(size = 10)) +
              scale_fill_manual(values = palette[c(3,5)]) +
              labs(y = "Size (Total/Non-Total)") +
              theme(axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(size = 6))
            
            # Non
            
            size_fraction_non_1 <- size_fraction_non[which(size_fraction_non$species %in% summary_tab.f2.X_1$species),]
            
            size_plot_non_1 <- ggplot(size_fraction_non_1, aes(x = species, y = perc , fill = fraction_name )) + 
              geom_bar(stat = "identity", width = 1) + 
              coord_flip() +
              theme_classic() +
              scale_fill_brewer(palette = "Dark2", limits = c("pico", "pico_nano", "nano", "nano_micro", "micro")) +
              theme(legend.position="none", text = element_text(size = 10)) +
              labs(y = "Size fraction") +
              theme(axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(size = 6))
            
            size_fraction_non_2 <- size_fraction_non[which(size_fraction_non$species %in% summary_tab.f2.X_2$species),]
            
            size_plot_non_2 <- ggplot(size_fraction_non_2, aes(x = species, y = perc , fill = fraction_name )) + 
              geom_bar(stat = "identity", width = 1) + 
              coord_flip() +
              theme_classic() +
              scale_fill_brewer(palette = "Dark2", limits = c("pico", "pico_nano", "nano", "nano_micro", "micro")) +
              theme(legend.position="none", text = element_text(size = 10)) +
              labs(y = "Size fraction") +
              theme(axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(size = 6))
            
    
          # Combining 
          
          
          sum1 <- ggarrange(reads_plot_1,
                            asv_plot_1,
                            depth_plot_1,
                            size_plot_tot_1,
                            size_plot_non_1,
                            occ_plot_1,
                            nrow = 1,
                            widths = c(1.5,1,1,1,1,1))
         
           sum2 <- ggarrange(reads_plot_2,
                            asv_plot_2,
                            depth_plot_2,
                            size_plot_tot_2,
                            size_plot_non_2,
                            occ_plot_2,
                            nrow = 1,
                            widths = c(1.5,1,1,1,1,1))
      
          
          ggsave("suppmat/Sum1.pdf", sum1, width = 20, height = 15)
          ggsave("suppmat/Sum2.pdf", sum2, width = 20, height = 15)
          
      
      
      
      
      # We get the legend for colored columns
      
      depth_plot.leg <- ggplot(depth_tall, aes(x = species, y = perc , fill = depth_level )) + 
        geom_bar(stat = "identity", width = 1) + 
        coord_flip() +
        theme_classic() + theme(text = element_text(size = 6), legend.position = "bottom") +
        scale_fill_brewer(palette = "Set1", labels = c("Euphotic", "Surface")) +
        labs(fill = "Depth") +
        theme_void()
      
      
      size_plot_non.leg <- ggplot(size_fraction_non, aes(x = species, y = perc , fill = fraction_name )) + 
        geom_bar(stat = "identity", width = 1) + 
        coord_flip() +
        theme_classic() +
        scale_fill_brewer(palette = "Dark2", limits = c("pico", "pico_nano", "nano", "nano_micro", "micro")) +
        theme(legend.position="bottom", text = element_text(size = 6)) +
        labs(fill = "Size Fraction") +
        theme_void()
      
      
      
      size_plot_tot.leg <- ggplot(size_fraction_total, aes(x = species, y = perc , fill = fraction_name)) + 
        geom_bar(stat = "identity", width = 1) + 
        coord_flip()  +
        theme_classic() +
        theme(legend.position="bottom", text = element_text(size = 6)) +
        scale_fill_manual(values = palette[c(3,5)], labels = c("Non-total", "Total")) +
        labs(fill = "Size Type (Total / Non-Total)") +
        theme_void()
      
      ggsave("suppmat/leg1.pdf", depth_plot.leg)
      ggsave("suppmat/leg2.pdf", size_plot_non.leg)
      ggsave("suppmat/leg3.pdf", size_plot_tot.leg)
      
      
    

# ========================================================================================






