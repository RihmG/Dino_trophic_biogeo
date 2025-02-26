
####### ------- BILANSITES --------

# This script takes in input a dataframe summarizing all the samples present in one unique coordinates site
# and creates a new abundance matrix averaging the abundance of reads for each species according to the number of 
# samples within each unique site


# Libraries 
library(dplyr)
library(tidyverse)
library(fossil)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)


# Reducing dataset with sites = unique Latitude Longitude couple

#First, we find the locations of unique coordinates couples, then we attribute the number of samples per location


#-------- 1. Data  -------
  
    
  sites2<-read.table("Lat_Long_Proj_ProjD_NbS_List2_2023_2", sep = "\t", h=T)# from a .sh script (-->README)
  sites2<-sites2[-1,]
  rownames(sites2) <- 1:nrow(sites2)
  
  #data0401.2 <- read_tsv("data0401_2.tsv")
  data0401.2 <- read_csv("data0401_2.csv")
  
  matsp2 <-fossil::create.matrix(as.data.frame(data0401.2), 
                        tax.name = "species",
                        locality = "file_code", 
                        abund.col = "n_reads",
                        abund = TRUE)
  
  matsp2<-as.data.frame(matsp2)# matrix transposition for next functions

  #mapping with number of samples per location 
  worldmap <- ne_countries(scale = "medium", returnclass = "sf")
  
  ggplot()+
    geom_sf(data = worldmap, color = "grey", fill = "grey") +
    theme(panel.background = element_rect(fill = "aliceblue")) +
    geom_point(data = sites2, aes(x = longitude, y = latitude, color = Nbsamples)) +
    scale_color_viridis(direction = -1) #+
  #xlim(-80,-60) +
  #ylim(-70,-60) 

  ggplot()+
    geom_sf(data = worldmap, color = "grey", fill = "grey") +
    theme(panel.background = element_rect(fill = "aliceblue")) +
    geom_point(data = sites2, aes(x = longitude, y = latitude, size = Nbsamples), color = 'blue', pch = 21)
  
  #with projects 
  ggplot()+
    geom_sf(data = worldmap, color = "grey", fill = "grey") +
    theme(panel.background = element_rect(fill = "aliceblue")) +
    geom_point(data = sites2, aes(x = longitude, y = latitude, color = project))

  
# ----- 2. Annual Treatment ------  
  
  #creating ID for each sites   
  sites2$site_ID<-paste("site",as.numeric(rownames(sites2)), sep="")
  
  #we remove the "." in sample names 
  sites2$listSamples <- gsub("\\.", "-", sites2$listSamples)
  
  #function for taking each samples per site location in sites1 and give a sub-abundance matrix associated 
  selectsample2<-function(i){
    list1<-strsplit(sites2$listSamples[i], split  = "|", fixed = T)
    print(list1)
    collabel <- which(gsub("\\.", "-", colnames(matsp2)) %in% unlist(list1))
    matsp2 %>% dplyr::select(collabel)
  }
  
  selectsample2(3)
  
  #now we create a list of sub-community matrices per sample for each site 
  listmat2<-lapply(as.list(1:length(as.data.frame(t(sites2)))), selectsample2)
  
  
  #Check if the number of sample per matrices in listmat is the same as the number calculated in site1
  testlist2<-function(i){
    return(isTRUE(length(listmat2[[i]])==sites2$Nbsamples[i]))
  }
  
  testlistr2<-lapply(as.list(1:895), testlist2)
  #looks ok
  
  
  #we create a list of one matrix per site, by combining the different matrices of the different samples
  #In order to cope with multi-sampling for the same location, we average the abundance 
  
  bilansitesmoy2<-function(i){
    mat.i<-as.data.frame(listmat2[[i]])
    mat.i$site<-rowSums(mat.i)/sites2$Nbsamples[i]
    mat.i2<-mat.i[,-c(1:sites2$Nbsamples[i])]
    mat.i2<-as.data.frame(mat.i2)
    colnames(mat.i2)<-paste("site",i)
    rownames(mat.i2)<-rownames(matsp2)
    return(as.data.frame(mat.i2))
  }
  
  
  bilansitesmoy2(9)# Test
  
  bilansites_listmoy2<-lapply(as.list(1:895), bilansitesmoy2)
  
  #we create a final raw abundance matrix with all the sites and species  
  matfinalemoy2<-do.call(cbind, bilansites_listmoy2)
  
  saveRDS(matfinalemoy2, "data_out/mothermat.RDS")
  

  

# ---- 3. Month Treatment ------

  # We do the same as before, but we divide the samples depending on the month they were sampled

  key_month <- data0401.2[, c(1,3,15,20,21,45)]
  key_month$month <- sub(".*-(\\d+)-.*", "\\1", key_month$date)
  key_month$year <- substr(key_month$date, 1, 4)
  key_month$file_code <- gsub("\\.", "-", key_month$file_code)
  key_month$month_num <- as.numeric(key_month$month)

# We check if there is only one month associated with one sample  
  
  list_code <- unique(key_month$file_code) # 1883
  sink("file_code.txt") # save the output 

# 
  for (code in list_code) 
  { 
    sub_key <- key_month[which(key_month$file_code == code),] # We create sub dataframes per code
    if (length(unique(sub_key$month)) == 1){ # we check if there is more than one month associated with one file code in the sub dataframe
      cat(paste0("On est bon pour ", code, "\n"))
    } else {
      cat(paste0("vérifie ", code, "\n"))
    }
  }
  sink(type = "output") 
  
  # There is only one file_code corresponding to a sampling month
  # So we can create a dataframe with one file_code associated with one month 
  month_key <- key_month %>%
    distinct(file_code, .keep_all = TRUE)
  
  month_key <- month_key[,c(1,7)]

  # Now we associate the month with the sample in the previous file
  
  sites3 <- sites2
  
  sites3 <- sites3 %>%
    separate_rows(listSamples, sep = "\\|") %>%
    left_join(month_key, by = c("listSamples" = "file_code"))

  
  # We have now a file that associate the site to the month
  sites_month <- sites3[,c(1,2,7,8)]
  sites_month$month_num <- as.numeric(sites_month$month)
  sites_month <- na.omit(sites_month)
  sites_month <- distinct(sites_month) # 983 


  # We can now create monthly presence/absence datasets (abundance matrix, vectors)
  
  # We have now reference datasets  
  
  
  if(!dir.exists(paste0("month_occ/sites/")))
  {
    dir.create(paste0("month_occ/sites/"), recursive = T)
  }
  
  if(!dir.exists(paste0("month_occ/matrix/")))
  {
    dir.create(paste0("month_occ/matrix/"), recursive = T)
  }


  for (i in 1:12)
  {
    sites_month.i <- sites_month[sites_month$month_num == i,]
    sites_month.i <- left_join(sites_month.i, sites2[,c(5:7)], by = "site_ID")
    saveRDS(sites_month.i, paste0("month_occ/sites/month_occ_", i, ".RDS"))
  }



    # Now we create abundance matrices for each month, 
    
    for (m in 1:12)
    {
      
      site.i <- as.data.frame(readRDS(paste0("month_occ/sites/month_occ_", m, ".RDS"))) # each sites present 
      
      matsp.i <- as.data.frame(create.matrix(as.data.frame(key_month[key_month$month_num == m,]), tax.name = "species",
                                             locality = "file_code", 
                                             abund.col = "n_reads",
                                             abund = TRUE))
      
      #function for taking each samples per site location in site.i and give a sub-abundance matrix associated 
      selectsample.i<-function(i){
        list1<-strsplit(site.i$listSamples[i], split  = "|", fixed = T)
        #print(list1)
        collabel <- which(gsub("\\.", "-", colnames(matsp.i)) %in% unlist(list1))
        matsp.i %>% dplyr::select(collabel)
      }
      
      #now we create a list of sub-community matrices per sample for each site 
      listmat.i <-lapply(as.list(1:length(as.data.frame(t(site.i)))), selectsample.i)
      
      
      # We find the REAL number of samples that belong to the month m
      vec_Nbsamples <- c()
      for (i in 1:nrow(site.i))
      {
        list1<-strsplit(site.i$listSamples[i], split  = "|", fixed = T)
        #print(list1)
        newNbsamples <- length(which(gsub("\\.", "-", colnames(matsp.i)) %in% unlist(list1)))
        vec_Nbsamples <- c(vec_Nbsamples, newNbsamples)
      }
      
      site.i$NewNbsamples <- vec_Nbsamples
      
      #Check if the number of sample per matrices in listmat.i is the same as the number calculated in site.i
      testlist.i<-function(i){
        return(isTRUE(length(listmat.i[[i]])==site.i$NewNbsamples[i]))
      }
      
      testlistr.i<-lapply(as.list(1:nrow(site.i)), testlist.i)
      which(testlistr.i == FALSE) 
      
      
      #we create a list of one matrix per site, by combining the different matrices of the different samples
      #In order to cope with multi-sampling for the same location, we average the abundance 
      
      bilansitesmoy.i<-function(i){
        mat.i<-as.data.frame(listmat.i[[i]])
        mat.i$site<-rowSums(mat.i)/site.i$NewNbsamples[i]
        mat.i2<-mat.i[,-c(1:site.i$NewNbsamples[i])]
        mat.i2<-as.data.frame(mat.i2)
        #colnames(mat.i2)<-paste("site",i)
        colnames(mat.i2) <- site.i$site_ID[i]# we give the name of the corresponding absolute site to the column 
        rownames(mat.i2)<-rownames(matsp.i)
        return(as.data.frame(mat.i2))
      }
      
      bilansites_listmoy.i<-lapply(as.list(1:nrow(site.i)), bilansitesmoy.i)
      
      
      #we create a final raw abundance matrix with all the sites and species  
      matfinalemoy.i<-do.call(cbind, bilansites_listmoy.i)
      
      saveRDS(matfinalemoy.i, paste0("month_occ/matrix/mothermat_", m, ".RDS"))
      
      if(m == 12) { # We close the loop manually 
        break
      }
    }
    
      

# =======================================================================================

