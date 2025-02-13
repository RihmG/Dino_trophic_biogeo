
##  ------- TROPHIC ANNOTATION / ABUNDANCE MATRIX BUILDING ----------

# This script takes raw data extracted from metaPR2 (abundance of ASVs of dinoflagellates species worldwide)
# and trophic annotation from Schneider et al. 2020 and creates abundance matrix corresponding to species in 
# sites of interest.

# Note: The Update part of the script correspond to the addition of trophic information from Mitra's 
# mixotrophy DataBase (2023)

# Libraries 
library(tidyverse)
library(funrar)
library(fossil)
library(vegan)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(raster)


#------ 1. Data Importation  -------

  #metaPR2 data, downloaded 25/11/2022
  
  data2511<- read.table(file = 'metapr2_ASVs_selected_abundance_Dinophyceae_2022-11-25.tsv',
                        sep = '\t', header = TRUE)
  
  #Schneider et al. (2020) database about trophic modes for dinoflagellates
  trophic_modes<-read.csv(file = 'trophic_mode.csv')
    
  
  #------ Counting species/genus in dataset ----------
  
  #Number of species in species column  
  data2511 %>% 
    group_by(species) %>%
    summarise(no_rows = length(species)) -> data2511_sp #df with species 
  
  #Number of genera in genus column 
  data2511 %>% 
    group_by(genus) %>%
    summarise(no_rows = length(genus)) -> data2511_gn #df with genera 



#------ 2. First Assignment  -------

  #before Joining, checking of the text format
  
  if(length(pat2 <- grep("_sp\\.", data2511_sp$species))) #looking at the rank where there is the "_sp." pattern
    # the \\ turns the . into a regular character
    pat2# pattern rank in the dataframe 
  
  pattern2 <- data2511_sp$species[pat2] #rank/species name correspondence 
  pattern2 <- as.data.frame(pattern2) #
  
  #Building of jointing dataset
  trophic_join <- trophic_modes[,c(1:4)]
  trophic_join <- trophic_join[,-2]
  
  
  #First joint on the species field basis 
  
  join1 <-left_join(data2511, trophic_join, by = c("species" = "ScientificName"))
  
  #How much NA after the first joint ? 
  join1 %>% 
    group_by(Trophy) %>%
    summarise(no_rows = length(Trophy)) # --> 94 403 NA observations 
  
  #Now we assign a trophic mode on the Genus field basis 
  #However, in the Schneider 2020 dataset, some genera belong to different trophic modes
  #We first reduce the dataset like this : 
  
  #  - if the genus match with one unique trophic mode, this mode is assigned to the genus
  #  - if the genus match with at least two different trophic modes, the genus has no matching trophic mode
  #(i.e not used for the following analysis)
  
  trophic_join2<-cbind(trophic_join, trophic_modes$Genus)
  
  trophic_gen <- trophic_join2[,2:4] #we first create a df with only genus 
  
  #renaming the genus column 
  trophic_gen %>%
    rename(genus = `trophic_modes$Genus`) -> trophic_gen
  
  # annotation ------
  
      #initialisation
      levels(as.factor(trophic_gen$genus))->trophgenus #vector with all possible genera 
      levels(as.factor(trophic_gen$Trophy))->trophTrophy #same for trophic modes 
      levels(as.factor(trophic_gen$typeMX))->trophMX #same for types of mixotrophy 
      
      trophgenus
      trophTrophy
      trophMX
      
      vectorposition=c()
      repetitionNumberGenus=c()
      FirstMode=c()
      newvectorMode=c()
      GenusMode=c()
      tpMX=c()#We store for each for loop iteration the mixotrophy type value corresponding to the genus and
      #from Schneider et al. 2020 dataset
      
      for (i in 1:length(trophgenus)){
        vectorposition=which(trophic_gen$genus==trophgenus[i])# position in trophic_gen of genus i in trophgenus
        repetitionNumberGenus=length(vectorposition)#the length of vectorposition allows to know the number of repetition of for loop for one genus
        FirstMode=trophic_gen$Trophy[vectorposition[1]]#identifying the first trophic type for each genus  
        GenusMode=1#index which tell if there is only one trophic type per genus (=1) or several (=0)
        for (j in 1:repetitionNumberGenus){ #for each genus, we read all the lines and so all possible trophic types
          GenusMode=as.numeric(FirstMode==as.factor(trophic_gen$Trophy[vectorposition[j]]))*GenusMode#assigning 1 to GenusMode while trophic types for a same genus are 
        }# the same as the first trophic type met by genus; otherwise, GenusMode becomes 0
        if (GenusMode==1){newvectorMode=c(newvectorMode,FirstMode)#if GenusMode = 1, trophic types are the same for a same genus, so we stock in the new vector the initial trophic type
        tpMX=c(tpMX,trophic_gen$typeMX[vectorposition[j]])}
        else {newvectorMode=c(newvectorMode,NA)#if GenusMode = 0, there are several trophic type per genus, then we stock NA in the new vector
        tpMX=c(tpMX,NA)}#same for type of mixotrophy
      }
      
      new_trophic_gen=matrix(c(trophgenus,newvectorMode,tpMX),length(trophgenus),3)
      new_trophic_gen<-as.data.frame(new_trophic_gen)#matrix with only one trophic type per genus 
      
  
  #================
  #we test if mixotrophic have homogeneous types of mixotrophy by genus 
  test.troph<-trophic_modes[trophic_modes$Trophy=='mixoplankton',]
  test.troph<-test.troph[,c(3,4,12)]
  #only Gymnodinium could be CM or pSNCM, but anyway in the new_trophic_gen it has no trophic type cause it could be
  #mixo, hetero or auto depending on the species
  #================
  
  
  #renaming columns 
  names(new_trophic_gen)[1] <- "genres"
  names(new_trophic_gen)[2] <- "Trophy"
  names(new_trophic_gen)[3] <- "typeMX"
  
  saveRDS(new_trophic_gen, "data_out/schneider_genus.RDS")
  
  # --> we now have a reliable dataset for assigning trophic types by genera on MetPR2 dataset 
  
  #on réalise une jointure des données MetaPR2 (lignes non assignées par espèces = join_gen) avec new_trophic_gen



## ------ 3. Second Assignment --------

  #we extract all the non-assigned lines from the first assignment based on species(=NA lines)
  join_gen2<-subset(join1, is.na(join1[,c("Trophy")]))
  
  #we make the assignment on the genus field 
  
  join_gen2<-join_gen2[,-c(48,49)]#removing the two Trophy and typeMX columns
  
  join2<-left_join(join_gen2, new_trophic_gen, by = c("genus" = "genres"))
  
  join2 %>% 
    group_by(Trophy) %>%
    summarise(no_rows = length(Trophy)) 

#==================

    #Creation of dframes with Trophic information 
  
  #from the first assignment
  sub_sp2<-subset(join1, !is.na(join1[,48]))
  
  #from the 2nd assignment
  sub_gn2<-subset(join2, !is.na(join2[,48]))
  #==================
  
  #Dataframe with missing trophic information after the two assignments
  
  sub_non2<-subset(join2, is.na(join2[,48]))
  
  ## ------ 4. Assignment case-by-case -------
  
  #Non-assigned species 
  sub_non2 %>% 
    group_by(species) %>%
    summarise(no_rows = length(species)) -> sub_non_sp2
  
  #trying to find "/" patterns 
  where_is_slash2<-grep("/", data2511_sp$species)
  where_is_slash2
  data2511_sp$species[where_is_slash] #Heterocapsa_nei/rotundata
  
  #Heterocapsa assignment
  sub_non2$Trophy = ifelse(sub_non2$species == "Heterocapsa_nei/rotundata", "mixoplankton", NA)
  sub_non2$typeMX = ifelse(sub_non2$species == "Heterocapsa_nei/rotundata", "CM", NA)
  
  sub_last1 <- subset(sub_non2, !is.na(sub_non2[,48]))
  
  sub_non3 <- subset(sub_non2, is.na(sub_non2[,48]))
  
  
  #Warnowaiceae are mostly heterotrophic 
  # ("Molecular phylogeny of ocelloid-bearing dinoflagellates (Warnowiaceae) as inferred from SSU and LSU rDNA sequences")
  # Hoppenrath et al. (2009)
  
  
  sub_non3$Trophy = ifelse(sub_non3$family == "Warnowiaceae", "protozooplankton", NA)
  
  sub_last2<-subset(sub_non3, !is.na(sub_non3[,48]))
  
  sub_non4 <- subset(sub_non3, is.na(sub_non3[,48]))
  

  # FINALLY we merge the assigned datasets 
  metatroph <- rbind(sub_last1, sub_last2, sub_gn2, sub_sp2)
  
  #Then we merge with all the remaining not-attributed data 
  datatroph <- rbind(metatroph, sub_non4)
  
  trophic.sp <- datatroph[,c(45,48,49)]
  trophic.sp <- trophic.sp[!duplicated(trophic.sp), ]
  summary(as.factor(trophic.sp$Trophy))


# ------- 5. Depth influence ---------

  # Check if there is potential effect of depth in community structure 
  
  sub_surf<-data2511[data2511$depth_level == "surface",]
  sub_euph<-data2511[data2511$depth_level == "euphotic",]
  
  ggplot(data2511, aes(x = longitude, y = depth, color = as.factor(project)))+
    geom_point() +
    scale_y_reverse()
  

  #map of spatial distribution of samples
  worldmap <- ne_countries(scale = "medium", returnclass = "sf")
  
  #on va extraire les localisations des samples 
  
  loc2<-data2511[,c("file_code", "longitude","latitude" )]
  
  loc2 %>%
    group_by(file_code) %>%
    slice(1) -> loc2
  
  
  coords2<-cbind(loc2$longitude, loc2$latitude)
  sploc2<-SpatialPointsDataFrame(coords2, loc2) #"sp" object
  proj4string(sploc2) <- CRS("+init=epsg:4326")
  
  #turn to sf object
  sfloc2<-st_as_sf(
    loc2,
    coords = c("longitude","latitude"),
    remove = TRUE,
    na.fail = TRUE,
    sf_column_name = NULL, 
    crs = 4326
  )
  
  #plot
  ggplot()+
    geom_sf(data = worldmap, color = "grey", fill = "grey") +
    theme(panel.background = element_rect(fill = "aliceblue")) +
    geom_point(data = subset(data2511, depth_level == 'euphotic'), aes(x = longitude, y = latitude, color = depth_level))
  
  ggplot()+
    geom_sf(data = worldmap, color = "grey", fill = "grey") +
    theme(panel.background = element_rect(fill = "aliceblue")) +
    geom_point(data = subset(data2511, depth_level == 'surface'), aes(x = longitude, y = latitude))

  #we remove samples which are below MLD 
    
  #import MLD layer  
  mld<-read.csv("Env/woa18_A5B7_M00an04.csv", header = F, sep = ",", quote = "\"",
                dec = ".")
  
  colnames(mld) <- c("latitude", "longitude", "value")
  mld <- mld[-c(1,2),]
  
  mld$latitude <-as.numeric(mld$latitude)
  mld$longitude <- as.numeric(mld$longitude)
  mld$value <- as.numeric(mld$value)
  
  
  #turn .csv into raster layer
  
  mld.sp2<-mld
  coordinates(mld.sp2) <- ~ longitude + latitude 
  
  ras_dom<-raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
                  crs="+proj=longlat +datum=WGS84 +no_defs ",
                  resolution = c(1/4, 1/4), vals=NA)
  
  mld.rast <- rasterize(mld.sp2, ras_dom, "value", update = TRUE) 
  
  
  #now we want to extract values of MLD raster from sites location 
  
  #we create a sp object with sample location 
  samples.spatial <- data2511[,c(1,20,21)]
  coordinates(samples.spatial) <- ~ longitude + latitude 
  
  #we extract the mld value for the ASV/samples dataset
  samples.mld<-as.data.frame(raster::extract(mld.rast, samples.spatial))
  samples.mld %>%
    rename(mld = `raster::extract(mld.rast, samples.spatial)`) -> samples.mld
  
  #merge the data
  data0401<-data2511
  data0401<-cbind(data0401, samples.mld)
  
  
  #now we want to remove all the mld values with NA, cause it's too close from the coast
  data0401.1<-data0401[!is.na(data0401$mld),]

  #now we remove all the rows with depth > mld 
  
    #first we assign a 0 to each NA depth value
    data0401.1[is.na(data0401.1)] <- 0
    
    data0401.2<-subset(data0401.1, depth < mld)
    
    #write_tsv(data0401.2, file = "data0401_2.tsv")
    write_csv(data0401.2, file = "data0401_2.csv")
    
    data0401.2 %>%
      group_by(species) %>%
      summarise()-> zzz

    

# ------ UPDATE 06/03 -------- 

    # New mixotrophic database 
    
    mdb <- read.csv2("dino_mdb.csv")
    
    # removing "*" in the end of TypeMX
    mdb$typeMX <-  gsub("\\*$", "", mdb$typeMX)
    
    # adding "_" in species names
    mdb$Species.Name <- gsub(" ", "_", mdb$Species.Name)
    
    # Now we compare with Schneider 2020, but only with Dinophyceae 
    schneider <- trophic_modes[trophic_modes$Class == 'Dinophyceae',]
    #schneider_mixo <- schneider[schneider$Trophy == 'mixoplankton',]
    
    mdb <- mdb[,c(1,3,4)]
    schneider <- schneider[,c(1,3,4)]
    
    mdb$Species.Name %in% schneider$ScientificName
    
    # we now spot species in schneider which are present in mdb 
    schneider[which(schneider$ScientificName %in% mdb$Species.Name),]
    
    compare <- left_join(schneider, mdb, by = c("ScientificName" = "Species.Name"))
    
    compare %>%
      rename(
        TrophyS = Trophy.x,
        typeMXS = typeMX.x,     
        TrophyM = Trophy.y,
        typeMXM = typeMX.y
      ) -> compare
    
    compare <- compare[!is.na(compare$ScientificName),]
    
    # we spot species present in mdb but not in schneider 
    mdb_only <- mdb[!mdb$Species.Name %in% schneider$ScientificName,]
    
    schneider$ScientificName %in% mdb_only$Species.Name 
    
    # --> looks like no species from mdb not present in schneider are in common with our database from metapr2
    
    # we will now create a new dataset based on mdb corrections 
    
    mdb_correction <- function(i){
      row = compare[i,]
      if (is.na(row$TrophyM))
      {
        row$trophy = row$TrophyS
        row$typemx = row$typeMXS
      }else {
        row$trophy = row$TrophyM
        row$typemx = row$typeMXM
      }
      return(row)
    }
    
    mdb_correction.list <- lapply(1:nrow(compare), mdb_correction)
    mdb_correction.df <- do.call(rbind, mdb_correction.list)
    
    
    
    # we remove spaces in order to compare 
    mdb_correction.df$TrophyS <- gsub("\\s+", "", mdb_correction.df$TrophyS)
    mdb_correction.df$TrophyM <- gsub("\\s+", "", mdb_correction.df$TrophyM)
    
    # here are the species which trophic type changed between schneider and mdb 
    mdb_correction.df[which(mdb_correction.df$TrophyS != mdb_correction.df$TrophyM),]
    
    
    # now we try to spot species in mdb not in schneider and potentially in metapr2 dataset
    new_mixo <- mdb_only[which(mdb_only$Species.Name %in% evenness_sites_troph$species),]
    
    saveRDS(new_mixo, "data_out/new_mixo.RDS")
    

# ===============================================================================================






