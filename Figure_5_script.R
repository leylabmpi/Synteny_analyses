#requirements
library(tidyverse)
library(ggplot2)
library(effsize)
library(tmaptools)

#################################
# DESCRIPTION:  
# This is a file that concentrates the analyses of the Hubif metagenomes, for the syntracker paper
# The analysis is conducted as follows:
#
#   1. data import
#        a. metadata is imported (two files)
#        b. Syntracker output is imported and combined (this is a list, each element contains the SynTrakcer analysis of a single species)
#        Note: SynTRACKER uses a metadata file, however, I had to add a new metadata file later: Here I combined the 
#               new metadata with the final dataframe, holding the APSS values, etc. Thereofore, some of the actions are not conventional. 

#   2. correcting/adding to the metadata already contained in the SynTracker output
#   3. Sub selecting N regions/pairwise-comparison, calculating APSS. Converting the list into a SINGLE dataframe.





#################
#    input data
#################
cophylogeny_a<-readRDS(file= "/ebio/abt3_projects/Strain_tracking_synteny_blocks/envs/R_env/Cophylogeny/cophylogeny_list.rds")
cophylogeny_b<-readRDS(file= "/ebio/abt3_projects/Strain_tracking_synteny_blocks/envs/R_env/Cophylogeny/cophylogeny_srgs_10_to_40_MAGS.rds")
cophylogeny_c<-readRDS(file= "/ebio/abt3_projects/Strain_tracking_synteny_blocks/envs/R_env/Cophylogeny/cophylogeny_2_10.rds")

metadata<-read.csv(file="/ebio/abt3_projects/Strain_tracking_synteny_blocks/analysis/cophylogeny/cophylogeny_metadata_mod.csv", sep=";", header = TRUE)
new_metadata_Liam_currated<-read.table(file = "/ebio/abt3/henav/combined_LFCurated_01.07.2020.txt", sep="\t", header = TRUE)
names(cophylogeny_a)<-pathnames
names(cophylogeny_b)<-pathnames_10_40
names(cophylogeny_c)<-pathnames_2_10
cophy_out<-c(cophylogeny_a,cophylogeny_b, cophylogeny_c)

####################
# functions
####################


#specifiy Mother/infant 
add_groupE_mod<-function(ungroupdf, names, GroupE) {
  print(names)
  ungroupdf$GroupE1<-unlist(GroupE[ungroupdf$sample1])
  ungroupdf$GroupE2<-unlist(GroupE[ungroupdf$sample2])
  ungroupdf$is.same.GroupE<-ungroupdf$GroupE1==ungroupdf$GroupE2
  return(ungroupdf)
}

  
  # function to do the sub sampling of regions/pairwise comparison
  subsample_regions_mod<-function(big_organized_dfs) {
    if nrow(big_organized_dfs)
    set.seed(43)
    biggest_group<-max(big_organized_dfs %>% 
                         group_by(sample1, sample2, GroupA1, GroupA2, GroupB1, GroupB2, GroupC1, GroupC2, GroupD1, GroupD2, GroupE1, GroupE2, is.same.GroupA, is.same.GroupB, is.same.GroupC, is.same.GroupD, is.same.GroupE) %>% 
                         mutate(regions = n()) %>%
                         ungroup() %>%
                         pull(regions))
    ifelse(biggest_group >= 80, 
           newdf<-big_organized_dfs %>% 
             group_by(sample1, sample2, GroupA1, GroupA2, GroupB1, GroupB2, GroupC1, GroupC2, GroupD1, GroupD2, GroupE1, GroupE2, is.same.GroupA, is.same.GroupB, is.same.GroupC, is.same.GroupD, is.same.GroupE) %>% 
             filter(n() > 79) %>%
             sample_n(80) %>%
             mutate(regions = n()) %>%
             summarise(average_score=mean(syn_score), compared_regions=mean(regions)), 
           newdf<-"remove this one")
    return(newdf)
  }
  
  #function to add the SRG name to a column in the dataframe for this SRG (each SRG is a list element)
  add_names_to_col<-function(subs2, names) { #function to add element names as column in the df
    subs2<-subs2 %>% mutate(Species = names)
  }
  
  
  #calculate qvals, for is.same.group false/true for example same location vs. different
  calculate_qvals_and_effect_size<-function(sharing_df) {
    qval_df<-sharing_df %>% 
      group_by(Species) %>% 
      summarise(species_size=n(), 
               # effect_size= ,
                price_avg = ifelse(n_distinct(is.same.current_province)==2, 
                                                     wilcox.test(average_score[is.same.current_province=="FALSE"], average_score[is.same.current_province=="TRUE"], alternative = "less")$p.value, 
                                                     0-10)) %>%
      mutate(pval=ifelse(price_avg>0, price_avg,"NA")) %>%
      select(!price_avg)
    qval_df$qval<-p.adjust(qval_df$pval, method = "BH", n = sum(qval_df$pval!="NA"))
    return(qval_df)
  }
  
  #function to add stars for significance levels
  add_stars<-function (sharing_df, qval_df) { #Function to add stars correspoinding to qvals + merge with the sharing df's
    sharing_df<- left_join(sharing_df, qval_df, by="Species") %>%
      mutate(stars = ifelse(is.na(qval), " ", 
                            ifelse(qval<0.0000005, "****",
                              ifelse(qval<0.00005, "***",
                                 ifelse(qval<0.005, "**",
                                    ifelse(qval<0.05, "*","NS"))))))
    return(sharing_df)
  }
  
  # function to plot the synteny scores of within province vs. Between province, for each SRG with wilcoxon pval hows only significant SRGs (wilcoxon)
  ploting_location_sharing<-function(stars_df ) { 
    stars_df<-stars_df %>% filter(stars != " " , stars != "NS")
    stars_df$Species<- reorder(stars_df$real_Species,-stars_df$qval)
    ggplot(stars_df, aes(x=Species, y=average_score, fill=factor(is.same.current_province))) + 
      geom_boxplot(outlier.size=0) +
      ylab("Synteny Score") + xlab("Species") +
      geom_text(data=stars_df, aes(y=1.03, label=stars), col='black', size=6) +
      coord_flip() + 
      ggtitle("Synteny scores between and within provinces (excluding MIPs)") +
      guides(fill=guide_legend(title=NULL)) +
      #geom_hline(yintercept =  0.94, linetype = "dashed", color = "red", size =0.8) +
      scale_fill_manual(values=c("#999999", "#56B4E9"), labels=c("Different", "Same"))
  }
  

  
  #################################
  #    analyisis starts here      #
  #################################
  #################################
  #.  part 1: modify metadata already contained in the syntracker output
  #################################
  
  # add GroupE to the tables = infant/Mother classification for each sample
  GroupE=list()
  GroupE[as.character(metadata$Sample)]<-as.character(metadata$Type)
  cophylogeny_out_filtered<-mapply(add_groupE_mod, cophy_out, names(cophy_out), MoreArgs =  list(GroupE), SIMPLIFY = F)
  names(cophylogeny_out_filtered)<-names(cophy_out)
  
  
  # next, filter genomes with less than 100 regions comparisons, remove empty list elements
  for (i in length(cophylogeny_out_filtered):1) { 
    if (class(cophylogeny_out_filtered[[i]]) != "data.frame") {
      cophylogeny_out_filtered[[i]]<-NULL
    }
    if (nrow(cophylogeny_out_filtered[[i]]) < 100) {
      cophylogeny_out_filtered[[i]]<-NULL
    }
  }
  
  
  #######
  # PART 3: sub select n regions/pairwise, convert to a single DF
  #######
  
  #subsample n regions/pairwise comparions
  subsampled_cophylogeny<-mapply(subsample_regions_mod, cophylogeny_out_filtered)
  sc_copy<-subsampled_cophylogeny
 
   #remove SRGs which are empty OR have <25 pairwise comparisons
  for (i in length(subsampled_cophylogeny):1) { 
    if (subsampled_cophylogeny[[i]]=="remove this one") {
      #print(subsampled_cophylogeny[[i]])
      print(names(subsampled_cophylogeny)[i])
      subsampled_cophylogeny[[i]]<-NULL
    }
    else if (nrow(subsampled_cophylogeny[[i]])<10) {
      print(names(subsampled_cophylogeny)[i])
      subsampled_cophylogeny[[i]]<-NULL
    }
  }
  
  subsampled_cophylogeny_with_SRG_name<-mapply(add_names_to_col, subsampled_cophylogeny, names(subsampled_cophylogeny), SIMPLIFY = F) #add srg names
  subsdf<-bind_rows(subsampled_cophylogeny_with_SRG_name) #turn to one big dataframe
  subsdf<-subsdf %>% filter(is.same.GroupA == "FALSE") # remove mother-infant pairs
  
  
  #################
  #   PART 4: Geo locating: adding geolocating to the new metadata file,
  #################
  
  new_metadata<- new_metadata_Liam_currated %>% 
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Baden-Wuerttenburg" , "Baden-Wuerttemberg")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Daclac" , "dak lak")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Bacgiang" , "Bac Giang")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Laichau" , "Lai chau")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Namdinh" , "Nam dinh")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Nghean" , "Ng hean")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Ninhbinh" , "Ninh binh")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Oggooue-Ivindo" , "Ogooue Ivindo")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Oggooue-Lolo" , "Ogooue Lolo")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Phutho" , "Phu tho")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Quangninh" , "Quang ninh")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Thaibinh" , "Thai binh")) %>%
    mutate(CurrState.Province = replace(CurrState.Province, CurrState.Province == "Vinhphuc" , "Vinh phuc")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Caugiay" , "Cau giay")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Chuaboc" , "Chua boc")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Daton" , "Da ton")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Donganh" , "Dong anh")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Dongda" , "Dong da")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Ducthuan" , "Duc thuan")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Gialam" , "Gia lam")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Haibatrung" , "Hai ba trung")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Hoangcau" , "Hoang cau")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Hoangliet" , "Hoang liet")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Hoangmai" , "Hoang mai")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Kienhung" , "Kien hung")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Kimma" , "Kim ma")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Lethanhnghi" , "Le thanh nghi")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Loduc" , "Lo duc")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Mailam" , "Mai lam")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Namtuliem" , "Nam tu liem")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Ninhnhat" , "Ninh nhat")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Nioye2" , "Nioye")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Tathanhoai" , "Ta thanh oai")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Thanhtri" , "Thanh tri")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Thanhxuan" , "Thanh xuan")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Truongdinh" , "Truong dinh")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Vutongphan" , "Vu tong phan")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Walddorfhaeslach" , "Walddorfhäslach")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Yenso" , "Yen so")) %>%
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Bactuliem" , "Bac tu liem")) %>% 
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Weil_im_Schoenbuch" , "Weil_im_Schonbuch")) %>% 
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Trungliet" , "Trung liet")) %>% 
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Nurtigen" , "Nurtingen")) %>% 
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Nuertingen" , "Nurtingen")) %>% 
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Moessingen" , "Mossingen")) %>% 
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Korntal-Muenchingen" , "Korntal-Munchingen")) %>% 
    mutate(CurrVillage.City = replace(CurrVillage.City, CurrVillage.City == "Bactuliem" , "Bac Tu Liem")) %>% 
    unite( CurrLocation, CurrVillage.City, CurrCountry, sep = ", ", remove = FALSE, na.rm = FALSE)  
  
  #find geo coordinates for each towns
  
  locations <- as.character(levels(as.factor(new_metadata$CurrLocation))) #create a df with provinces
  geo_locations<-geocode_OSM(locations, as.data.frame = TRUE, keep.unfound = TRUE) 
  geo_locations <- geo_locations %>% select(query,lat,lon)
  
  #manually add missing locations
  manual_locations<- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(manual_locations)<-c("query","lat","lon")
  manual_locations[1,] <- c("Chiengcoi, Vietnam",21.318699508920968, 103.91343354292853)
  manual_locations[2,] <- c("Chiengden, Vietnam",21.394628514537096, 103.84041799611192)
  manual_locations[3,] <- c("Chiengsang, Vietnam",21.08110537566546, 104.23366144594792)
  manual_locations[4,] <- c("Chiengxom, Vietnam",21.37992358419129, 103.91941682421505)
  manual_locations[5,] <- c("Huala, Vietnam" ,21.27858142118783, 103.88756283178886)
  manual_locations[6,] <- c("Paga, Gabon",-2.0647368181093086, 10.616409334618066)
  manual_locations[7,] <- c("Soga, Gabon", -1.2321163375181936, 11.546449308755625)
  manual_locations[8,] <- c("Tonthattung, Vietnam", 21.37992358419129, 105.81071419574725)
  manual_locations[9,] <- c("Yennguu, Vietnam" ,20.951550644578518,105.83946358091217)
  manual_locations[10,] <- c("NA, Gabon"  ,NA,NA)
  manual_locations[11,] <- c("NA, Vietnam"  ,NA,NA)
  manual_locations[12,] <- c("NA, NA"  ,NA,NA)
  manual_locations[13,] <- c("Truong dinh, Vietnam",20.995296367750456, 105.84807012309781)
  
  #replace some missing location with the manual location
  geo_locations[geo_locations$query %in% manual_locations$query, ] <- manual_locations
 
  # add geo locations to the new_metadata table
  new_metadata_backup<-new_metadata
  new_metadata <- new_metadata %>% select(ParticipantID, SamplingLocation_1, CurrLocation,CurrCountry,CurrState.Province,CurrVillage.City, Age_years, Age_months, Expat)
  new_metadata<-left_join(new_metadata,geo_locations, by=c("CurrLocation"="query"))
  new_metadata$Expat<-replace(new_metadata$Expat, new_metadata$Expat=="Expat", 1) 
  ##########
  #   part 5: combine the info in the new_metadata file with the pairwise comparisons dataframe 
  ##########
  
  # create and populate named lists of participantID -> different metadata fields

  #first, declare the lists
  current_country<-list() # current country
  current_province<-list() #  state/province
  current_town<-list() # town/village
  Expat<-list() # expat yes/no
  current_lon<-list() #geo coordiantes
  current_lat<-list() 

  # Fill the lists 
  current_country[as.character(new_metadata$ParticipantID)]<-as.character(new_metadata$CurrCountry)
  current_province[as.character(new_metadata$ParticipantID)]<-as.character(new_metadata$CurrState.Province)
  current_town[as.character(new_metadata$ParticipantID)]<-as.character(new_metadata$CurrVillage.City)
  Expat[as.character(new_metadata$ParticipantID)]<-as.character(new_metadata$Expat)
  current_lon[as.character(new_metadata$ParticipantID)]<-as.character(new_metadata$lon)
  current_lat[as.character(new_metadata$ParticipantID)]<-as.character(new_metadata$lat)
  
  # assigne to subsdf (complete table holding all the filtered and relevant pairwise comparisons, with APSS values)
  # but first, add columns to the table
  tmp_df<-data.frame("sample1" = character(), 
                     "sample2"= character(),
                     "current.country.1" = character(), 
                     "current.country.2" = character(), 
                     "current.province.1" = character(), 
                     "current.province.2" = character(),
                     "current.town.1" = character(), 
                     "current.town.2"  = character(), 
                     "current.long.1"  = numeric(), 
                     "current.long.2" = numeric(), 
                     "current.lat.1"  = numeric(), 
                     "current.lat.2" = numeric(),
                     "Expat.1"=character(), 
                     "Expat.2" =character(), 
                     "is.same.current.country"=logical(), 
                     "is.same.current_province"=logical(), 
                     "is.same.current.town"=logical(),
                     "are.any.expats"=logical(), 
                     "Geo.Dist"=numeric())
  
  for (i in 1:nrow(subsdf)) {
    sample1<-as.character(subsdf$sample1[i]) 
    sample2<-as.character(subsdf$sample2[i])
    current.country.1<-current_country[sample1][[1]]
    current.country.2<-current_country[sample2][[1]]
    current.province.1<-current_province[sample1][[1]]
    current.province.2<-current_province[sample2][[1]]
    current.town.1<-current_town[sample1][[1]]
    current.town.2<-current_town[sample2][[1]]
    current.long.1<-current_lon[sample1][[1]]
    current.long.2<-current_lon[sample2][[1]]
    current.lat.1<-current_lat[sample1][[1]]
    current.lat.2<-current_lat[sample2][[1]]
    Expat.1<-Expat[sample1][[1]]
    Expat.2<-Expat[sample2][[1]]
    is.same.current.country<-current.country.1 == current.country.2
    is.same.current.province<-current.province.1 == current.province.2 
    is.same.current.town<-current.town.1 == current.town.2 
    are.any.expats<- ifelse(Expat.1==1 | Expat.2==1, "TRUE","FALSE")
    
    temprow<-data.frame(sample1,sample2,current.country.1,current.country.2,current.province.1,current.province.2,
                        current.town.1,current.town.2,current.long.1,current.long.2,
                        current.lat.1,current.lat.2, Expat.1,Expat.2,
                        is.same.current.country,is.same.current.province ,
                        is.same.current.town, are.any.expats)
    tmp_df<-rbind(tmp_df, temprow)
  }   
  
  full_subsdf<-cbind(subsdf, tmp_df) # this is the full table, with current locations, etc. 
  colnames(full_subsdf)[1]<-"sample1"
  colnames(full_subsdf)[2]<-"sample2"
  
  
  full_subsdf$current.lat.1<-as.numeric(full_subsdf$current.lat.1)   
  full_subsdf$current.lat.2<-as.numeric(full_subsdf$current.lat.2) 
  full_subsdf$current.long.1<-as.numeric(full_subsdf$current.long.1)
  full_subsdf$current.long.2<-as.numeric(full_subsdf$current.long.2)
  
  
  full_subsdf<-full_subsdf %>% rowwise() %>%
    mutate(geodist=approx_distances(c(current.long.1,current.lat.1), c(current.long.2,current.lat.2), 
                                    target = "km", projection = "EPSG:4326"))
  
  
  
  ####################################
  
  #########################
  # prepare a box plot of synteny scores: WITHIN-PROVINCE vs. BETWEEN-PROVINCES (IN THE SAME COUNTRY) 
  # Each species is a double box. 
  

  
  # calculate q-values for differences between groups within each species, add stars for significance. 
  
  qval_df_same_location<-mapply(calculate_qvals_by_location, list(subsdf), SIMPLIFY = F)
  stars_df_same_location<-mapply(add_stars, list(subsdf),qval_df_same_location, SIMPLIFY = F)
  stars_df_same_location[[1]]$real_Species<-unlist(oldname_to_species[ stars_df_same_location[[1]]$Species])
  
  #do the actual plotting:
  
  star_plots_by_location<-mapply(ploting_location_sharing, stars_df_same_location, SIMPLIFY = F)

  