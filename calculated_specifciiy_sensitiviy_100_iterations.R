# This is a script to determine the APSS threshold between same and different strains, and to check the sensitivity and specificity of the tool 
# The input: 
# use the 100 dfs per species/sampling depth combination to run 100 ROC curves per combination
# store the best TH in a vector: average the TH per combination
# use the average TH as was done previously (old script)

################
#   requirements
################

library(pROC)

################
#.  functions
################

rocplots1<-function(dframes, names) { #function to rocplot and calculate optimal scores, return the optimal tresholds
  youden_vec<-c() #vector holding 100 youden best thresholds
  for (i in 1:100) {
    rocinf<-roc(dframes[[i]]$is.same.GroupA, dframes[[i]]$average_score, levels = c("FALSE", "TRUE"))
    best_th<-coords(rocinf, x="best", transpose = TRUE, best.method = "youden")   #find the best threshold 
    youden_vec<-c(youden_vec, best_th[1])
  }
  average_best_th<-mean(youden_vec)
  return(average_best_th)
} 



rocplots2<-function(dframes, names, tresholds) { #user provides threshold values. Returns the specificity/sensitivity.  
  number_regions<-as.character(dframes[[1]]$compared_regions[1])
  Best_threshold<-as.numeric(tresholds[number_regions])
  spec.sens.df<-data.frame (specificity=integer(), #DF holding 100  best specificity/sensitiviy values
                           sensitivity=integer()) 
  for (i in 1:100) { #repeat for each of the 100 dataframes per species/sampling depth combination
    print(i)
    rocinf<-roc(dframes[[i]]$is.same.GroupA, dframes[[i]]$average_score, levels = c("FALSE", "TRUE"))
    best_th<-coords(rocinf, x=Best_threshold, input="threshold", transpose = TRUE)   #find the best threshold 
    spec.sens.df[i,]<-best_th[2:3]
  }
  ave_values<-data.frame(specificity=mean(spec.sens.df$specificity), sensitivity=mean(spec.sens.df$sensitivity), Regions=number_regions)
  return(ave_values)
}  

rocplots3<-function(dframes, Best_threshold) { #user provides threshold values. Returns the specificity/sensitivity.  
  rocinf<-roc(dframes[[i]]$is.same.GroupA, dframes[[i]]$average_score, levels = c("FALSE", "TRUE"))
  best_th<-coords(rocinf, x=Best_threshold, input="threshold", transpose = TRUE)   #find the best threshold 
} 



####################
# Script starts here
####################

#########################################
# 1.  finding best optimal synteny score threshold using the training set: 
#########################################

# REMOVE tables without enough same-person and different person comparisons (at least 3 of each type for training, 2 for testing)

for (i in 1:length(reduced_trainlist1)) { #goes over each "species"
  for (j in length(reduced_trainlist1[[i]]):1) { # goes over each "sampling depth"
    # check only the first of 100 dataframes for each combination, 
    # as all 100 have the same number of cases
    testval<-as.data.frame(reduced_trainlist1[[i]][[j]][[1]]) %>% filter(is.same.GroupA == TRUE) %>% nrow() # checks the first of 100 dataframes in each combination, as all are the same, except average synteny score
    testval2<-as.data.frame(reduced_trainlist1[[i]][[j]][[1]]) %>% filter(is.same.GroupA == FALSE) %>% nrow()
    if (testval<3 | testval2<3 ) {
      reduced_trainlist1[[i]][[j]]<-NULL
    }
  }
}


out1<-list()
counter<-0
for (i in 1:length(reduced_trainlist1)) {
  out1[[i]]<-mapply(rocplots1, reduced_trainlist1[[i]] ,names(reduced_trainlist1[[i]]), SIMPLIFY = F)
  counter<-counter+1
}
names(out1)<-names(reduced_trainlist1)

#remove emtpy list elements 
for (element in 1:length(out1)) {
  ifelse(length(out1[[element]])==0, out1[element]<-NULL, "this is OK")
}
#calculate average threshold per given sampling depth. 

all_th<-as.data.frame(t(bind_cols(out1)))
colnames(all_th)<-c("threshold")
all_th<-all_th %>% rownames_to_column() %>% 
  rowwise %>% 
  mutate("subsample"=str_split(rowname, "\\...|\ ")[[1]][2]) %>%
  select(-rowname) %>%
  group_by(subsample) %>%
  summarise("Mean_th" = mean(threshold))

#convert the all_th dataframe to list: subsample_value=>synteny threshold
thresholds_list<-list() 
thresholds_list[as.character(all_th$subsample)]<-all_th$Mean_th



#########################################
#.  Use the thresholds from the training set on the test set. 
# Report sensitivity/specificity for the test set: 
#########################################



for (i in 1:length(reduced_testlist1)) { #goes over each "species"
  for (j in length(reduced_testlist1[[i]]):1) { # goes over each "sampling depth"
    # check only the first of 100 dataframes for each combination, 
    # as all 100 have the same number of cases
    testval<-as.data.frame(reduced_testlist1[[i]][[j]][[1]]) %>% filter(is.same.GroupA == TRUE) %>% nrow() # checks the first of 100 dataframes in each combination, as all are the same, except average synteny score
    testval2<-as.data.frame(reduced_testlist1[[i]][[j]][[1]]) %>% filter(is.same.GroupA == FALSE) %>% nrow()
    if (testval<3 | testval2<3 ) {
      reduced_testlist1[[i]][[j]]<-NULL
    }
  }
}


outtest1_temp<-list()
outtest1<-list()
outtest2<-list()

#1. make a list with regions subsampled => best threshold values (based on average best thresholds from the training set)
#2. Make a second list, based on specific thresholds (for each species/sample combination, use the corresponding TH)

for (i in 1:length(reduced_testlist1)) { 
  outtest1_temp[[i]]<-mapply(rocplots2, reduced_testlist1[[i]] ,names(reduced_testlist1[[i]]), MoreArgs =list(thresholds_list) , SIMPLIFY = F)
  outtest1[[i]]<-as.data.frame(bind_rows(outtest1[[i]]))
}
names(outtest1)<-names(reduced_testlist1)

for (i in 1:length(outtest1)) {
  outtest1[[i]]<-as.data.frame(bind_rows(outtest1[[i]]))
}

spec.sens.averages<-as.data.frame(bind_rows(outtest1))

all_th_test<-spec.sens.averages %>% 
  rowwise %>% 
  group_by(Regions) %>%
  mutate("Mean_Sensitivity" = mean(sensitivity), "Mean_Specificity" = mean(specificity),  
         "std_Sensitivity" = sd(sensitivity), "std_Specificity" = sd(specificity),  
         "ymin" = Mean_Specificity-std_Specificity,"ymax" = Mean_Specificity+std_Specificity, 
         "xmin" = Mean_Sensitivity-std_Sensitivity,"xmax" = Mean_Sensitivity+std_Sensitivity, 
         "xmaxcor"=ifelse(xmax>1, 1, xmax), "ymaxcor"=ifelse(ymax>1, 1, ymax))

# plot 
pdf(file = "specificity_sensitivity_based_on_average_th_with_EB_based_on_100_itereations.pdf", width=6.5, height = 5)
all_th_test %>% 
  mutate(Regions = factor(Regions, levels=c("15","20","40","60","80","120", "160", "200"))) %>%
  ggplot(aes(x=Mean_Sensitivity, y=Mean_Specificity,col=Regions))+
  geom_point(aes(size=1)) + scale_color_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin = ymin,ymax = ymaxcor),alpha=0.3,width=0.005) + 
  geom_errorbarh(aes(xmin = xmin,xmax = xmaxcor),alpha=0.3, height=0.005) 
dev.off()

pdf(file = "specificity_sensitivity_based_on_average_th_based_on_100_itereations.pdf", width=6.5, height = 5)
all_th_test %>% 
  mutate(Regions = factor(Regions, levels=c("15","30","40","60","80","120", "160", "200"))) %>%
  ggplot(aes(x=Mean_Sensitivity, y=Mean_Specificity,col=Regions, palette="Dark2"))+
  geom_point() + scale_color_brewer(palette="Dark2") 
dev.off()

# Do the same, this time with paired testing/training values per species and sample depth
# First, Create a list of thresholds per each combination

for (i in 1:length(out1)) { 
  print(i)
  outtest2[[i]]<-list()
  for (j in 1:length(out1[[i]])) {
    print(j)
    print(names(out1)[i])
    print(names(out1[[i]])[j])
    print(names(reduced_testlist1[names(out1)])[i])
    print(names(reduced_testlist1[[names(out1)[i]]][names(out1[[i]][j])]))
    ifelse(is.null(reduced_testlist1[[names(out1)[i]]][[names(out1[[i]][j])]]), 
           print("skip"),
           outtest2[[i]][[j]]<-mapply(rocplots3, reduced_testlist1[[names(out1)[i]]][names(out1[[i]][j])] , MoreArgs =list(out1[[i]][[j]][[1]]) , SIMPLIFY = F)
    )
  }  
}

# convert the list to a table
paired_thersholds<-as.data.frame(t(bind_cols(outtest2)))
colnames(paired_thersholds)<-c("threshold", "specificity", "sensitivity")
paired_thersholds<-paired_thersholds %>% rownames_to_column() %>% 
  rowwise %>% 
  mutate("subsample"=str_split(rowname, "\\...|\ ")[[1]][2]) %>%
  select(-rowname) %>%
  group_by(subsample) %>%
  mutate("Mean_Sensitivity" = mean(sensitivity), "Mean_Specificity" = mean(specificity), "Mean_Threshold" = mean(threshold), 
         "std_Sensitivity" = sd(sensitivity), "std_Specificity" = sd(specificity), "std_Threshold" = sd(threshold), 
         "ymin" = Mean_Specificity-std_Specificity,"ymax" = Mean_Specificity+std_Specificity, 
         "xmin" = Mean_Sensitivity-std_Sensitivity,"xmax" = Mean_Sensitivity+std_Sensitivity, 
         "xmaxcor"=ifelse(xmax>1, 1, xmax), "ymaxcor"=ifelse(ymax>1, 1, ymax), "Mean_th" = mean(threshold))
#plot
pdf(file = "specificity_sensitivity_based_on_paired_th_with_EB.pdf", width=6.5, height = 5)
paired_thersholds %>% 
  mutate(Regions_subsampled = factor(subsample, levels=c("15","20","30","40","60","80","100","120","140", "160", "180","200"))) %>%
  ggplot(aes(x=Mean_Sensitivity, y=Mean_Specificity,col=Regions_subsampled, palette="Dark2"))+
  geom_point(aes(size=1)) + scale_color_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin = ymin,ymax = ymaxcor),alpha=0.3,width=0.005) + 
  geom_errorbarh(aes(xmin = xmin,xmax = xmaxcor),alpha=0.3, height=0.005) 
dev.off()

pdf(file = "specificity_sensitivity_based_on_paired_th.pdf", width=6.5, height = 5)
paired_thersholds %>% 
  mutate(Regions_subsampled = factor(subsample, levels=c("15","20","30","40","60","80","100", "120","140", "160", "180","200"))) %>%
  ggplot(aes(x=Mean_Sensitivity, y=Mean_Specificity,col=Regions_subsampled, palette="Dark2"))+
  geom_point() + scale_color_brewer(palette="Dark2") 
dev.off()
