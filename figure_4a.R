#This is a script to generate scatter plot of popANI vs APSS, with density plots of each distance.

# also create a venn diagram of instrain/syntracker intersections
# input is the pairwise popANI and APSS tables for the huBif data set (processed inStrain and SynTracker output)

library(venneuler)
library(gridExtra)

pairwise_comp_ani<-new_instrain_out %>% select(genome, name1, name2, popANI, gtdbtk_phylum,gtdbtk_class,gtdbtk_family, gtdbtk_genus)
pairwise_comp_apss<-stars_df_same_family[[1]] %>% ungroup() %>% select(Species, sample1, sample2, average_score)
colnames(pairwise_comp_ani)<-c("Species","sample1","sample2", "popANI", "Phylum", "Class","Family","Genus")


pairwise_comp_ani<-pairwise_comp_ani %>% mutate(Species_2=str_c(
  str_split(Species,"_", simplify = T)[, 1], 
  "__",
  str_split(Species,"_", simplify = T)[, 3], 
  "__",
  str_split(Species,"_", simplify = T)[, 5])) 


# find the pairwise comparisons that are detected by both instrain and syntracker
# as the order of sample1-sample2 could differ in the two tables, I inner-join the tables 
# using the two options, then i combine the DFs. 

overlapping_comparisons_1<-inner_join(pairwise_comp_ani,pairwise_comp_apss, by=c("Species_2"="Species","sample1"="sample1", "sample2"="sample2"))

overlapping_comparisons_2<-inner_join(pairwise_comp_ani,pairwise_comp_apss, by=c("Species_2"="Species","sample1"="sample2", "sample2"="sample1")) 

overlapping_comparisons<-rbind(overlapping_comparisons_1,overlapping_comparisons_2) 

top_5_apss<-quantile(overlapping_comparisons$average_score, 0.95)
top_5_ani<-quantile(overlapping_comparisons$popANI, 0.95)

overlapping_comparisons <- overlapping_comparisons %>%
  mutate(dotcol= if_else(average_score>top_5_apss & popANI>top_5_ani,"purple", 
                         if_else(average_score>top_5_apss & popANI<top_5_ani,"red", 
                                 if_else(average_score<top_5_apss & popANI>top_5_ani,"blue","dark grey")
                         )
  ))



# start with figure
empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

#top plot: APSS density+line at median APSS
median_APSS <- overlapping_comparisons %>% 
  pull(average_score) %>% 
  median() %>%
  signif(4)

plot_top <- ggplot(overlapping_comparisons, aes(average_score, fill="red")) + 
  geom_density(alpha=.5) + 
  labs(x = NULL, y = "Density") +
  geom_vline(xintercept=median_APSS, color="red") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(legend.position = "none") 

#side plot: popANI density+line at median popANI
#first, find median ani
median_ani <- overlapping_comparisons %>% 
  pull(popANI) %>% 
  median() %>%
  signif(4)
#now do the plot
plot_side<-ggplot(overlapping_comparisons, aes(popANI)) + 
  geom_density(alpha=.3,fill="blue") + 
  labs(x = NULL, y = "Density") +
  geom_vline(xintercept=median_ani, color="blue") +
  coord_flip() +
  theme_bw() + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position = "none")
# calculate correlation
cor_coef<-cor(overlapping_comparisons$average_score, overlapping_comparisons$popANI, method = "kendall") %>% signif(4)

#Main plot, x-y scatter
scatter<-ggplot(data=overlapping_comparisons, aes(x=average_score, y=popANI, colour = dotcol)) + 
  geom_point() +
  xlim(0.5,1)+ylim(0.965,1) +
  scale_colour_identity() +
  labs(x = "APSS", y = "popANI") +
  annotate(geom="text", x=0.56, y=0.9975, label=paste(cor_coef, sep = ": "), color="red") +
  theme_bw() 

grid.arrange(plot_top, empty, scatter, plot_side, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))


#create venn diagram
vd <- venneuler(c(Instrain=nrow(pairwise_comp_ani), SynTracker=nrow(pairwise_comp_apss), "Instrain&SynTracker"=nrow(overlapping_comparisons))
                )
plot(vd)


