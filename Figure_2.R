# script to comapre the APSS based trees in figure 2 to the MASH-based trees published by Abram et al. 

# requirements
library(tidyverse)
library(ggtree)
library(ape)
library(ggplot2)
library(TreeDist)
library(phangorn)
library(TreeTools)

#functions: 
#function to turn the final table, with average pairwise synteny scores,  to a symetric distance matrix. 
dist_mats<-function(dfss) { #input is a list of tables, output of SynTracker
  print("runing dist mats !")
  
  #select relevant columns:
  dfs<-dfss %>% ungroup() %>% 
    select(sample1,sample2,average_score)  
  dfs<-as.data.frame(dfs) #create a second df, with the same values but the samples are in reverse order (to create a symetric matrix)
  dfs2<-data.frame(dfs$sample2,dfs$sample1, dfs$average_score)
  colnames(dfs2)<-c("sample1","sample2", "average_score")
  
  newdfs<-rbind(dfs, dfs2) %>% #bind the two dataframes
    spread(key = sample2, value = average_score) #convert to wide format
  row.names(newdfs)<-newdfs$sample1
  newdfs$sample1<-NULL #remove column with sample names
  for (i in 1:nrow(newdfs)) { #add 1 for each comparison of a sample against itself
    for (j in 1:ncol(newdfs)) {
      if (i==j) {
        newdfs[i,j]<-1
      }
    }
  }
  return(newdfs)
}

#importing data
metadata<-read.csv(file = "10_per_group_metadata.tab", sep="\t", header = TRUE)
fulltree_mash<-read.tree(file="Supplementary_Data_Fig_3c_10k_newick.nwk") #import the MASH based full tree from Abram et al 2021
md5_to_gca<-read.table(file = "Supplementary_Data_Fig_3c_md5_gca_list.csv", header = T, sep = ",") #same, imported from Abram et al
output_object<-readRDS(file = "e.coli_single_genomes_output") #this is the output of the synteny analysis of the 10 genomes per phylogroup.

#convert md5-gca to list, change tip labels from MD5 to gca
dictionary<-list()
dictionary[md5_to_gca$MD5]<-md5_to_gca$Assembly.Accession
for (i in 1:length(dictionary)) { #change tip labels
  fulltree_mash$tip.label[fulltree_mash$tip.label==md5_to_gca[i,1]]<-dictionary[md5_to_gca[i,1]][[1]]
}

#proning ANI tree to keep relevnt nodes. 
merged_identifiers<-metadata[,-(2)]
newtree<-keep.tip(fulltree_mash, merged_identifiers$Sample) # subset the full mash based tree to keep only the 140 genomes in the synteny analysis. 

#calculate synteny trees, plot against ani tree

mat20<-dist_mats(output_object$blastcmddb_output$`pairwise comparison tables`$`blastcmddb_output 20 regions/pairwise`)
mat30<-dist_mats(output_object$blastcmddb_output$`pairwise comparison tables`$`blastcmddb_output 30 regions/pairwise`)
mat40<-dist_mats(output_object$blastcmddb_output$`pairwise comparison tables`$`blastcmddb_output 40 regions/pairwise`)
mat60<-dist_mats(output_object$blastcmddb_output$`pairwise comparison tables`$`blastcmddb_output 60 regions/pairwise`)
mat80<-dist_mats(output_object$blastcmddb_output$`pairwise comparison tables`$`blastcmddb_output 80 regions/pairwise`)
mat100<-dist_mats(output_object$blastcmddb_output$`pairwise comparison tables`$`blastcmddb_output 100 regions/pairwise`)
mat200<-dist_mats(output_object$blastcmddb_output$`pairwise comparison tables`$`blastcmddb_output 200 regions/pairwise`)

mat.list<-list(mat20,mat30,mat40,mat60,mat80,mat100,mat200) #change the distances of the tree, so 0 is the closest
names(mat.list)<-c('20','30','40','60','80','100','200')
matfunc<-function(mat) {abs(mat-1)}
synteny_matrices<-mapply(matfunc, mat.list, SIMPLIFY = F)


tree.func<-function(synteny_matrix){upgma(as.matrix(synteny_matrix))}#make the trees
upgma_trees<-mapply(tree.func,synteny_matrices, SIMPLIFY = F)


# calculate between-tree distances + pvalues .
calc_dist<-function(real_trees, newtree) {
  real_distance<-ClusteringInfoDistance(real_trees, newtree, normalize = TRUE) #note to self: check if I need to convert the format of newtree
  nRep <- 10000 # Use more replicates for more accurate estimate of expected value
  randomTrees <- lapply(logical(nRep), function (x) RandomTree(real_trees$tip.label))
  randomDists <- ClusteringInfoDistance(newtree, randomTrees, normalize = TRUE)
  nThisSimilar <- sum(randomDists < real_distance)
  pValue <- nThisSimilar / nRep
  return(c(pValue,real_distance))
  #return(real_distance)
}

my_pvals_and_dists<-mapply(calc_dist, upgma_trees, list(newtree))



# plot the trees. 
t1<-upgma_trees$`30`
t2<-newtree

merged_identifiers<-metadata[,-(2)]
colnames(merged_identifiers)<-c("label", "Group")

T1 <- ggtree(t1 ) +    
  theme_tree2(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) +   
  geom_tiplab()   


T2 <- ggtree(t2) +   
  theme_tree2(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) +   
  geom_tiplab(hjust =1) +   
  scale_x_reverse()     


d1 = T1$data[T1$data$isTip,]  
d1$x[] = 1  
d2 = T2$data[T2$data$isTip,]  
d2$x[] = 2  

TTcon <- rbind(d1, d2)  
TTconcon<-left_join(TTcon,merged_identifiers, by="label")
T1 = ggtree(t1, layout='roundrect') %<+%   
  merged_identifiers +
  theme_tree2(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) + 
  geom_tippoint(aes(label = Group, col=Group)) +
  xlab("synteny distance")# +
#geom_tiplab(align = TRUE) # if required +   xlim(0,5.2) 
T1mod = ggtree(t1) %<+%   
  merged_identifiers +
  theme_tree2(legend.position='left', plot.margin = unit(c(0,0,0,0),"cm")) + 
  geom_tippoint(aes(label = Group, col=Group)) 


T2 = ggtree(t2, layout='roundrect') %<+%   
  merged_identifiers +
  theme_tree2(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) +  
  geom_tippoint(aes(label = Group, col=Group)) +
  #geom_tiplab(hjust =1, align = TRUE) +   
  xlab("MASH-based distance") +
  scale_x_reverse() # if required add within scale_x_reverse() -> limits = c(5.2, 0))  



L1 = ggplot(TTconcon, aes(x = x, y = y, colour=Group, group = label)) + geom_line() +   
  theme_void() + theme(legend.position="none", plot.margin = unit(c(1,0,1,0),"cm"))  

cowplot::plot_grid(T1, L1 ,T2, nrow = 1, align = "hv")

# make scatter plot
forplot<-as.data.frame(t(my_pvals_and_dists))
forplot<-forplot %>% rownames_to_column()
colnames(forplot)<-c("Regions","pval","RFD")
corel<-cor(log(as.numeric(forplot$Regions)), forplot$RFD, method = "spearman")

p<-ggplot(forplot, aes(x = log(as.numeric(Regions)), y = RFD)) +
  theme_bw()+
  geom_point(color = "black",
             size = 3,
             show.legend = FALSE) +
  stat_smooth(method = "lm",
              col = "red",
              se = FALSE,
              size = 1.5,
              linetype = "dashed") +
  labs(x = "log(Subsampled Regions)", y="Robinson Foulds Distance", size = 6) +
  geom_text(x=4.5, y=0.34, label="ρ = -0.93", size = 6, fontface="italic")

grDevices::cairo_pdf("RFD.Regions.pdf")
p
dev.off()
