library(tidyverse)
library(DESeq2)
library(viridis)

#########################################
#load in raw counts from each dataset

#read in gene counts from Li et al. 2015
countdataLi <- read_csv("Li2015GeneCounts.csv")
colnames(countdataLi) <- c("GeneID","Li.S6.1","Li.S7.1","Li.S1.1","Li.S2.1",
                           "Li.S3.1","Li.S1.2","Li.S2.2","Li.S3.2","Li.S4.2",
                           "Li.S5.2","Li.S6.2","Li.S7.2","Li.S4.1","Li.S5.1")

#read in gene counts from Wang et al. 2014
countdataWang <- read_csv("Wang2014GeneCounts.csv")
names <- str_sub(colnames(countdataWang)[-1],31,-35)
colnames(countdataWang) <- c("GeneID", paste0("Wang.S", str_sub(names,6,-1), ".", str_sub(names,1,1)))

#combine into one dataset
countdataBoth <- as.data.frame(merge(countdataLi, countdataWang, by="GeneID")) %>%
  mutate(GeneID = str_sub(GeneID, 1, -9)) %>% 
  column_to_rownames("GeneID") %>% 
  as.matrix()

group <- factor(str_sub(colnames(countdataBoth), 1, -3))
coldata <- data.frame(row.names=colnames(countdataBoth), group) 

#make DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = countdataBoth, colData = coldata, design = ~0+group)
#run the pipeline
dds <- DESeq(dds)
#filter low counts; >= 100 eliminates ~half of genes
keep <- rowSums(counts(dds)) >= 100 
sum(keep)
dds <- dds[keep,]
resultsNames(dds)


#############################
#Fig 1A; PCA
#plots PCA of all 11 Wang et al. 2014 segments (4 reps each)

rld <- vst(dds, blind=FALSE)

##LI DATA
#pca <- plotPCA(rld[,1:14], intgroup="group", returnData = TRUE) %>% 
#mutate(Segment = as.numeric(str_sub(group.1, 5, -1)))

#Wang data
pca <- plotPCA(rld[,15:58], intgroup="group", returnData = TRUE) %>% 
  mutate(Segment = as.numeric(str_sub(group.1, 7, -1)))
percentVar <- round(100 * attr(pca, "percentVar"))

ggplot(pca, aes(PC1,PC2, color = Segment)) +
  geom_point(size = 4) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_gradientn(colors = c("#E1E858","#129B01","#204B13"), name = "", breaks=c(1:11)) +
  theme_minimal() +
  theme(legend.key.height= unit(1.3, 'cm'))

#############################
#Fig 1B; DEGs
#calculate pairwise DEGs between segments

#DEGs of all pairwise comparisons
#loops over each pair (55 in total), combines into one dataframe
#p<0.05, adjusted with BH; lfcShrink run after
allComp <- combn(1:11,2)
for(i in 1:ncol(allComp)){
  if(i==1){DEGs<-data.frame()}
  con <- c("group", paste0("Wang.S",allComp[1,i]), paste0("Wang.S",allComp[2,i]))
  res <- results(dds, contrast=con, independentFiltering=TRUE, 
                 alpha=0.05, pAdjustMethod="BH")
  res <- lfcShrink(dds, contrast=con, res=res, type="ashr")
  resTemp <- data.frame(res) %>% 
    rownames_to_column("GeneID")
  
  if(nrow(resTemp) != 0){
    resTemp$Comparison <- paste(con[2], con[3], sep = "_")
    DEGs <- rbind(DEGs, resTemp)
  }
}

#filter with stricter thresholds
DEGsFilt <- DEGs %>% 
  filter(padj < 0.001 & abs(log2FoldChange) > 1)
length(unique(DEGsFilt$GeneID))
#20415 genes before filtering
#13066 DEGs after filtering (padj < 0.001, lfc>1)

dgeList <- unique(DEGsFilt$GeneID)


#Fig 1B; PCC
#calculate Pearson Correlation Coefficient between each pair of segments

#function to average replicates
avg <- function(x){t(apply(x, 1, function(y) tapply(y, str_sub(colnames(x),1,-3), mean)))}

#get normalized Wang counts and average
hm <- counts(dds[,15:58],normalized=T) %>% 
  avg()

#get correlation matrix
hmCor <- cor(hm) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample1") %>% 
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "Cor") %>% 
  mutate(Sample1 = factor(Sample1, paste0("Wang.S",1:11)),
         Sample2 = factor(Sample2, paste0("Wang.S",1:11))) %>% 
  arrange(Sample1, Sample2)

#plot just PCC data
ggplot(hmCor, aes(Sample1, Sample2, fill = Cor, color = Cor)) +
  geom_tile() +
  scale_fill_viridis(option = "D",name = "", breaks = c(0.6,0.8,1)) +
  scale_color_viridis(option = "D",name = "", breaks = c(0.6,0.8,1)) +
   labs(y="Base                         Tip",x="Base                         Tip") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom")

#Fig 1B; PCC + DEGs
#make 1 heatmap with PCC on top and DEGs on bottom

#get filtered DEG counts and reshape to fit with PCC data
#add in rows with 0 DEGs
temp <- data.frame(Sample1 = c("Wang.S9", "Wang.S10"),
                   Sample2 = c("Wang.S10", "Wang.S11"),
                   n = 0)
degCor <- DEGsFilt %>% 
  dplyr::select(Comparison) %>% 
  group_by(Comparison) %>% 
  mutate(n = n()) %>% 
  distinct() %>% 
  separate(col = Comparison, into = c("Sample1", "Sample2"), sep = "_") %>% 
  mutate(Sample1 = factor(Sample1, levels = paste0("Wang.S",1:11)),
         Sample2 = factor(Sample2, levels = paste0("Wang.S",1:11))) %>% 
  rbind(temp) %>% 
  dplyr::select(Sample1,Sample2, n) %>% 
  mutate(n = -n) #swap sign so the heatmap is interpretable

#rescale to same range (for plotting)
library(scales)
degCor$n <- rescale(degCor$n, to = c(min(hmCor$Cor), max(hmCor$Cor)))

#merge and plot; all will be on same scale
degCor %>% 
  dplyr::rename("Cor" = "n") %>% 
  rbind(hmCor) %>% 
  distinct(Sample1, Sample2, .keep_all = TRUE) %>% 
ggplot(aes(Sample1, Sample2, fill = Cor, color = Cor)) +
  geom_tile() +
  scale_fill_viridis(option = "D",name = "", breaks = c(0.4,0.6,0.8,1), limits = c(0.4,1)) +
  scale_color_viridis(option = "D",name = "", breaks = c(0.4,0.6,0.8,1), limits = c(0.4,1)) +
   labs(y="Base                         Tip",x="Base                         Tip") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")


#get DEG legend for figure 
#ignore the heatmap, just use the legend
DEGsFilt %>% 
  dplyr::select(Comparison) %>% 
  group_by(Comparison) %>% 
  mutate(n = n()) %>% 
  distinct() %>% 
  separate(col = Comparison, into = c("Sample1", "Sample2"), sep = "_") %>% 
  mutate(Sample1 = factor(Sample1, levels = paste0("Wang.S",1:11)),
         Sample2 = factor(Sample2, levels = paste0("Wang.S",1:11))) %>% 
  rbind(temp) %>% 
  dplyr::select(Sample1,Sample2, n) %>% 
ggplot(aes(Sample1, Sample2, fill = n, color = n)) +
  geom_tile() +
  scale_fill_viridis(option = "D",name = "", direction = -1) +
  scale_color_viridis(option = "D",name = "", direction = -1) +
   labs(y="Base                         Tip",x="Base                         Tip") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right")


#########################################
#Fig 1C+D
#heatmap of all DEGs along with clustered profiles
#this section will hierarchically cluster the DEGs, define expression patterns using dynamicTreeCut,
#and plot line graphs of expression patterns

#start after DEGs have been calculated
#average reps and scale counts (zscore)
scaleAvg <- function(x){t(scale(apply(x, 1, function(y) tapply(y, str_sub(colnames(x),1,-3), mean))))}
hm <- counts(dds, normalized = T) %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% DEGsFilt$GeneID) %>% 
  dplyr::select(contains("Wang")) 
hm <- scaleAvg(hm)

#hierarchically cluster
dist <- dist(hm)
hmc <- hclust(dist)
#reorder heatmap based on Optimal Leaf Ordering (does not change tree, just the order of subtrees)
#not necessary for expression patterns, just to make the heatmap look more cohesive
library(seriation)
hmcOlo <- reorder(hmc, dist)

#split ordered DEGs into clusters (expression patterns)
library(dynamicTreeCut)
cut <- cutreeDynamic(dendro = hmcOlo, minClusterSize = 800,  
                     method = "hybrid", distM = as.matrix(dist), deepSplit = 0,
                     pamStage = T, pamRespectsDendro = T)

#get clusters, remove unassigned genes
Clusters <- data.frame("Gene" = rownames(hm), "Cluster" = cut)[hmcOlo$order,] %>% 
  filter(Cluster != 0)

#this block reorders the clusters (default ordering is by size; we want similarity),
#reshapes the dataframe for plotting, and calculates the mean of each cluster
order <- Clusters %>% 
  mutate(Cluster = factor(Cluster, levels = unique(Clusters$Cluster))) %>% 
  dplyr::count(Cluster) %>% 
  mutate(NewCluster = 1:(length(unique(Clusters$Cluster))))
clusterDEG <- data.frame(hm) %>% 
  rownames_to_column("Gene") %>% 
  merge(Clusters, by = "Gene") %>% 
  pivot_longer(!c("Gene","Cluster"), names_to = "Sample", values_to = "CPM") %>% 
  mutate(Sample = factor(Sample)) %>% 
  merge(order, by = "Cluster") %>% 
  dplyr::select(-Cluster) %>% 
  mutate(ClusterLabel = paste0("Cluster ", NewCluster, " (n=", n, ")"))
clusterMeans <- clusterDEG %>% 
  group_by(ClusterLabel,NewCluster,Sample) %>% 
  dplyr::summarise(mean = mean(CPM))
clusterNames <- unique(clusterDEG$ClusterLabel)
names(clusterNames) <- unique(clusterDEG$NewCluster)
lev = paste0("Wang.S",1:11)
ggplot() +
  geom_line(data = clusterDEG, aes(factor(Sample, levels = lev), CPM, group = Gene), alpha = 0.1) +
  geom_line(data = clusterMeans, aes(factor(Sample, levels = lev), mean, group = ClusterLabel), color = "red") +
  labs(title = paste0("Wang et al. DEGs (fdr<0.001,abs(log2FC)>1), n=",sum(order$n)), 
       y = "z-score",
       x= "") +
  facet_wrap(~NewCluster, labeller = as_labeller(clusterNames)) +
  theme(axis.text.x = element_text(angle = 90))

#make new list of DEGs assigned to clusters (those with low similarity to any cluster were excluded)
dgeList2 <- unique(clusterDEG$Gene)

#some of those clusters can be regrouped; easier to visualize and work with
clusterDEG2 <- clusterDEG %>% 
  mutate(Cluster = case_when(NewCluster %in% c(2,3) ~ 1,
                             NewCluster %in% c(1,4) ~ 2,
                             NewCluster %in% c(6) ~ 3,
                             NewCluster %in% c(5,7) ~ 4,
                             NewCluster %in% c(8,9,10) ~ 5))
clusterMeans2 <- clusterDEG2 %>% 
  group_by(Cluster,Sample) %>% 
  dplyr::summarise(mean = mean(CPM))
#give proper names
clusterNamesNew <- c("Blade Tip\n(n=2,036)",
                     "Mid-Blade\n(n=2,418)",
                     "Blade Base\n(n=1,573)",
                     "Maturation Zone\n(n=2,882)",
                     "Division and Elongation\nZone (n=3,565)")
names(clusterNamesNew) <- 1:5

#replot vertically stacked with better color scheme; 
#scales::squish allows us to color each line segment with the mean zscore
library(scales)
ggplot() +
  geom_line(data = clusterDEG2, aes(factor(Sample, levels = lev), CPM, group = Gene), alpha=0.05) +
  geom_line(data = clusterMeans2, aes(factor(Sample, levels = lev), mean, group = Cluster, color = mean), size = 1.2) +
  scale_color_viridis_c(limits = c(-1.5, 1), oob = scales::squish) +
  labs(y = "z-score",
       x= "Base                   Tip") +
  facet_wrap(~Cluster, labeller = as_labeller(clusterNamesNew), nrow = 5) + 
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.position = "none",
        strip.text = element_text(size=10),
        panel.grid.major.x = element_blank())


####################################
#Fig 1E
#enriched GO Terms of each cluster

library(topGO)
#read in GO term data (received from Phytozome)
GO2geneID <- readMappings("MSURiceGoTerms.txt")
geneNames <- names(GO2geneID)

#run loop to get enriched terms for genes in each cluster; save in one dataframe
for(i in 1:length(unique(clusterDEG2$Cluster))){
  if(i==1){Res <- data.frame()}
  
  #Read in list of genes from clusters
  genes <- filter(clusterDEG2,Cluster == i)
  geneList <- factor(as.integer(geneNames %in% unique(genes$Gene)))
  names(geneList) <- geneNames

  #run GO analysis, biological processes only
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, 
                annot = annFUN.gene2GO, gene2GO = GO2geneID, nodeSize = 5)
  
  #Run fisher test with default weighted algorithm and filter by pval
  resultFis <- runTest(GOdata, statistic = "fisher")
  sig <- data.frame(resultFis@score) %>% 
    filter(resultFis.score < 0.05)
  
  #Generate results table
  temp <- GenTable(GOdata, weightFisher = resultFis, orderBy = resultFis, topNodes = nrow(sig))
  temp$Cluster <- i
  Res <- rbind(Res, temp)
}
#make final table
Res$weightFisher <- gsub("< ", "", Res$weightFisher)
allRes <- data.frame("Term" = Res$Term,
                     "pval" = as.numeric(Res$weightFisher),
                     "Cluster" = Res$Cluster)
allRes2 <- allRes %>% 
  filter(pval < 0.001)



