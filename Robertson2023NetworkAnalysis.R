library(tidyverse)
library(DESeq2)

###########################
#remove DEGs only DE in segment 1; skews the whole network
DEGsFiltNoS1 <- DEGsFilt %>% 
  filter(str_sub(Comparison,1,8) != "Wang.S1_")
dgeListNetwork <- dgeList2[dgeList2 %in% DEGsFiltNoS1$GeneID]

###########################
#Network setup

#MSU and RAP ID conversion
#CisBP rice TF information is in RAP format
library(riceidconverter)
riceIDConv <- RiceIDConvert(rownames(countdataBoth), "MSU", "RAP")

#get TF data
#read in TF information from CISBP; save GeneID and TF family name
#1537 -> 1444 TFs after conversion
tfsCISBP <- read.csv("TF_Information.csv") %>% 
  dplyr::select(DBID, Family_Name) %>% 
  distinct() %>% 
  mutate(RAP = str_replace(DBID, "OS", "Os"),
         RAP = str_replace(RAP, "G", "g")) %>% 
  merge(riceIDConv, by="RAP", all.x=T) %>% 
  dplyr::select(MSU, Family_Name, -RAP, -DBID) %>% 
  dplyr::rename("GeneID" = "MSU") %>% 
  filter(!is.na(GeneID)) %>% 
  distinct()
#read in TF information from TFDB, remove those from CISBP
#additional 597 TFs
tfsTFDB <- read.csv("OsTFs.csv") %>% 
  dplyr::rename("GeneID" = "Gene_ID") %>% 
  distinct() %>% 
  dplyr::rename("Family_Name" = "Family") %>% 
  filter(!GeneID %in% tfsCISBP$GeneID)
tfs <- rbind(tfsCISBP, tfsTFDB)

#match TFs to motifs (also from CisBP)
#TFs can be associated with multiple motifs (usually from the same TF family)
motifs <- read.csv("TF_Information_all_motifs_plus.csv") %>% 
  dplyr::select(Motif_ID, DBID, Family_Name) %>% 
  distinct() %>% 
  filter(Motif_ID != ".") %>% 
  mutate(RAP = str_replace(DBID, "OS", "Os"),
         RAP = str_replace(RAP, "G", "g")) %>%  
  merge(riceIDConv, by="RAP", all.x=T) %>% 
  mutate(MSU = ifelse(str_sub(RAP,1,3) == "LOC", RAP, MSU)) %>% 
  dplyr::select(-RAP, -DBID) %>% 
  dplyr::rename("GeneID" = "MSU") %>% 
  filter(!is.na(GeneID))

#filter counts to include DEGs and all TFs (DEG or not)
#GENIE3 can identify regulatory TFs with different patterns than their targets,
#so the TFs don't necessarily have to be differentially expressed
rawCounts <- countdataBoth %>%
  as.data.frame() %>%  
  dplyr::select(paste0("Wang.S", rep(2:11,each=4), ".", 1:4),
                paste0("Li.S", rep(1:7,each=2), ".", 1:2)) %>% 
  filter(rownames(.) %in% rownames(dds)) %>% 
  filter(rownames(.) %in% dgeListNetwork | rownames(.) %in% tfs$GeneID) 

##########################################
#get promoter sequences for FIMO and SEA tools from meme-suite
#can run in browser or on command line

#these tools will retrieve the sequences upstream of the TSS of each gene
#FIMO is used to find matches of motifs in individual promoters
#SEA is used to find enriched motifs in a group of sequences (the expression clusters)

#load genome and index
library(Rsamtools)
OsGenome <- FaFile("OsGenome.fa")

#read gff3 with gene annotations
library(ape)
gff <- read.gff("Osativa.gff3") %>% 
  filter(type == "gene") %>% 
  #mutate(GeneID = str_sub(attributes, -10, -1)) #if removing LOC_ from GeneID
  mutate(GeneID = str_sub(attributes, -14, -1)) 

#convert gff to GRanges
library(GenomicFeatures)
geneAnnotations <- makeGRangesFromDataFrame(gff, keep.extra.columns = T)

#this first part is for Fig 1F; enriched motifs in each cluster
#get promoter sequence for all detected genes for control
allGenePromoters <- geneAnnotations[geneAnnotations@elementMetadata$GeneID %in% rownames(dds),]
promoterSeq <- extractUpstreamSeqs(OsGenome,allGenePromoters,width=500)
promoterSeq@ranges@NAMES <- allGenePromoters$GeneID
#writeXStringSet(promoterSeq, "allGenePromoters.fna") #write to fasta file

#loop to get promoters of DEGs in each cluster
for(i in 1:length(unique(clusterDEG2$Cluster))){
  clust <- clusterDEG2 %>% 
    filter(Cluster == i)
  geneList <- geneAnnotations[geneAnnotations@elementMetadata$GeneID %in% unique(clust$Gene),]
  promoterSeq <- extractUpstreamSeqs(OsGenome,geneList,width=500)
  promoterSeq@ranges@NAMES <- geneList$GeneID
  #writeXStringSet(promoterSeq, paste0("2023_05_SEA/Cluster",i,".Promoters.fa")) #write to fasta file
}

#this part is for the network; identify motifs in the promoter region of each gene
#we'll use this to narrow down the list of potential regulators for each target gene
#get promoters of all DEGs and TFs
DEGandTFPromoters <- geneAnnotations[geneAnnotations@elementMetadata$GeneID %in% rownames(rawCounts),]
promoterSeq <- extractUpstreamSeqs(OsGenome,DEGandTFPromoters,width=500)
promoterSeq@ranges@NAMES <- DEGandTFPromoters$GeneID
#writeXStringSet(promoterSeq, "DEGandTFPromoters.fna") #write to fasta file


#################################
#final setup for network prediction
#use output from FIMO

#load fimo data for regulators
#creates a dataframe of target GeneIDs and the enriched motifs in their promoters
fimo <- read_tsv("2023.04.20.fimoAllGenes500bp.tsv") %>% 
  filter(sequence_name %in% rownames(rawCounts)) %>% 
  #may have to mess around with the qvalue; there's a sweet spot between a sparse and dense network
  filter(`q-value` < 0.11) %>% 
  dplyr::select(motif_id,sequence_name) %>% 
  dplyr::rename("Motif_ID" = "motif_id",
                "TargetGeneID" = "sequence_name") %>% 
  distinct()

#merge with the motifs dataframe made earlier
#TFs are matched to target genes based on motif
fimoTargets <- fimo %>% 
  merge(motifs, by = "Motif_ID") %>% 
  dplyr::rename("RegulatorGeneID" = "GeneID") %>% 
  filter(RegulatorGeneID %in% rownames(rawCounts)) %>% 
  filter(TargetGeneID %in% rownames(rawCounts))


#convert that dataframe into a list to be supplied to GENIE3
#GENIE3 can use custom predictors for each target
#takes 2-4 mins; print progress (could split() be used here instead to speed it up?)
target <- unique(fimoTargets$TargetGeneID)
regulatorsList <- vector(mode = "list", length = length(target))
for(i in 1:length(target)){
  temp <- fimoTargets[fimoTargets$TargetGeneID == target[i],4]
  regulatorsList[i] <- list(unique(temp))
  if(i%%1000==0){print(i)}
}
names(regulatorsList) <- target

#plot distribution of # of regulators per target
#too few will create a sparse network; too many will be dense
#an average of 20-50 works well
data.frame(X = lengths(regulatorsList)) %>% 
  ggplot(aes(X)) +
  geom_histogram()

#refilter raw counts based on remaining genes
#some will have been lost if they had no sig. motifs in their promoter or no motif-TF matches
#10,471 -> 10,010 genes after filtering
filteredCounts <- rawCounts[rownames(rawCounts) %in% names(regulatorsList) | 
                            rownames(rawCounts) %in% fimoTargets$RegulatorGeneID,] %>% 
  as.matrix()


################################
#GENIE3
library(GENIE3)
library(igraph)
library(RCy3) #this package lets you load the network directly into Cytoscape

#run GENIE3 using custom regulator list per gene
#can take a while, especially with large numbers of genes or conditions
#to speed it up, use 1 less core than you have available 
#(i.e. if you have 4 cores on your computer, use 3, or R will start to lock up)
genie <- GENIE3(filteredCounts, targets=names(regulatorsList), regulators=regulatorsList, nCores = 3)

#get a table of regulator-target pairs and the weight 
genieLinkList <- getLinkList(genie) 
#you have to mess around with weight filter; there is no standard/significant number
#again, it's a balance between a sparse and dense network
genieLink <- genieLinkList %>% 
  filter(weight > 0.173) 

#use this to check number of nodes and edges; don't try loading a network with over 10k genes
paste0("Network with weight threshold of ", round(min(genieLink$weight),3), ", comprised of ",
       length(unique(c(genieLink$regulatoryGene, genieLink$targetGene))), " nodes and ",
       nrow(genieLink), " edges.")


#convert filtered network table to igraph & load into Cytoscape (Cytoscape needs to be open)
graph <- graph.data.frame(genieLink,directed = T)
createNetworkFromIgraph(graph, "TFNetwork")


#add metadata to network, including: cluster from Fig 1D, TF family, outdegree, etc
DEGClusters <- clusterDEG2 %>% 
  dplyr::rename("GeneID" = "Gene") %>% 
  dplyr::select(GeneID, Cluster) %>% 
  distinct()
deg <- igraph::degree(graph, mode = "out")
tfSelect <- data.frame(GeneID = V(graph)$name) %>% 
  mutate(TF = ifelse(GeneID %in% fimoTargets$RegulatorGeneID, "Yes", "No"),
         Degree = as.integer(deg)) %>% 
  merge(tfs, by="GeneID", all.x=T) %>% 
  merge(DEGClusters, by="GeneID", all.x=T) %>% 
  #merge(groups, by = "GeneID", all.x=T) %>% 
  column_to_rownames("GeneID")
#load into Cytoscape
loadTableData(tfSelect)




#########################
#other functions for exploring network


#get TFs of each cluster
TfDEGs <- tfs %>% 
  filter(GeneID %in% dgeList)
TfClusters <- clusterDEG %>% 
  dplyr::select(Gene, NewCluster) %>% 
  dplyr::rename("GeneID" = "Gene") %>% 
  distinct() %>% 
  merge(TfDEGs, by = "GeneID")
TfClusterCount <- TfClusters %>% 
  group_by(NewCluster, Family_Name) %>% 
  summarise(Count = n()) %>% 
  pivot_wider(., names_from = NewCluster, values_from = Count)
TfClusterCount[is.na(TfClusterCount)] <- 0
TfClusterAll <- TfClusters %>%
  group_by(NewCluster) %>% 
  summarise(Count = n())


#group by module (select TFs and targets in CytoScape)
files <- list.files("FinalModules")
mod <- data.frame()
for(i in 1:length(files)){
  temp <- read.csv(paste0("FinalModules/",files[i])) %>% 
    filter(selected == "true" | selected == "TRUE") %>% 
    mutate(Group = files[i])
  mod <- rbind(mod, temp)
}

groups <- mod %>% 
  dplyr::select(name, Group) %>% 
  rename("name" = "GeneID") %>% 
  mutate(Dup = duplicated(GeneID) | duplicated(GeneID, fromLast = T)) %>% 
  mutate(Group = ifelse(duplicated(GeneID), "Multiple", Group)) %>% 
  filter(Dup == F | (Dup == T & Group == "Multiple")) %>% 
  dplyr::select(-Dup) %>% 
  distinct()

#####################
#GO terms of new modules

Res <- data.frame()
mods <- unique(mod$Group)
for(i in 1:5){
  
  #Read in list of genes from clusters
  genes <- mod %>% 
    filter(Group == mods[i])
  geneList <- factor(as.integer(geneNames %in% unique(genes$id)))
  names(geneList) <- geneNames
  
  for(h in 2){
    ont <- c("MF","BP", "CC")
    #GO function; change "ontology" to MF, BP, or CC
    GOdata <- new("topGOdata", ontology = ont[h], allGenes = geneList, 
                  annot = annFUN.gene2GO, gene2GO = GO2geneID, nodeSize = 5)
    
    #Run fisher test with default weighted algorithm and filter by pval
    resultFis <- runTest(GOdata, statistic = "fisher")
    sig <- data.frame(resultFis@score) %>% 
      filter(resultFis.score < 0.05)
    
    #Generate results table
    temp <- GenTable(GOdata, weightFisher = resultFis, orderBy = resultFis, topNodes = nrow(sig), numChar = 100)
    temp$Module <- mods[i]
    Res <- rbind(Res, temp)
  }
}
Res$weightFisher <- gsub("< ", "", Res$weightFisher)
allRes <- data.frame("Term" = Res$Term,
                     "pval" = as.numeric(Res$weightFisher),
                     "Module" = Res$Module)
allRes <- allRes[!duplicated(allRes[,1:2]),]












