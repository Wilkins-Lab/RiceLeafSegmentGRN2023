library(tidyverse)
library(DESeq2)
library(GENIE3)
library(igraph)
library(Rgraphviz)
library(RCy3)

###########################
#remove S1 DEGs
DEGsFiltNoS1 <- DEGsFilt %>% 
  filter(str_sub(Comparison,1,8) != "Wang.S1_")
dgeListNetwork <- dgeList2[dgeList2 %in% DEGsFiltNoS1$GeneID]

###########################
#NETWORK
###########################
#MSU and RAP ID conversion
library(riceidconverter)
riceIDConv <- RiceIDConvert(rownames(countdataBoth), "MSU", "RAP")

#get TF data
#read in TF information from CISBP
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

#match TFs to motifs
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

#filter =counts to only include DEGs and all TFs (DEG or not)
rawCounts <- countdataBoth %>%
  as.data.frame() %>%  
  dplyr::select(paste0("Wang.S", rep(2:11,each=4), ".", 1:4),
                paste0("Li.S", rep(1:7,each=2), ".", 1:2)) %>% 
  filter(rownames(.) %in% rownames(dds)) %>% 
  filter(rownames(.) %in% dgeListNetwork | rownames(.) %in% tfs$GeneID) 

##########################################
#get promoter sequences for FIMO
#or enriched motifs in each DEG cluster

#load genome and index
library(Rsamtools)
OsGenome <- FaFile("OsGenome.fa")

#read gff3
library(ape)
gff <- read.gff("Osativa.gff3") %>% 
  filter(type == "gene") %>% 
  #mutate(GeneID = str_sub(attributes, -10, -1)) #if removing LOC_
  mutate(GeneID = str_sub(attributes, -14, -1)) 

#convert gff to GRanges
library(GenomicFeatures)
geneAnnotations <- makeGRangesFromDataFrame(gff, keep.extra.columns = T)

#all detected genes for control
allGenePromoters <- geneAnnotations[geneAnnotations@elementMetadata$GeneID %in% rownames(dds),]
promoterSeq <- extractUpstreamSeqs(OsGenome,allGenePromoters,width=500)
promoterSeq@ranges@NAMES <- allGenePromoters$GeneID
writeXStringSet(promoterSeq, "SEA5Clusters/allGenePromoters.fna")

#loop to get promoters for each cluster
for(i in 1:length(unique(clusterDEG2$Cluster))){
  clust <- clusterDEG2 %>% 
    filter(Cluster == i)
  geneList <- geneAnnotations[geneAnnotations@elementMetadata$GeneID %in% unique(clust$Gene),]
  promoterSeq <- extractUpstreamSeqs(OsGenome,geneList,width=500)
  promoterSeq@ranges@NAMES <- geneList$GeneID
  writeXStringSet(promoterSeq, paste0("2023_05_SEA/Cluster",i,".Promoters.fa"))
}

#all DEGs and tfs for FIMO
DEGandTFPromoters <- geneAnnotations[geneAnnotations@elementMetadata$GeneID %in% rownames(rawCounts),]
promoterSeq <- extractUpstreamSeqs(OsGenome,DEGandTFPromoters,width=500)
promoterSeq@ranges@NAMES <- DEGandTFPromoters$GeneID
writeXStringSet(promoterSeq, "2023.05.08.DEGandTFPromoters.fna")


#################################
#load fimo data for regulators
fimo <- read_tsv("2023.04.20.fimoAllGenes500bp.tsv") %>% 
  filter(sequence_name %in% rownames(rawCounts)) %>% 
  filter(`q-value` < 0.11) %>% 
  dplyr::select(motif_id,sequence_name) %>% 
  dplyr::rename("Motif_ID" = "motif_id",
                "TargetGeneID" = "sequence_name") %>% 
  distinct()

fimoTF <- fimo %>% 
  merge(motifs, by = "Motif_ID") %>% 
  dplyr::rename("RegulatorGeneID" = "GeneID") %>% 
  filter(RegulatorGeneID %in% rownames(rawCounts)) %>% 
  filter(TargetGeneID %in% rownames(rawCounts))

fimoTargets <- fimoTF


#find regulators for each target gene based on presence of motif
#takes 2-4 mins; print progress
target <- unique(fimoTargets$TargetGeneID)
regulatorsList <- list()
for(i in 1:length(target)){
  temp <- fimoTargets[fimoTargets$TargetGeneID == target[i],4]
  regulatorsList <- append(regulatorsList, list(unique(temp)))
  if(i%%1000==0){print(i)}
}
names(regulatorsList) <- target
#plot distribution of regulators per target
data.frame(X = lengths(regulatorsList)) %>% 
  ggplot(aes(X)) +
  geom_histogram()

#refilter
filteredCounts <- rawCounts[rownames(rawCounts) %in% names(regulatorsList) | rownames(rawCounts) %in% fimoTargets$RegulatorGeneID,] %>% 
  as.matrix()

genie <- GENIE3(filteredCounts, targets=names(regulatorsList), regulators=regulatorsList, nCores = 3)

genieLinkList <- getLinkList(genie) 
genieLink <- genieLinkList %>% 
  filter(weight > 0.173) 
length(unique(c(genieLink$regulatoryGene, genieLink$targetGene)))

#add data to graph
graph <- graph.data.frame(genieLink,directed = T)
createNetworkFromIgraph(graph, "TFNetwork")


DEGClusters <- clusterDEG2 %>% 
  dplyr::rename("GeneID" = "Gene") %>% 
  dplyr::select(GeneID, Cluster) %>% 
  distinct()

deg <- igraph::degree(graph, mode = "out")
tfSelect <- data.frame(GeneID = V(graph)$name) %>% 
  mutate(TF = ifelse(GeneID %in% fimoTF$RegulatorGeneID, "Yes", "No"),
         Degree = as.integer(deg)) %>% 
  merge(tfs, by="GeneID", all.x=T) %>% 
  merge(DEGClusters, by="GeneID", all.x=T) %>% 
  merge(groups, by = "GeneID", all.x=T) %>% 
  column_to_rownames("GeneID")
loadTableData(tfSelect)


#########################
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
write.csv(TfClusterCount, "2023.05.05.TfClusterCount.csv")


TfClusterAll <- TfClusters %>%
  group_by(NewCluster) %>% 
  summarise(Count = n())


#########################################
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

#write.csv(allRes, "2023.05.11.GRNModuleGoTerms.csv", row.names = F)











