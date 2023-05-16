
#pick up from Robertson2023NetworkAnalysis after loading FIMO data

#MYB or Gbox filters
motifsSub <- motifs %>% 
  filter(grepl("Myb/SANT", Family_Name))
  #filter(grepl("bZIP|bHLH", Family_Name))

fimoTF <- fimo %>% 
  merge(motifs, by = "Motif_ID") %>% 
  dplyr::rename("RegulatorGeneID" = "GeneID") %>% 
  filter(RegulatorGeneID %in% rownames(rawCounts)) %>% 
  filter(TargetGeneID %in% rownames(rawCounts)) %>% 
  filter(RegulatorGeneID %in% motifsSub$GeneID | TargetGeneID %in% motifsSub$GeneID | Motif_ID %in% motifsSub$Motif_ID)


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
  filter(weight > 0.14) #0.18 for Gbox #0.14 for Myb
length(unique(c(genieLink$regulatoryGene, genieLink$targetGene)))


#add data to graph
graph <- graph.data.frame(genieLink,directed = T)
createNetworkFromIgraph(graph, "MybNetwork")

DEGClusters <- clusterDEG %>% 
  dplyr::rename("GeneID" = "Gene",
                "Cluster" = "NewCluster") %>% 
  dplyr::select(GeneID, Cluster) %>% 
  distinct()

deg <- igraph::degree(graph, mode = "out")
tfSelect <- data.frame(GeneID = V(graph)$name) %>% 
  mutate(TF = ifelse(GeneID %in% fimoTF$RegulatorGeneID, "Yes", "No"),
         Degree = as.integer(deg)) %>% 
  merge(tfs, by="GeneID", all.x=T) %>% 
  merge(DEGClusters, by="GeneID", all.x=T) %>% 
  #add old modules
  merge(groups, by = "GeneID", all.x=T) %>% 
  column_to_rownames("GeneID")
loadTableData(tfSelect)



#########################################
#group by module
#2023GboxModules or 2023MybModules
files <- list.files("2023GboxModules")
mod <- data.frame()
for(i in 1:length(files)){
  temp <- read.csv(paste0("2023GboxModules/",files[i])) %>% 
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
for(i in 1:3){
  
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








