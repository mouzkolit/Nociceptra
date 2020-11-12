  library("biomaRt")
  nociceptor_expression <- read.csv("/home/physiologie/Dropbox/iPSC/Notebooks/Supercluster_iPSC_network_mRNA_all_zelllines.csv")
  mart_hsa <- useEnsembl(biomart="ensembl", dataset='hsapiens_gene_ensembl')
  entrezgene <- getBM(mart=mart_hsa,
                          attributes=c("entrezgene_id","external_gene_name"), 
                          filters="external_gene_name", 
                          values=nociceptor_expression$external_gene_name)
  
  gesamt.data <- merge(nociceptor_expression,entrezgene)
  library("DOSE")
  library("enrichplot")
  library("org.Hs.eg.db")

# provide an empty dataframe that can be used to cbind all the queried diseases 

data <- data.frame(
  ticker=character(),
  value=numeric(),
  date = as.Date(character()),
  stringsAsFactors=FALSE
)


for (i in 1:7){
  test <- gesamt.data[gesamt.data$hierachical_cluster == i,]
  disease <- enrichDO(test$entrezgene_id, ont = "DO", pvalueCutoff = 0.001, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 10000, qvalueCutoff = 0.2, readable = FALSE)
  cluster <- rep(i,nrow(disease@result))
  gesamt <- disease@result
  gesamt$cluster <- cluster
  data <- rbind(data, gesamt)
}


# make it for each individual cluster to detect potential codriving genes that drive diseases

color <- data.frame(
  ticker=character(),
  value=numeric(),
  date = as.Date(character()),
  stringsAsFactors=FALSE
)


for (t in unique(gesamt.data$cluster)){
  print(t)
  test_color <- gesamt.data[gesamt.data$cluster == t,]
  disease_color <- enrichDO(test_color$entrezgene_id, ont = "DO", pvalueCutoff = 0.001, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 10000, qvalueCutoff = 0.2, readable = FALSE)
  cluster_color <- rep(t,nrow(disease_color@result))
  color_enr <- disease_color@result
  color_enr$cluster <- cluster_color
  color <- rbind(color, color_enr)
}


for (i in 1:7){
  enrichment_plot <- gesamt.data[gesamt.data$hierachical_cluster == i,]
  disease1 <- enrichDGN(test$entrezgene_id)
  
  
  edox <- setReadable(disease1, 'org.Hs.eg.db', 'ENTREZID')
  p1<-heatplot(edox)
  plot(p1)
}

data_end <- data[data$p.adjust<= 0.05,]
color_end <- color[color$p.adjust <= 0.05,]

write.csv(data_end,"Supercluster_disease_enrichment.csv")
write.csv(color_end,"Color_cluster_disease_enrichment.csv")