library(openxlsx)
library(Rcpp)
library(umap)
library(makedummies)
library(Polychrome)
library(dplyr)
library(ggplot2)
library(dbscan)


## Loading representative tss data -----------------
Rep_tss_path <- "source_data/representative_tss.txt"
Rep_tss <- read.table(Rep_tss_path, stringsAsFactors=F)
Rep_tss <- Rep_tss[-which(Rep_tss[,1] == "chrM"),]


## Loading core promoter elements--------------------------
CPE <- read.table("source_data/CPE.txt", stringsAsFactors=F)

## Loading DPPs ---------------------------
# DPSs of genes meet 1st energetic criteria
DPPs_E1 <- read.xlsx("source_data/DPPs_1th.xlsx", rowNames=T)
DPPs_E2 <- read.xlsx("source_data/DPPs_2th.xlsx", rowNames=T)
DPPs_E3 <- read.xlsx("source_data/DPPs_3th.xlsx", rowNames=T)

# conbine the all data
data_all <- rbind(DPPs_E1, DPPs_E2, DPPs_E3)

# Obtain core less gene, TATA box, and Inr containg genes
core_less <- names(which(apply(CPE, 1, sum) == 0 ) )
TATA_box <- rownames(CPE[which(CPE[,"TATA_8mer"] == 1),])
Inr <- rownames(CPE[which(CPE[,"Inr_6mer"] == 1),])


## make a binary matrix----------------
count_matrix <- data_all[,grep("count",
colnames(data_all))]
rownames(count_matrix) <- rownames(data_all)
# converting each data
for(i_row in 1:nrow(count_matrix)){
  for(i_col in 1:ncol(count_matrix)){
    if(count_matrix[i_row, i_col] == "YES"){
      count_matrix[i_row, i_col] <- 1
    }else{
      count_matrix[i_row, i_col] <- 0
    }
  }
}


## Make rank features of energy ----------------
hit_1st <- hit_2nd <- rep(0, nrow(count_matrix))
names(hit_1st) <- names(hit_2nd) <- rownames(count_matrix)
hit_1st[which(rownames(count_matrix) %in% rownames(DPPs_E1))] <- 1
hit_2nd[which(rownames(count_matrix) %in% rownames(DPPs_E2))] <- 1

# combine the data
data_comb_dumm_for_umap <- cbind(count_matrix, hit_1st, hit_2nd)
data_comb_dumm_for_umap_Coreless <- cbind(count_matrix[which(rownames(count_matrix) %in% core_less),],
hit_1st[which(names(hit_1st) %in% core_less)],
hit_2nd[which(names(hit_2nd) %in% core_less)])


## UMAP --------------------
# set up parameters
umap_params <- umap.defaults
umap_params$n_neighbors <- 15
umap_params$random_state <- 42


## Converting data type ------
# all genes
data_num <- t(apply(data_comb_dumm_for_umap, 1,
  function(x) as.numeric(x)))

# core-less genes
data_num_CL <- t(apply(data_comb_dumm_for_umap_Coreless, 1,
  function(x) as.numeric(x)))

## removing NA data in columns-------
# all genes
rem_col <- which(apply(data_num, 2, sum) == 0)
if(length(rem_col)){
  data_num <- data_num[,-rem_col]
}

# core-less genes
rem_col_CL <- which(apply(data_num_CL, 2, sum) == 0)
if(length(rem_col_CL)){
  data_num_CL <- data_num_CL[,-rem_col_CL]
}


## removing NA data in rows-------
# all data
rem_row <- which(apply(data_num, 1, sum) == 0)
if(length(rem_row)){
  data_num <- data_num[-rem_row,]
}

# core-less gene data
rem_row_CL <- which(apply(data_num_CL, 1, sum) == 0)
if(length(rem_row_CL)){
  data_num_CL <- data_num_CL[-rem_row_CL,]
}

## perform UMAP ------
# all genes
ump_out_comb <- umap(data_num ,config=umap_params)
Data5_com <- ump_out_comb$layout

# core-less genes
ump_out_comb_CL <- umap(data_num_CL ,config=umap_params)
Data5_com_CL <- ump_out_comb_CL$layout

## perform DBSCAN -------------------------------
# all genes
data_2 <- data.frame(Data5_com[,1],Data5_com[,2])
res <- dbscan(x = data_2, eps = 1.5, minPts = 5)
data_dbscan <- data_2 %>% mutate(cluster = res$cluster)
data_dbscan$cluster <- factor(data_dbscan$cluster)

# core-less genes
data_2_CL <- data.frame(Data5_com_CL[,1],Data5_com_CL[,2])
res_CL <- dbscan(x = data_2_CL, eps = 1.5, minPts = 5)
pos_3 <- which(res_CL$cluster== 3); pos_1 <- which(res_CL$cluster== 1); pos_2 <- which(res_CL$cluster== 2)
res_CL$cluster[pos_3] <- 1; res_CL$cluster[pos_1] <- 2; res_CL$cluster[pos_2] <- 3
data_dbscan_CL <- data_2_CL %>% mutate(cluster_CL = res_CL$cluster)
data_dbscan_CL$cluster_CL <- factor(data_dbscan_CL$cluster_CL)


## plot the data --------------------------------

# plot all data

ggplot(data = data_dbscan, mapping = aes(x = Data5_com[,1], y = Data5_com[,2] ,colour = cluster),family ="Helvetica") +
    geom_point(size=0.3,alpha = 0.3) +   # 散布図に可視化
    xlab("UMAP coordinate 1") +
    ylab("UMAP coordinate 2") +
    theme(
      panel.background = element_rect(fill = "white", colour = "black", size = 1.0),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      legend.key = element_blank(),
      text = element_text(size = 15),
      #axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
    ggsave("umap_all_genes.pdf",  width = 4, height =4 ,device = "pdf")

# save the data
write.xlsx(data_dbscan, file="umap_all.xlsx", rowNames =T)

# plot TATA-box ------------
ggplot(data = data.frame(Data5_com[intersect(rownames(Data5_com) ,TATA_box),]),
mapping = aes(x = Data5_com[intersect(rownames(Data5_com) ,TATA_box),1],y = Data5_com[intersect(rownames(Data5_com) ,TATA_box),2] ),
family ="Helvetica") +
geom_point(size=1,alpha = 0.3,shape=2) +
xlab("UMAP coordinate 1") +
ylab("UMAP coordinate 2") +
theme(
  panel.background = element_rect(fill = "white", colour = "black", size = 1.0),
  panel.grid = element_blank(),
  strip.background = element_blank(),
  legend.key = element_blank(),
  text = element_text(size = 15),
  #axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank()
)
# save
ggsave("umap_tata.pdf",  width = 4, height =4 ,device = "pdf")
write.xlsx(data.frame(Data5_com[intersect(rownames(Data5_com) ,TATA_box),]), file="umap_tata.xlsx", rowNames =T)

# plot INR ------------------
ggplot(data = data.frame(Data5_com[intersect(rownames(Data5_com) ,Inr),]),
mapping = aes(x = Data5_com[intersect(rownames(Data5_com) ,Inr),1], y = Data5_com[intersect(rownames(Data5_com) ,Inr),2] ),
family ="Helvetica") +
geom_point(size=0.3,alpha = 0.3,shape=3) +
xlab("UMAP coordinate 1") +
ylab("UMAP coordinate 2") +
theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.0),
panel.grid = element_blank(),
strip.background = element_blank(),
legend.key = element_blank(),
text = element_text(size = 15),
#axis.title = element_blank(),
axis.text = element_blank(),
axis.ticks = element_blank()
)
# save
ggsave("umap_inr.pdf",  width = 4, height =4 ,device = "pdf")
write.xlsx(data.frame(Data5_com[intersect(rownames(Data5_com) ,Inr),]), file="umap_inr.xlsx", rowNames =T)

# plot core-less --------------
ggplot(data = data_dbscan_CL, mapping = aes(x = Data5_com_CL[,1], y = Data5_com_CL[,2] ,colour = cluster_CL),family ="Helvetica") +
    geom_point(size=0.3,alpha = 0.3) +   # 散布図に可視化
    xlab("UMAP coordinate 1") +
    ylab("UMAP coordinate 2") +
    theme(
      panel.background = element_rect(fill = "white", colour = "black", size = 1.0),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      legend.key = element_blank(),
      text = element_text(size = 15),
      #axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )

# save
ggsave("umap_core-less_genes.pdf",  width = 4, height =4 ,device = "pdf")
write.xlsx(data_dbscan_CL, file="umap_core-less.xlsx", rowNames =T)

## saving data ----------------------------------
# all data
for(i_write in 1:max(as.numeric(data_dbscan$cluster))){
  genes_each <- sapply(rownames(data_dbscan)[which(as.character(data_dbscan$cluster)==i_write)], function(x) unlist(strsplit(x, "\\.|\\|"))[1])
  write.table(genes_each, file=sprintf("gene_all%d.txt", i_write), col.names=F, row.names=F, quote=F)
}

# core-less data
for(i_write in 1:max(as.numeric(data_dbscan_CL$cluster_CL))){
  genes_each <- sapply(rownames(data_dbscan_CL)[which(as.character(data_dbscan_CL$cluster_CL)==i_write)], function(x) unlist(strsplit(x, "\\.|\\|"))[1])
  write.table(genes_each, file=sprintf("gene_coreless%d.txt", i_write), col.names=F, row.names=F, quote=F)
}


## for fig. 5g calculation ----------------------

# cluster 1 ----
cluster1_E1st <- intersect(names(which(hit_1st==1)), rownames(data_dbscan[which(data_dbscan$cluster==1),]))
# length(cluster1_E1st) >>> 1296

cluster1_E2nd <- intersect(names(which(hit_2nd==1)), rownames(data_dbscan[which(data_dbscan$cluster==1),]))
# length(cluster1_E2nd) >>> 369

cluster1_E3rd  <- intersect(names(which((hit_2nd!=1)&(hit_1st!=1))), rownames(data_dbscan[which(data_dbscan$cluster==1),]))
# length(cluster1_E3rd) >>> 179


# cluster 2 ----
cluster2_E1st <- intersect(names(which(hit_1st==1)), rownames(data_dbscan[which(data_dbscan$cluster==2),]))
# length(cluster2_E1st) >>> 9283

cluster2_E2nd <- intersect(names(which(hit_2nd==1)), rownames(data_dbscan[which(data_dbscan$cluster==2),]))
# length(cluster2_E2nd) >>> 16

cluster2_E3rd  <- intersect(names(which((hit_2nd!=1)&(hit_1st!=1))), rownames(data_dbscan[which(data_dbscan$cluster==2),]))
# length(cluster2_E3rd) >>> 978


# cluster 3 ---- denote as cluster 4 in the figures
cluster3_E1st <- intersect(names(which(hit_1st==1)), rownames(data_dbscan[which(data_dbscan$cluster==3),]))
# length(cluster3_E1st) >>> 829

cluster3_E2nd <- intersect(names(which(hit_2nd==1)), rownames(data_dbscan[which(data_dbscan$cluster==3),]))
# length(cluster3_E2nd) >>> 194

cluster3_E3rd  <- intersect(names(which((hit_2nd!=1)&(hit_1st!=1))), rownames(data_dbscan[which(data_dbscan$cluster==3),]))
# length(cluster3_E3rd) >>> 85

# cluster 4 ---- denote as cluster 3 in the figures
cluster4_E1st <- intersect(names(which(hit_1st==1)), rownames(data_dbscan[which(data_dbscan$cluster==4),]))
# length(cluster4_E1st) >>> 4

cluster4_E2nd <- intersect(names(which(hit_2nd==1)), rownames(data_dbscan[which(data_dbscan$cluster==4),]))
# length(cluster4_E2nd) >>> 2352

cluster4_E3rd  <- intersect(names(which((hit_2nd!=1)&(hit_1st!=1))), rownames(data_dbscan[which(data_dbscan$cluster==4),]))
# length(cluster4_E3rd) >>> 19
