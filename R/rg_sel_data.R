#' @title Data-based reference genes selection.
#' @description Automatically select reference genes from data, based on the extent to which they are stably expressed in batches and differently expressed across batches.
#'
#' @param obj A Seurat object with batch information.
#' @param batch The column name for batch information in obj.
#' @param dim The number of PCA used in RunUMAP and FindNeighbors.
#' @param resolution The resolution of clusters in FindClusters.
#' @param p The top p genes are retained when selecting invariable genes in clusters and across clusters, with p ranging from 0 to 1.
#' @param k The top k genes are retained when selecting the most differently expressed genes across batches.
#'
#' @return A vector of reference genes.
#' @export
#'
#' @import Seurat
#' @importFrom dplyr %>% group_by summarise
#' @importFrom RVAideMemoire mood.medtest
#' @importFrom stats na.omit sd wilcox.test

rg_sel_data <- function(obj, batch, dim=20, resolution=0.2, p=0.02, k=50){
  obj$batch <- obj@meta.data[,batch]
  obj.list <- SplitObject(obj, split.by = "batch")
  batch_name <- unique(obj$batch)

  # 1. clustering
  print("Clustering for each batch")
  seu_subset.list <- list()
  for(bat in 1:length(batch_name)){
    seu <- obj.list[[bat]]
    seu <- seu %>%
      NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
      ScaleData(features = rownames(seu))
    seu <- seu %>%
      RunPCA(features = VariableFeatures(object=seu)) %>%
      RunUMAP(reduction = "pca", dims = 1:dim) %>%
      FindNeighbors(dims = 1:dim) %>%
      FindClusters(resolution = resolution)
    seu$cluster <- seu@meta.data[,paste0(DefaultAssay(obj),"_snn_res.",resolution)]
    seu_subset.list[[batch_name[bat]]] <- seu
  }

  # 2. stably expressed in clusters
  print("Computing stably expressed genes in clusters")
  gene <- rownames(GetAssayData(obj,layer="counts"))
  gene_count <- matrix(0, nrow = length(gene), ncol = length(batch_name)) # whether each gene is invariable in each batch
  rownames(gene_count) <- gene
  colnames(gene_count) <- batch_name
  for(bat in 1:length(batch_name)){
    seu_subset <- seu_subset.list[[batch_name[bat]]]
    num_clus <- length(unique(seu_subset$cluster))
    data_all <- GetAssayData(seu_subset,layer="counts")
    gene_res <- matrix(0, ncol = num_clus, nrow = nrow(data_all)) # in the current batch, whether each gene is invariable in each cluster
    rownames(gene_res) <- rownames(data_all)
    for(i in 1:num_clus){
      seu_tmp <- subset(seu_subset,subset=cluster==(i-1))
      data <- as.matrix(GetAssayData(seu_tmp,layer="counts"))
      clus_mean <- apply(data,1,mean)
      clus_sd <- apply(data,1,sd)
      clus_div <- clus_sd/clus_mean
      gene_rank <- rownames(data)[order(clus_div)]

      N = floor(nrow(data)*p)
      gene_sel <- gene_rank[1:N] # invariable genes in the current cluster
      gene_res[rownames(data_all)%in%gene_sel,i] <- 1
    }
    gene_res <- as.data.frame(gene_res)
    gene_res$count <- apply(gene_res,1,sum)

    gene_batch <- c() # invariable genes in the current batch
    N = floor(nrow(data_all)*p)
    for(num in seq(num_clus,1,-1)){
      gene_batch <- append(gene_batch, rownames(gene_res[gene_res$count==num,]))
      if(length(gene_batch)>=N){
        break
      }
    }
    gene_count[gene%in%gene_batch,bat] <- 1
  }
  gene_count <- as.data.frame(gene_count)
  gene_count$count <- apply(gene_count,1,sum)
  gene_final <- c() # the final invariable genes
  N = floor(length(gene)*p)
  for(num in c(length(batch_name),length(batch_name)-1)){ # tolerance: one batch
    gene_final <- append(gene_final, rownames(gene_count[gene_count$count==num,]))
    if(length(gene_final)>=N){
      break
    }
  }
  gene_in <- gene_final

  # 3. stably expressed across clusters
  print("Computing stably expressed genes across clusters")
  gene_count <- matrix(0, nrow = length(gene), ncol = length(batch_name))  # whether each gene is invariable in each batch
  rownames(gene_count) <- gene
  colnames(gene_count) <- batch_name
  for(bat in 1:length(batch_name)){
    seu_subset <- seu_subset.list[[batch_name[bat]]]
    data <- GetAssayData(seu_subset,layer="counts")
    data <- as.data.frame(data)
    gene_name <- rownames(data)
    num_clus <- length(unique(seu_subset$cluster))

    SI <- matrix(NA,nrow = length(gene_name),ncol = 3)
    colnames(SI) <- c("mean","sd","SI")
    rownames(SI) <- gene_name
    for(i in 1:length(gene_name)){
      dat <- data.frame(exp=t(data[i,]),cluster=seu_subset$cluster)
      colnames(dat) <- c("exp","cluster")
      dat_mean <- dat %>%
        group_by(cluster) %>%
        summarise(mean=mean(exp))
      SI[i,1] <- mean(dat$exp)
      SI[i,2] <- sd(dat_mean$mean)
    }
    SI[,3] <- SI[,2]/SI[,1]
    gene_res <- data.frame(V1=SI[,3])
    gene_rank <- rownames(gene_res)[order(gene_res$V1)]


    N = floor(nrow(data)*p)
    gene_sel <- gene_rank[1:N] # invariable genes in the current batch
    gene_count[gene%in%gene_sel,bat] <- 1
  }
  gene_count <- as.data.frame(gene_count)
  gene_count$count <- apply(gene_count,1,sum)
  gene_final <- c()  # the final invariable genes
  N = floor(length(gene)*p)
  for(num in c(length(batch_name),length(batch_name)-1)){
    gene_final <- append(gene_final, rownames(gene_count[gene_count$count==num,]))
    if(length(gene_final)>=N){
      break
    }
  }
  gene_across <- gene_final

  # 4. take intersection
  candidate_gene <- intersect(gene_across, gene_in)
  candidate_gene <- candidate_gene[startsWith(candidate_gene,"MT")==FALSE]

  # 5. differently expressed across batches
  print("Computing differentially expressed genes across batches")
  batch_num <- length(batch_name)
  batch_median <- matrix(0, nrow=length(candidate_gene),ncol=1)
  rownames(batch_median) <- candidate_gene
  for(i in 1:(batch_num-1)){
    for(j in (i+1):batch_num){
      batch1 <- batch_name[i]
      batch2 <- batch_name[j]

      seu_subset <- seu_subset.list[[batch1]]
      data <- as.data.frame(GetAssayData(seu_subset,layer="counts"))
      data <- data[candidate_gene,]
      dat1 <- data.frame(data=t(data),type=seu_subset$cluster)
      num1 <- length(unique(seu_subset$cluster))

      seu_subset <- seu_subset.list[[batch2]]
      data <- as.data.frame(GetAssayData(seu_subset,layer="counts"))
      data <- data[candidate_gene,]
      dat2 <- data.frame(data=t(data),type=seu_subset$cluster)
      num2 <- length(unique(seu_subset$cluster))

      for(gene in 1:length(candidate_gene)){
        for(n1 in 1:num1){
          for(n2 in 1:num2){
            cur_gene1 <- dat1[dat1$type==(n1-1),gene]
            cur_gene2 <- dat2[dat2$type==(n2-1),gene]
            response <- c(cur_gene1,cur_gene2)
            fact <- c(rep(batch1,length(cur_gene1)),rep(batch2,length(cur_gene2)))
            med_test <- data.frame(response=response,fact=fact)
            med_test <- na.omit(med_test)
            if(nrow(med_test[med_test$fact==batch1,])>30 & nrow(med_test[med_test$fact==batch2,])>30
               & length(unique(med_test$response))>1){
              med_test$fact <- as.factor(med_test$fact)
              test <- mood.medtest(response~fact,data=med_test)
              if(test$p.value<=0.05){
                batch_median[gene,] <- batch_median[gene,]+1
              }
            }
          }
        }
      }
    }
  }
  median_p <- as.data.frame(batch_median)

  # 6. select the top k genes as the final reference genes
  rg <- c()
  for(num in sort(unique(median_p$V1),decreasing = TRUE)){
    rg <- append(rg, rownames(median_p)[median_p$V1==num])
    if(length(rg)>=k){
      break
    }
  }
  return(rg)
}

