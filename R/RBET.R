#' @title RBET.
#' @description Reference-based batch effect testing.
#'
#' @param obj A Seurat object of single cell data.
#' @param rgene A vector of reference genes.
#' @param batch The column name for batch information in obj.
#' @param k The sub-sampling size for RBET.
#'
#' @return The batch effect of dataset.
#' @export
#'
#' @import Seurat
#' @import uwot
#' @importFrom dplyr %>%

RBET <- function(obj, rgene, batch, k=50){
  obj$batch <- obj@meta.data[,batch]
  data <- GetAssayData(obj) %>% as.data.frame()

  hgene <- c()
  ngene <- c()
  for(gene in rgene){
    if(gene%in%rownames(data)){
      hgene <- append(hgene,gene)
    }else{
      ngene <- append(ngene,gene)
    }
  }
  if(length(ngene)==1){
    print(paste0("Reference gene ", ngene, " is not found."))
  }else if(length(ngene)>1){
    print(paste0("Reference genes ", paste0(ngene, collapse = ", "), " are not found."))
  }
  rgene <- hgene

  data <- data[rgene,] %>% t() %>% as.data.frame()
  data <- as.data.frame(uwot::umap(data))
  data$batch = obj$batch
  colnames(data) = c("umap1", "umap2", "batch")

  batchname <- unique(data$batch)
  batchnum <- length(batchname)
  res <- c()
  for(i in 1:(batchnum-1)){
    for(j in (i+1):batchnum){
      batch1 <- batchname[i]
      batch2 <- batchname[j]
      x <- data[, 1:2][which(data[, 'batch'] == batch1),]
      y <- data[, 1:2][which(data[, 'batch'] == batch2),]
      res <- append(res, MAC2(x, y, k))
    }
  }
  res_all <- mean(res)
  return(res_all)
}

