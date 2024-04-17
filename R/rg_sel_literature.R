#' @title Literature-based reference genes selection.
#' @description Select reference genes based on the extent to which they are differently expressed across batches.
#'
#' @param obj A Seurat object with batch information.
#' @param batch The column name for batch information in obj.
#' @param celltype The column name for cell type information in obj.
#' @param rgene Housekeeping genes validated by previous literature.
#'
#' @return A vector of reference genes.
#' @export
#'
#' @import Seurat
#' @importFrom dplyr %>%
#' @importFrom RVAideMemoire mood.medtest

rg_sel_literature <- function(obj, batch, celltype, rgene){
  data <- GetAssayData(obj,layer="counts") %>% as.data.frame()
  batch_info <- obj@meta.data[,batch]
  celltype_info <- obj@meta.data[,celltype]

  hgene <- c()
  for(gene in rgene){
    if(gene%in%rownames(data)){
      hgene <- append(hgene,gene)
    }else{
      print(paste0("The requested gene ", gene, " is not found."))
    }
  }
  data <- data[hgene,]
  dat <- data.frame(data=t(data),batch=batch_info,type=celltype_info)
  info_batch <- unique(batch_info)
  info_celltype <- unique(celltype_info)

  batch_median <- matrix(NA, nrow=length(hgene),ncol=length(info_celltype))
  rownames(batch_median) <- hgene
  colnames(batch_median) <- info_celltype

  for(gene in 1:length(hgene)){
    for(type in 1:length(info_celltype)){
      response <- c()
      fact <- c()
      for(bat in 1:length(info_batch)){
        cur_gene <- dat[dat$batch==info_batch[bat]&dat$type==info_celltype[type],gene]
        cur_gene <- log(cur_gene+1)
        response <- append(response,cur_gene)
        fact <- append(fact,c(rep(info_batch[bat],length(cur_gene))))
      }
      if(length(unique(response))>1){
        fact <- as.factor(fact)
        test <- mood.medtest(response~fact)
        batch_median[gene,type] <- test$p.value
      }
    }
  }
  median_p <- batch_median<=0.05
  median_p[is.na(median_p)] <- FALSE
  median_p <- as.data.frame(median_p)
  median_p$count <- apply(median_p,1,sum)

  rg <- c()
  for(num in sort(unique(median_p$count),decreasing = TRUE)){
    rg <- append(rg, rownames(median_p)[median_p$count==num])
    if(length(rg)>=2){
      break
    }
  }
  return(rg)
}

