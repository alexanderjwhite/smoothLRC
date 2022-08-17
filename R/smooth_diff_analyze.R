smooth_diff_analyze <- function(input, labels = NULL, method = "t-test", assay_type = "logcounts", alpha = 0.05){


  if(method == "t-test"){
    differ <- t.test
  } else if (method == "wilcox"){
    differ <- wilcox.test
  } else {
    stop("Please select an appropriate method.")
  }

  smoothlrc_labels <- input@colData$smooth_cluster
  if(is.null(smoothlrc_labels) & is.null(labels)){
    stop("Labels are not provided and smoothLRC cluster labels do not exist.")
  } else if(is.null(labels)){
    labels <- smoothlrc_labels
  }
  feature_ind <- NULL
  unique_labels <- unique(labels)
  data <- SummarizedExperiment::assay(input, assay_type)
  for (i in unique_labels){
    print(i)
    clust_ind <- which(labels==i)
    nonclust_ind <- which(labels!=i)
    for (j in 1:nrow(data)){
      print(j)
      x <- data[j,clust_ind]
      y <- data[j,nonclust_ind]
      comparison <- differ(x,y)
      pval <- ifelse(is.nan(comparison$p.value),1,comparison$p.value)
      if(pval <= alpha){
        feature_ind <- c(feature_ind,j)
      }

    }
  }

}


design <- model.matrix(~labels, data = data_f)
fit <- lmFit(data_f, design)
contrast.matrix <- makeContrasts("labels", levels = design)
fit2C <- contrasts.fit(fit, contrast.matrix)
fit2C <- eBayes(fit2C)
topTable(fit2C)
