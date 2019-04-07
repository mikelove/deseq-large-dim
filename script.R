library(DESeq2)

# this will eventually become a special DESeq() call...
deseq_large_dim <- function(dds) {
  dds <- estimateSizeFactors(dds)
  dds <- getBaseMeansAndVariances(dds)
  dispersions(dds) <- DESeq2:::roughDispEstimate(y=counts(dds,normalized=TRUE),
                                                 x=model.matrix(design(dds),colData(dds)))
  dispersions(dds) <- pmax(dispersions(dds), 1e-8)
  res <- DESeqResults(mcols(dds)[,"baseMean",drop=FALSE])
  ape <- lfcShrink(dds, coef=2, type="apeglm", res=res, svalue=TRUE,
                   apeAdapt=FALSE, no.shrink=TRUE, returnList=TRUE)

  # TODO:
  # populate the mcols of dds with LFCs, SEs, p-values, etc.

  # TODO:
  # not clear to me yet how the user should interact with the 'dds' to get results
  # a generic results() call is not acceptable for this case
  # maybe a special type of results() call will be needed
  
  out <- list(log2FoldChange=ape$fit$map,
              lfcSE=ape$fit$sd,
              pvalue=2 * ape$fit$fsr)
  out
}

dmr <- function(x) exp(rnorm(length(x),log(0.05),0.5))
dds <- makeExampleDESeqDataSet(n=100,m=500,dispMeanRel=dmr)
fit <- deseq_large_dim(dds)
