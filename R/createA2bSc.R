#' Create a2bSc object
#'
#' @return a2bSc object
#' @export
#'
#' @examples createA2bSc()
#'
createA2bSc = function(){
  load("data/scSeq/06.AnnCluster/scAnn.rda")
  ## get expression and feature information
  obj = new("a2bsc")
  obj@Exp = sc@assays$RNA@counts
  obj@pData = sc@meta.data
  obj@fData = data.frame(gene = sc@assays$RNA@counts@Dimnames[[1]])

  ## get dimension reduction information
  obj@dimRe = new("DimRe")
  obj@dimRe@umap = sc@reductions$umap@cell.embeddings %>% data.frame()
  obj@dimRe@tnse = sc@reductions$tsne@cell.embeddings %>% data.frame()

  ## get annotation information
  ### main annotation
  obj@Ann@main@markers = read.table("data/scSeq/06.AnnCluster/markers.xls", header = T)
  obj@Ann@main@ann = data.frame(id = colnames(sc), cluster = sc@meta.data$seurat_clusters,
                                cellType = sc@active.ident,
                                group = sc@meta.data$orig.ident)

  ### T cell annotation
  load("data/scSeq/07.SubCluster/TcellAnn.rda")
  obj@Ann@Tcell@markers = read.table("data/scSeq/07.SubCluster/Tmarkers.xls", header = T)
  obj@Ann@Tcell@ann = data.frame(id = colnames(scT), cluster = scT@meta.data$seurat_clusters,
                                 cellType = scT@active.ident,
                                 group = sc@meta.data$orig.ident[
                                   match(colnames(scT), colnames(sc))])

  return(obj)
}
