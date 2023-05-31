
#' Calculte differential expression genes
#'
#' @param counts count data
#' @param coldata meta information for sample
#' @param ident.1 group1
#' @param ident.2 group2
#' @param group group name
#'
#' @return deg result data.frame
#' @export
#'
#' @examples CalDEG(counts, coldata, ident.1 = "WT", ident.2 = "KO", group = 'groupType')
#'
#'
#'
CalDEG = function(counts, coldata, ident.1 = "WT", ident.2 = "KO", group = 'groupType',
                  core = 1){

  require(limma)
  require(future)
  plan("multicore", workers = core)

  print("create seurat obj...")

  seu = Seurat::CreateSeuratObject(counts = counts, meta.data = data.frame(coldata))

  print("calculate dge...")
  dge.cluster =  Seurat::FindMarkers(seu, ident.1 = ident.1, ident.2 = ident.2, group.by = group)
  dge.cluster = data.frame(dge.cluster)
  dge.cluster$gene = row.names(dge.cluster)

  return(dge.cluster)
}
