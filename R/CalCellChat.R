#' CellChat analysis
#'
#' @param data.input expression matrix
#' @param meta meta information
#' @param group choose the group of sample
#' @param db choose database: mouse:CellChatDB.mouse;
#' @param core multi-process
#'
#' @return cellchat object
#' @export
#'
#' @examples
#'
CellChatAnalysis = function(data.input, meta, group="cellType",
                            db = CellChatDB.mouse, core = 1){

  require(CellChat)
  options(stringsAsFactors = FALSE)

  data.input = normalizeData(data.input)
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = group)

  # Add cell information into meta slot of the object (Optional)
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = group) # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

  # Set the ligand-receptor interaction database
  CellChatDB <- db # use CellChatDB.mouse if running on mouse data
  # showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChatDB

  # Preprocessing the expression data for cell-cell communication analysis

  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multicore", workers = core) # do parallel
  #> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
  #> explicitly specify either 'multisession' or 'multicore'. In the current R
  #> session, 'multiprocess' equals 'multisession'.
  #> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
  #> processing ('multicore') is not supported when running R from RStudio
  #> because it is considered unstable. For more details, how to control forked
  #> processing or not, and how to silence this warning in future R sessions, see ?
  #> parallelly::supportsMulticore
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  # cellchat <- projectData(cellchat, PPI.human)

  # Compute the communication probability and infer cellular communication network

  cellchat <- computeCommunProb(cellchat)
  #> triMean is used for calculating the average gene expression per cell group.
  #> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
  #> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  # Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)


  # Calculate the aggregated cell-cell communication network

  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}
