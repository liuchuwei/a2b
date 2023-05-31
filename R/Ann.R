#' Get summary annotation information
#'
#' @param ann annotation information
#'
#' @return summary annotation information
#' @export
#'
#' @examples GetAnnSumm(ann)
#'
GetAnnSumm = function(ann){
  cluster_summ = ann[,c(2, 3)]
  cluster_summ = ann[!duplicated(main_ann_sum$cluster),] %>% arrange(cluster)

  sample_summ = table(ann$cellType, ann$group) %>% data.frame()
  sample_summ = reshape2::dcast(sample_summ, Var1 ~ Var2)

  sample_summ = table(ann$cellType, ann$group) %>% data.frame()
  sample_summ = reshape2::dcast(sample_summ, Var1 ~ Var2)

  raw_summ = table(ann$cluster, ann$group) %>% data.frame()
  raw_summ = reshape2::dcast(raw_summ, Var1 ~ Var2)

  raw_summ = table(ann$cluster, ann$group) %>% data.frame()
  raw_summ = reshape2::dcast(raw_summ, Var1 ~ Var2)

  out = list(cluster_summ = cluster_summ,
             sample_summ = sample_summ,
             raw_summ = raw_summ)
}
