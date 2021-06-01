library(fgsea)

load("/cloud/project/kegg_definitions_list.RData")
fgsea(
  pathways = kegg_definitions_grouped_list,
  stats    = gene_exp_vector,
  # minSize  = 2,
  # maxSize  = 500,
  eps = 0.0
)


