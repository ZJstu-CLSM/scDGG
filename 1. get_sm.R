suppressMessages(library(parallelDist))
suppressMessages(library(SNFtool)) # SNF;spectralClustering
suppressMessages(library(GSEABase)) # getGmt, load pathway information
suppressMessages(library(AUCell)) # AUCell, pathway scoring method
suppressMessages(library(Seurat))
library(corpcor)


#Load Gene Matrix
load_matrix_for_GSE <- function(path){
  expr_matrix = read.csv(path,check.names=FALSE, row.names=1)
  expr_matrix = as.matrix(expr_matrix)
  return(expr_matrix)
}


# load pathway
load_pathway <- function(path,name){
  gSet=getGmt(paste(path,name,sep='\\'))
  return(gSet)
}


get_mat_path <- function(mat_gene, paPath, human_pathway_name) {
  gSet = load_pathway(paPath, paste(human_pathway_name,'.gmt',sep=''))
  return (gSet)
}


get_gene_coexpression_matrix <- function(mat_gene) {
  # Standardized gene expression data
  mat_gene_scaled = scale(mat_gene)

  # Calculation of Pearson correlation coefficient matrix between genes
  cor_matrix = cor(mat_gene_scaled, method = "pearson")

  print(dim(cor_matrix))

  return(cor_matrix)
}


get_coexpression_matrix <- function(expr_matrix, pathway_genes) {

  # Screening for genes within pathways
  if (is.list(pathway_genes)) {
    pathway_genes <- unlist(pathway_genes)
  }
  pathway_genes <- as.character(pathway_genes)
  print(pathway_genes)
  # Remove potential spaces and standardize case to ensure exact matches
  # pathway_genes <- trimws(toupper(pathway_genes))
  pathway_genes <- trimws(pathway_genes)

  #colnames(expr_matrix) <- toupper(colnames(expr_matrix))
  selected_genes <- intersect(colnames(expr_matrix), pathway_genes)

  #  Extraction of expression subsets of related genes
  expr_submatrix <- expr_matrix[, selected_genes]

  expr_scaled <- scale(expr_submatrix)
  cor_matrix <- cor(expr_scaled, method = "pearson")
  cor_matrix[lower.tri(cor_matrix)] <- t(cor_matrix)[lower.tri(cor_matrix)]
  diag(cor_matrix) <- 1

  return(cor_matrix)
}


main <- function(data_path, data_num=1, view_num=4, pathway_name) {
  for(i in 1:data_num) {
    mat_name = paste('data', '.csv', sep='')
    mat_gene = load_matrix_for_GSE(paste(data_path, mat_name, sep='\\'))
    # mat_gene = t(mat_gene)
    coexp_matrices = {}
    for(j in 1:view_num) {
      save_path_W = paste(data_path, paste(paste('sm_', i, '_', sep=''), j, '.csv',sep=''), sep='\\')
      save_path_mat_path = paste(data_path, paste(paste('mat_path_', i, '_', sep=''), j, '.csv',sep=''), sep='\\')

      gSet = get_mat_path(mat_gene, paPath, pathway_name[j])

      pathway_genes <- geneIds(gSet)
      coexp_matrix <- get_coexpression_matrix(mat_gene, pathway_genes)

      coexp_matrices[j] = coexp_matrix
      print(coexp_matrix)
      write.csv(coexp_matrix, save_path_W)
    }


    mat_sm_gene = get_gene_coexpression_matrix(mat_gene)
    write.csv(mat_sm_gene, paste(data_path, paste('sm_', i, '_', view_num+1, '_.csv', sep=''), sep='\\'))
  }
}


paPath = "C:\\Users\\Administrator\\Desktop\\StudyResource\\scPML\\scPML-main\\demo\\pathway_data\\human"
mouse_pathway_name = c('KEGG', 'Reactome', 'Wikipathways', 'biase')
human_pathway_name = c('KEGG', 'Reactome', 'inoh', 'pathbank','panther','humancyc')
# human_pathway_name = mouse_pathway_name
args = commandArgs(trailingOnly = TRUE)
# base_path = args[[1]]
base_path="F:\\StudyResource\\Project\\data\\AD-brain"

pathway_name = human_pathway_name
# ref
data_path = paste(base_path, sep='\\')
main(data_path, data_num=1, view_num=6, pathway_name)
