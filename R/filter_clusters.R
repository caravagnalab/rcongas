# #' Filter output clusters
# #'
# #' @description From a CONGAS+ fit object remove clusters that do not pass filters.
# #' These include number of cells and abundance of cluster (mixing proportions). This function the assigns the 
# #' cells back to the most similar remaining cluster. In case the inference was performed using MAP we choose the cluster with
# #' the most similar CNV profile (based on euclidean distance), otherwise we assign the cell to the cluster with second highest 
# #' probability and renormalize the z.
# #'
# #'
# #' @param x Rcongas object
# #' @param ncells minimum size of a valid cluster expressed as absolute number of cells
# #' @param abundance minimum size of a cluster expressed as a percentage of the total cell number 
# #'
# #' @return
# #' @export
# #'
# #' @examples
# #'
# #' x = Rcongas::example_object
# #'
# #' print(x)
# #'
# #' # Equivalent filters for this model
# #' x %>% filter_clusters(abundance = .5)
# #'
# #' x %>% filter_clusters(ncells = 150)
# #'
# #' 


# filter_clusters <- function() {
  
# }