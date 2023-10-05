#' Build Clone Tree form a CONGAS+ object Using medicc2 
#'
#' This function builds a clone tree using the medicc2 tool. The function takes care of input file preparation from CONGAS* results, 
#' running medicc2, cleaning up input files after execution, and plotting the tree using the `ggtree` package.
#'
#' @param x A CONGAS+ object with a valid fit.
#' @param medicc_path The file path of the `medicc2` executable. If not supplied, the function will attempt to find it in the system path.
#' @param result_dir The directory to save the results of medicc2. Defaults to "./medicc_results".
#' @param input_dir The directory where input files for medicc2 will be saved. Defaults to "./medicc_input".
#' @param clean_inputs Logical, should the function delete the input files after the execution of medicc2? Defaults to TRUE.
#' @param plot Logical, should the function plot the resulting tree? Defaults to TRUE.
#' 
#' @return A list containing the newick tree (`tree`) and the plot (`plot`) if `plot` is set to TRUE.
#' 
#' @examples
#' \dontrun{
#' build_clone_tree(x, medicc_path = "/path/to/medicc2")
#' }
#' 
#'
#' @export

build_clone_tree <- function(x, medicc_path = NULL, result_dir = "./medicc_results", input_dir = "./medicc_input", clean_inputs = TRUE, plot = TRUE) {
  
  if(is.null(medicc_path)) {
    if(!is_bin_on_path("medicc2")){
      cli::cli_alert_danger("Please provide a valid path for medicc2 executable!")
      stop()
    }
    medicc_path  <- system("which medicc2", intern = T)
  }
  
  
  input_path <- create_input_for_medicc(x, input_dir = input_dir)
  
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
  
  system(glue::glue("{medicc_path} {input_path} {result_dir} --total-copy-numbers --input-allele-columns value"))
  
  if(clean_inputs) system(glue::glue("rm -rf {input_dir}"))
  
  nw_tree <- ggtree::read.tree(file.path(result_dir, "medicc_input_congas_final_tree.new"))
  nw_plot <- NULL
  if(plot){
    
    nw_plot <- plot_tree_congas(nw_tree)
  }
  
  return(list(tree = nw_tree, plot = nw_plot))
    
  
}


is_bin_on_path = function(bin) {
  exit_code = suppressWarnings(system2("command", args = c("-v", bin), stdout = FALSE))
  return(exit_code == 0)
}


create_input_for_medicc <- function(x, input_dir = "./medicc_input") {
  
  dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
  
  input_table <- get_CNA(x) %>% deidify() %>% 
    rename(sample_id = cluster, chrom = chr, start = from, end = to)
  
  input_path <- file.path(input_dir, "medicc_input_congas.tsv")
  
  write.table(input_path,x =  input_table, sep = "\t", row.names = FALSE)
  
  return(input_path)
}


plot_tree_congas <- function(nw_tree) {
  
  ggtree::ggtree(nw_tree) + ggtree::geom_tiplab() + ggtree::geom_nodepoint(color = "steelblue") + 
    ggtree::geom_tippoint() + ggtree::theme_tree2() + ggtree::geom_rootedge()
  
}
