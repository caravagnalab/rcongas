

#' A simple function to divide the genome in segments
#'
#' The functions first finds the maximum number of segments of a given length for each chromosome in the
#' human genome (hg38 assembly, onlys chrs 1:22, X and Y) then it sample the number of actual splits
#' in a range (0,MAX). In the end it adds some stochasticity to the segment ends and starts.
#'
#'
#'
#' @param approx_length Approximate segment legnth. Integer
#'
#' @return A data frame with genomic segments (with 3 cols chr, from, to)
#'
#'
#'

divide_genome_in_segments <- function(approx_length = 1e+8){


  genome <- hg38_karyo
  good_chr <- names(genome) %in% c(1:22, "X")
  genome <- genome[good_chr]
  max_segments <- (genome %/% approx_length) + 1
  segments <- vapply(X = max_segments, FUN = sample.int, size=1, FUN.VALUE = integer(1))
  breaks <- mapply(genome, segments + 1, FUN = function(x,y) seq(1,x,length.out = y))
  rand_breaks <- mapply(breaks, names(breaks),FUN = function(x,y) randomize_breaks(x,y, length = approx_length, chr_lengths = genome))
  splitted_breaks_df <- mapply(rand_breaks,names(rand_breaks), FUN = create_df_from_range, SIMPLIFY = F)
  coord_df <- Reduce(splitted_breaks_df, f = rbind2, init = NULL)
  coord_df <-  coord_df %>% dplyr::mutate(from = as.integer(from), to = as.integer(to))
  return(coord_df)

}


#' Create a dataframe from a breakpoint list
#'
#' Convert a chr and a breakpoint vector in a segments data frame
#'
#' @param breaks vector of breakpoints for a given chromosome
#' @param name chhr name
#'
#' @return a data frame with genomic segments (with 3 cols chr, from, to)
#'

create_df_from_range <- function(breaks, name){
  l <- length(breaks)
  df <- data.frame(chr = rep(name, l-1),from = breaks[1:(l-1)], to = breaks[2:l])
  return(df)
}

#' Add stocasticity in the brakpoints
#'
#' Sample the noise from a uniform distribution [-(seg_length * shuffle_perc), (seg_length * shuffle_perc)]
#' and add it to the breakpoint
#'

randomize_breaks <- function(breaks, chr,length, chr_lengths,shuffle_perc = 0.1){
  l <- length(breaks)
  if(l <= 2)
    return(breaks)
  shuffle_int <- round(runif(n = l-2, min = -(length * shuffle_perc), max = (length * shuffle_perc)))
  return(c(breaks[1], ifelse(breaks[2:(l-1)] + shuffle_int  > chr_lengths[chr], chr_lengths[chr], breaks[2:(l-1)] + shuffle_int), breaks[l]))
}

#'Generete the ploidy for a segment
#'
#'@param length integer, length of the sample (basically it should be always one)
#'@param probs vector with the same length as karyo, the probability of picking a given CNV value
#'@param karyo vector of integer, CNV states allowed
#'
#'@return an integer, the CNV state for tthat segment
#'
generate_ploidy <- function(length, probs, karyo){
  return(sample(karyo, length, prob = probs, replace = T))
}

#' Generat the ploidy given a parent clone
#'
#'
#' @param paret_ploidy vector of integers, the ploidy of the parent clone
#' @param spots integer, how many sites to be changed
#' @param random boolean, should the spots be sampled randomly
#' @param non_random_spots vector of integers,
#' @param changes vector of integers, allowed values for deletion and amplification
#'
#' @return vector of integers, ploidy profile of the child subclone
#'
#'

generate_ploidy_parent <- function(parent_ploidy, spots, random, non_random_spots, changes) {
  breaks <- seq_along(parent_ploidy)
  if(random)
    mutation_spots <- sample(breaks, size = spots)
  else
    mutation_spots <- non_random_spots
  parent_ploidy <- change_spot(parent_ploidy, mutation_spots, changes)
  return(parent_ploidy)
}


#' Add a deletion or an insertion in a given segment
#'
#'
#' @param parent_ploidy vector of integers, CNV value in the parent clone
#' @param mutation_spots vector of integers, index of the segmnts to be changed in the parent
#' @param changes vector of integers, allowed values for deletion and amplification
#'
#' @return the vector with the new subclonal ploidy
#'
#'

change_spot <- function(parent_ploidy, mutation_spots, changes){
  for(i in mutation_spots){
    parent_ploidy[i] <- parent_ploidy[i] + sample(changes, size = 1)
  }
  return(parent_ploidy)
}


#'Generat the ploidy given a parent clone (wih a finer control)
#'
#'It is basically equal to {\link{generate_ploidy_parent}} but provides the direct controll of the segments
#'in which we want to place the mutation
#'
#'
#'@param spots_deletion segments where we want to insert a single chromosome deletion
#'@param spots_insertion segments where we want to insert a single or double (randomly sampled)
#'  chromose amplification
#'@param spots_O segments where we wanto to insert a double chromosome deletion
#'@inheritParams generate_ploidy_parent
#'
#'@return vector of integers, ploidy profile of the child subclone
#'
generate_ploidy_parent_fine <- function(parent_ploidy, spots_deletion, spots_insertions, spots_O){
  breaks <- seq_along(parent_ploidy)
  parent_ploidy <- change_spot(parent_ploidy, spots_deletion, c(-1))
  parent_ploidy <- change_spot(parent_ploidy, spots_insertions, c(1,2))
  parent_ploidy <- change_spot(parent_ploidy, spots_O, c(-2))

  return(parent_ploidy)
}


#' Generate a random tree
#'
#'@param K number of subclones in the population
#'@param unrelated_prob probability of the emergence of a separate population
#'
#'
#'@return a data frae with two columns indicating the cluster and its parent cluster
#'@export


generate_random_tree <- function(K, unrelated_prob = 0){
  res <- matrix(nrow = K, ncol = 2)
  res[,1] <- seq_len(K)
  for(k in seq_len(K)){
    if(k==1){
      res[k,2] <- res[k,1]
    } else if(rbernoulli(1,p = unrelated_prob)){
      res[k,2] <- res[k,1]
    } else{
      res[k,2] <- sample(size = 1,x = res[1:(k-1),1])
    }
  }
  res <- as.data.frame(res)
  colnames(res) <- c("cluster", "parent")
  return(res)
}


generate_baf <-  function(cnv_df, ploidy_to_baf = list("0" = 0, "1" = 0, "2" = 0.5, "3" = 0.33, "4" = c(0.5,0.25), "5" = c(0.2,0.4), "6" = c(0.17,0.33,0.5))
                         ){

  pld <- cnv_df %>% dplyr::select(matches("ploidy") & !matches("real"))
  res <- data.frame(matrix(nrow = nrow(pld), ncol = ncol(pld)))
  for(i in seq_len(ncol(res)))
    res[,i] <- sapply(ploidy_to_baf[as.character(pld[,i])], function(x) sample(x = x, size = 1))
  colnames(res) <- paste0("baf", 1:ncol(res))
  return(res)

}

#'
#' @param K integer, number of clusters to simulate
#' @param tree a data frame with two columns, the first with the index of each cluster and the second with its parent
#' @param spots number of segment to change when using random mode
#' @param probs probabilities of picking a given CN value from the karyotype, needs to have length = length(karyo)
#' @param random shold the mutated segments among parent-child subclones be choosen at random?
#' @param changes what are th allowed changes in ploidy between parent-child subclnes when fine = FALSE
#' @param fine allow a finer controll on the different segments
#' @param karyo integer vector, alloed CN values
#' @param type one of c("parent", "random"), indicates if the clusters have an evolutionary relionship or they are unrelated
#' @param non_random_spots if random = FALSE, wihch spotss should harbor a mutation?
#'@param spots_deletion segments where we want to insert a single chromosome deletion
#'@param spots_insertion segments where we want to insert a single or double (randomly sampled)
#'  chromose amplification
#'@param spots_O segments where we wanto to insert a double chromosome deletion
#'@param approx_length approximate segment length
#'
#'@return df with segmentss and ploidy information

generate_cluster_ploidy_df <- function(K= 2, spots = 5,  probs = c(0.0,1,0.0,0.0),div_factor_dist = 1e6,
                                       random = TRUE, changes = c(-1,1,2), fine = FALSE, karyo = 1:4, type="parent", non_random_spots = c(17,18),
                                       spots_deletion = c(17,18), spots_insertions = c(8), spots_O =  c(),
                                       tree = generate_random_tree(K), segments = divide_genome_in_segments(), baf = TRUE){
  length <-  nrow(segments)



  res <- matrix(nrow = length, ncol = K)
  if(type == "unrelated"){
    for(j in seq_len(ncol(res))){
      res[,j] <- generate_ploidy(length, probs, karyo)
    }
  } else if(type == "parent"){
    for(clusters in seq_len(nrow(tree))){

      if(tree[clusters,1] == tree[clusters,2]){
        res[,clusters] <- generate_ploidy(length, probs, karyo)
      } else {
        if(!fine)
          res[,clusters] <- generate_ploidy_parent(res[,tree[clusters,2]], spots = spots ,random = random, changes = changes, non_random_spots = non_random_spots)
        else
          res[,clusters] <- generate_ploidy_parent_fine(res[,tree[clusters,2]], spots_deletion = spots_deletion, spots_insertions = spots_insertions, spots_O = spots_O)
      }
    }

  }

  colnames(res) <- paste0("ploidy", seq_len(K))

  cnv_df <- cbind2(segments, res)
  cnv_df <- cnv_df %>% mutate(dist = to - from) %>% mutate(mu = rnbinom(nrow(cnv_df),mu = dist/div_factor_dist + 3, size = 6))

  if(baf)
    cnv_df <- cbind(cnv_df,generate_baf(cnv_df))

  return(cnv_df)

}



trivial_hap <-  function(cnv_df){

  tot_genes <- sum(cnv_df$mu)
  A <- rep(TRUE, tot_genes)
  B <- rep(FALSE, tot_genes)
  return(data.frame(A,B))
}


#' Run a complete simulation
#'
#'This function runs a n entre simulation, it first segments the genome, then it assigns copy number to each clone based on the given tree
#'
#'
#'
#' @param ncells integer, number of cells to simulate
#' @param props proportions of each cluster, needs to have length = K
#' @param theta_shape shape parameter for the Gamma modelling the library factor
#' @param theta_rate rate parameter for the Gamma modelling the library factor
#' @return a CNVSimulation object
#' @export
#'
#' @examples
#'
run_simulation <- function(cnv_df, ncells = 1000,props = c(0.8,0.2), K = 2,theta_loc = 23000, theta_shape = 0.3,  model_matrix = NULL, coeff = c(0),class_rand = c(2,2), k_rand = 2,random_cov = FALSE, het_snps_for_gene = 1,
                           prob_random = c(0.99,0.1), random_class_dim = c(0.2,0.1),perc_dropout = 0, baf_var = 20, gamma_shape = 0.45, gamma_rate = 0.05, monoallelic_rate = 0.8,
                           haplotypes = trivial_hap(cnv_df), max_snps = 5, fast = FALSE, gamma_shape_fast = 6, gamma_rate_fast = 6, theta_shape_fast = 6,theta_rate_fast = 4){

  if(length(props) != K){
    stop("vector of cellular proportions must have the same length of the number of clusters")
  }

  # generate number of genes and prepare result matrices

  props <-  sort(props, decreasing = T)
  baf <-  cnv_df %>% select(matches("baf"))
  res <- matrix(nrow = ncells, ncol = nrow(cnv_df))
  cnv_mat <- matrix(nrow = ncells, ncol = nrow(cnv_df))

  if(fast){
    theta_vec <- rgamma(ncells,theta_shape_fast, theta_rate_fast)
  }
  else{
    theta_vec <- round(rlnorm(ncells,log(theta_loc), theta_shape))
  }

  clust_prop_idx <- round(cumsum(c(0,ncells * props)))
  tot_genes <- sum(cnv_df$mu)
  gene_counts <- matrix(nrow = ncells, ncol = tot_genes)
  snps <-  rpois(tot_genes, het_snps_for_gene)
  perc_read <- snps / max_snps
  minor_freq <- matrix(0L,nrow = ncells, ncol = tot_genes)
  gene_cumsum <-  cumsum(c(0,cnv_df$mu))
  genes_to_seg <- vector(mode = "list", length = nrow(cnv_df))
  baf_mat <- matrix(nrow = ncells, ncol = tot_genes)

  seg_names <-   paste(cnv_df$chr, cnv_df$from, cnv_df$to, sep = ":")
  cell_names <-  paste0("cell", seq_len(nrow(res)))
  gene_names <-  paste0("gene",seq_len(tot_genes))

  names(genes_to_seg) <- seg_names


  # build the linear model

    if(!is.null(model_matrix)){

      model_eff <- model_matrix %*% coeff

    } else {
      model_eff <- matrix(1, nrow = ncells, ncol = tot_genes)
    }

    if(random_cov){

      rand_matrix <- matrix(0, nrow = ncells, ncol =k_rand)
      for(k in seq_len(k_rand)){
        rand_matrix[,k] <-  ifelse(rbernoulli(ncells, p = random_class_dim[k]),1,0)

      }
      coeff_rand <- matrix(sapply(1:(tot_genes * k_rand),function(x) sample(x = c(0, runif(1)), prob = c(0.99,0.01), size = 1, replace = T)),nrow=k_rand,ncol=tot_genes)
      random_eff <-  rand_matrix %*% coeff_rand

    } else{
      random_eff <-  matrix(1, nrow = ncells, ncol = tot_genes)

    }

  drp <- matrix(rbernoulli(tot_genes * ncells, 1-perc_dropout), nrow = ncells, ncol = tot_genes)
  if(fast){
    gene_means <-  rgamma(tot_genes, gamma_shape_fast, gamma_rate_fast)
  } else{
    gene_means <-  rgamma(tot_genes, gamma_shape, gamma_rate)
  }
  weights <-  drp * gene_means

  # generate gene counts


  for(j in seq_len(ncol(res))){
    genes_curr <- gene_names[csidx(gene_cumsum,j)]
    genes_to_seg[[j]] <- genes_curr
    for(i in seq_len(K)){

      actual_ploidy <- cnv_df %>%  dplyr::select(paste0("ploidy", i)) %>%  unlist()

      genes_idx <-  csidx(gene_cumsum,j)
      cells_idx <-  csidx(clust_prop_idx,i)
      genes_in_seg <-  length(genes_idx)
      cells_in_cluster <-  length(cells_idx)

      if(fast){
        gene_counts[cells_idx, genes_idx] <-
          rpois(genes_in_seg * cells_in_cluster,
                ( actual_ploidy[j] * theta_vec[cells_idx] *
                    random_eff[cells_idx,genes_idx] *
                    model_eff[cells_idx,genes_idx] * weights[cells_idx,genes_idx])
                + 1e-8)


      } else{

        gene_counts[cells_idx, genes_idx] <-
            rpois(genes_in_seg * cells_in_cluster,
                  ( actual_ploidy[j] *
                      random_eff[cells_idx,genes_idx] *
                      model_eff[cells_idx,genes_idx] * weights[cells_idx,genes_idx])
                  + 1e-8)

      }

      minor_freq[cells_idx, genes_idx] <-
        rbeta(genes_in_seg * cells_in_cluster,
              baf[j,i] * 300,
              (1-baf[j,i]) * 300)

      cnv_mat[cells_idx,j] <- rep(x = actual_ploidy[j], cells_in_cluster)
      baf_mat[cells_idx, genes_idx] <- sapply(genes_idx, function(x) {
                                                                  rep(x = baf[j,i], cells_in_cluster)})

    }
  }

  if(!fast){

    gene_counts <-
      extraDistr::rmvhyper(length(theta_vec),gene_counts,theta_vec)

  }




  colnames(gene_counts) <- gene_names
  rownames(gene_counts) <- cell_names


  minor_freq <-  ifelse(minor_freq > 0.5, 1 - minor_freq, minor_freq)
  mask <- rbernoulli(length(minor_freq),monoallelic_rate)

  minor_freq[mask] <- rbinom(length(minor_freq),1, minor_freq)[mask]


  minor_counts <- round(gene_counts * minor_freq)

  major_counts <- round((gene_counts) - minor_counts)


  colnames(major_counts) <- gene_names
  rownames(major_counts) <- rownames(res)

  colnames(minor_counts) <- gene_names
  rownames(minor_counts) <- rownames(res)


  BAF_df <- colSums(minor_counts) / (colSums(minor_counts) + colSums(major_counts))

  BAF_df_seg <-  vector(length = length(genes_to_seg))

  B_allele <-  matrix(0L,nrow = ncells, ncol = tot_genes)
  A_allele <- matrix(0L,nrow = ncells, ncol = tot_genes)

  B_allele[,haplotypes$B] <- major_counts[,haplotypes$B]
  B_allele[,!haplotypes$B] <-  minor_counts[,!haplotypes$B]
  A_allele[,haplotypes$A]  <- major_counts[,haplotypes$A]
  A_allele[,!haplotypes$A] <-  minor_counts[,!haplotypes$A]

  colnames(A_allele) <- gene_names
  rownames(A_allele) <- rownames(res)

  colnames(B_allele) <- gene_names
  rownames(B_allele) <- rownames(res)

  for(j in seq_len(length(genes_to_seg))) {

    genes <- genes_to_seg[[j]]
    res[,j] <- rowSums(gene_counts[, genes])
    BAF_df_seg[j] <- mean(BAF_df[genes])
  }


  colnames(res) <- seg_names
  colnames(cnv_mat) <- seg_names
  rownames(res) <- cell_names
  rownames(cnv_mat) <- cell_names
  cnv_df$ploidy_real <-cnv_df  %>% select(matches("ploidy")) %>% as.matrix(.) %*% props
  cnv_df$baf_real <- cnv_df  %>% select(matches("baf")) %>% as.matrix(.) %*% props


  clust_ids <-  rep(seq_len(K), round(props * ncells))
  clust_ids <- as.data.frame(clust_ids)
  colnames(clust_ids) <- "cluster_id"
  rownames(clust_ids) <- rownames(res)






  return(structure(list(counts = res, cnv_mat = cnv_mat, cnv= cnv_df, clust_ids = clust_ids,
                        gene_counts = gene_counts, allelic_freq = list(
                        B_allele_counts = B_allele,
                        A_allele_counts = A_allele,
                        BAF = BAF_df, BAF_seg = BAF_df_seg),
                        covariates = list(known = model_eff, random = random_eff), params = list(theta_loc = theta_loc, theta_shape = theta_shape,
                                                                                                 gamma_shape = gamma_shape, gamma_rate = gamma_rate
                                                                                                 )
                        , nCount_RNA = theta_vec ,norm_factors = scale(theta_vec) * mean(gene_means), from_genes_to_seg = genes_to_seg),
                        class = "CNVSimulation"))
}




run_simulation_generative <- function(cnv_df, ncells = 1000, props = c(0.8,0.2), K = 2 ,theta_rate = 8, theta_shape = 8, nbinom = FALSE, size = 6, perc_genes = 1){




  res <- matrix(nrow = ncells, ncol = nrow(cnv_df))
  cnv_mat <- matrix(nrow = ncells, ncol = nrow(cnv_df))
  theta_vec <- rgamma(nrow(res),theta_shape, theta_rate)



  clust_prop_idx <- round(cumsum(c(0,ncells * props)))

  for(j in seq_len(ncol(res))) {

    noise <- rnorm(1,0, 0.1)
    for(i in seq_len(K)){

      cells_idx <-  csidx(clust_prop_idx,i)
      actual_ploidy <- cnv_df %>%  select(paste0("ploidy", i)) %>%  unlist()
      if(i > 1) {
        parent_ploidy <- cnv_df %>%  select(paste0("ploidy", i-1)) %>%  unlist()
      } else {
        parent_ploidy <-  actual_ploidy
      }
      actual_sum <-  sum(actual_ploidy *cnv_df$mu) / sum(cnv_df$mu)
      parent_sum <- sum(parent_ploidy * cnv_df$mu) / sum(cnv_df$mu)
      mu_k <- ifelse(cnv_df$mu[j] > 0, cnv_df$mu[j], 1)
      if(nbinom){

        subclone_counts <- rnbinom(round(ncells * props[i]),mu =  (mu_k * perc_genes * theta_vec * (actual_ploidy[j] + noise) + 1e-8) / actual_sum, size = size)
        clone_counts <-  rnbinom(round(ncells * props[i]),mu =  (mu_k * (1-perc_genes) * theta_vec * (parent_ploidy[j] + noise) + 1e-8) / parent_sum, size = size)

      } else {
        subclone_counts <- rpois(round(ncells * props[i] ), (mu_k * perc_genes * theta_vec * (actual_ploidy[j] + noise) + 1e-8) / actual_sum)
        clone_counts <- rpois(round(ncells * props[i] ), (mu_k * (1-perc_genes) * theta_vec * (parent_ploidy[j] + noise) + 1e-8) / parent_sum)


      }
      cnv_mat[cells_idx,j] <- rep(actual_ploidy[j], length(cells_idx))
      res[cells_idx,j] <- subclone_counts + clone_counts
      subclone_id <- rep(x = actual_ploidy[j], round(ncells * props[i]))
    }

  }

  clust_ids <-  rep(seq_len(K), round(props * ncells))
  clust_ids <- as.data.frame(clust_ids)



  colnames(res) <- paste(cnv_df$chr, cnv_df$from, cnv_df$to, sep = ":")
  colnames(cnv_mat) <- colnames(res)
  rownames(res) <- paste0("cell", seq_len(nrow(res)))
  rownames(cnv_mat) <- rownames(res)
  colnames(clust_ids) <- "cluster_id"
  rownames(clust_ids) <- rownames(res)
  cnv_df$ploidy_real <-cnv_df  %>% select(matches("ploidy")) %>%   as.matrix(.) %*% props
  data <-  list(counts = res, cnv_mat = cnv_mat, cnv= cnv_df, clust_ids = clust_ids,
                          params = list(theta_rate = theta_rate, theta_shape = theta_shape, theta = theta_vec
                          )
  )
  x <-structure(list(data = data, reference_genome = "hg38"),
  class = "CNVSimulation")

  return(x)
}



write_results <- function(df, out_prefix, sep = ",") UseMethod(write_results)

write_results.CNVSimulation <- function(df, new_dir = FALSE,dir_name,  out_prefix, sep = ",") {

  if(new_dir){
    dir.create(file.path(".", new_dir), showWarnings = FALSE)
    out_prefix <- file.path(".", new_dir, out_prefix)
  } else {
    out_prefix <-  file.path(".", out_prefix)
  }

  write.table(x = df$cnv_df , file = paste0(out_prefix,"_cnv.csv"), sep = sep, row.names = FALSE)
  write.table(x = df$counts , file = paste0(out_prefix,"_data.csv"), sep = sep)
  write.table(x = df$cnv_mat , file = paste0(out_prefix,"_cnv_mat.csv"), sep = sep)
  write.table(x = df$clust_ids , file = paste0(out_prefix,"_clust_ids.csv"), sep = sep)
  write.table(x = df$gene_counts , file = paste0(out_prefix,"_gene_counts.csv"), sep = sep)

}


`[.CNVSimulation` <- function(x, i, j) {
    x$data$cnv <- x$cnv[j,]
    x$data$counts <- x$counts[i,j]
    x$data$clust_ids <-  x$clust_ids[i,, drop = FALSE]
    x$data$cnv_mat <- x$cnv_mat[i,j]

    return(x)
}



