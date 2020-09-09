filter_bins <- function(conf, nmin_genes = 1, dim = 1000000){

  inp$cnv

}


binned_mean <- function(x, annot, binSize = 100, FUN = mean){

  if(binSize > nrow(x)){
    warning(paste0("Bin too large using ", nrow(x), " instead."))
    binSize <- nrow(x)
  }

  rest <- nrow(x)%%binSize
  res <- matrix(NA, nrow = floor(nrow(x)/binSize) + ifelse(rest == 0, 0, 1), ncol = ncol(x))
  rnames <- vector()
  k <- 1
  for(i in seq(1,nrow(x),binSize)){
    for(j in seq_len(ncol(x))){
      if(i + binSize - 1 <= nrow(x)){
        res[k,j] <- FUN(x[i:(i+binSize-1),j])

      }else{
        res[k,j] <- FUN(x[(i-(binSize-rest)):(i+rest-1),j])
      }
    }
    if(i + binSize - 1 <= nrow(x))
      rnames[k] <- paste0(annot$chromosome_name[i],":",annot$start_position[i], ":", annot$end_position[i + binSize-1])
    else
      rnames[k] <- paste0(annot$chromosome_name[i],":",annot$start_position[i], ":", annot$end_position[i + rest-1])

    k <- k + 1
  }
  colnames(res) <- colnames(x)
  rownames(res) <- rnames
  return(res)
}

binned_recalibration <- function(x, annot, binSize = 100, FUN = mean){

  if(binSize > ncol(x)){
    warning(paste0("Bin too large using ", ncol(x), " instead."))
    binSize <- ncol(x)
  }

  rest <- nrow(x)%%binSize
  res <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  rnames <- vector(length = nrow(x))
  for(i in seq(binSize/2, nrow(x) - binSize/2)){

    res[i,] <- apply(x[(i - binSize/2):(i+binSize/2-1),],2,FUN)
    rnames[i] <- paste0(annot$chromosome_name[i],":",annot$start_position[i], ":", annot$end_position[i])
  }
  colnames(res) <- colnames(x)
  rownames(res) <- rnames
  return(res)
}


correct_bins <- function(bins, filt = 2.0e+07, hard = FALSE){

  if(hard){
    bins <- bins %>% filter(bins$dist > filt)
  }
  res <- data.frame()
  i <- 1
  while(i < nrow(bins)){

    tmp <- bins[i,]
    k <- 0
    while(bins$tot[i] == bins$tot[i+1] & bins$chr[i] == bins$chr[i+1] & i < nrow(bins)){
      tmp <- rbind2(tmp, bins[i+1,])
      i <- i + 1
      k <- k + 1
    }
    new <- data.frame(bins$chr[i], bins$from[i-k], bins$to[i], bins$tot[i])
    res <- rbind2(res,new)
    i <- i + 1
  }

  colnames(res) <- c("chr","from", "to","tot")
  return(res)
}

custom_binned_apply <- function(x, annot, bins,FUN = mean, binDims = 100){

  res <- matrix(NA ,ncol = ncol(x))
  rnames <- vector()
  genesBin <- vector()
  k <- 1

  for(i in seq_len(nrow(bins))){
    mask <- ifelse(annot$start_position > bins$from[i] & annot$end_position < bins$to[i], TRUE, FALSE)
    genBin <- sum(mask)
    if(genBin == 0){
      res <- rbind(res, rep(0, ncol(x)))
      rnames[k] <- paste0(bins$chr[i],":",bins$from[i], ":", bins$to[i])
      genesBin[k] <- 0
      k <- k + 1
      next
    }

    if(genBin  < binDims)
      NbinDims <- genBin
    else
      NbinDims <- binDims


    for(bin in seq(1,genBin,NbinDims)){

      end <- bin + (NbinDims - 1)
      if(end > genBin ) next
      toBind <- apply(x[which(mask)[bin:(bin + (NbinDims - 1))],],2 ,FUN , na.rm = T)
      res <- rbind(res,toBind)

      if(bin == 1 & NbinDims == binDims)
        rnames[k] <- rnames[k] <- paste0(annot$chromosome_name[which(mask)[bin]],":",bins$from[i], ":", annot$end_position[which(mask)[end]])
      else if(bin == 1 & NbinDims != binDims)
        rnames[k] <- rnames[k] <- paste0(annot$chromosome_name[which(mask)[bin]],":",bins$from[i], ":", bins$to[i])
      else
        rnames[k] <- rnames[k] <- paste0(annot$chromosome_name[which(mask)[bin]],":",annot$start_position[which(mask)[bin]], ":", annot$end_position[which(mask)[end]])

      genesBin[k] <- NbinDims
      k <- k + 1

    }

    rest <- genBin %% NbinDims

    if(rest != 0){
      toBind <- apply(x[which(mask)[(genBin - rest + 1):genBin],], 2, FUN, na.rm = T)
      res <- rbind(res,toBind)
      rnames[k] <- paste0(annot$chromosome_name[which(mask)[(genBin - rest + 1)]],":",annot$start_position[which(mask)[(genBin - rest + 1)]], ":", bins$to[i])
      genesBin[k] <- rest
      k <- k + 1
    }
  }

  res <- res[-1,]
  colnames(res) <- colnames(x)
  rownames(res) <- rnames
  return(list(matrix = res, binDims = genesBin))

}


custom_fixed_binned_apply <- function(x, annot, bins,FUN = mean){

  if(is.null(x)) return(list(matrix = NULL, binDims = NULL))
  res <- matrix(NA ,ncol = ncol(x), nrow = nrow(bins))
  rnames <- vector(length = nrow(bins))
  expressed_genes <- matrix(NA ,ncol = ncol(x), nrow = nrow(bins))
  genesBins <- vector(length = nrow(bins))

  for(i in seq_len(nrow(bins))){
    mask <- annot$start_position > bins$from[i] & annot$end_position < bins$to[i]
    genBin <- sum(mask)
    rnames[i] <- paste0(bins$chr[i],":",bins$from[i], ":", bins$to[i])
    genesBins[i] <- genBin

    if(genBin == 0) {
      res[i,] <- rep(0, ncol(x))
    }else{
      expressed_genes[i,] <- apply(x[which(mask),],2 ,function(x) return(sum(x != 0, na.rm = T)))
      res[i,] <- apply(x[which(mask),],2 ,FUN , na.rm = T)
    }
  }

  colnames(res) <- colnames(x)
  rownames(res) <- rnames
  colnames(expressed_genes) <- colnames(res)
  rownames(expressed_genes) <- rownames(res)
  return(list(matrix = res, binDims = expressed_genes, binDims_fixed = genesBins))
}


getChromosomeDF <- function(df, online = FALSE, genome = "hg38", chrs = c(1:22,"X")){

  if(online & !require("biomaRt")) stop("Please install biomaRt if you want to use the online functionality")

  gene_names <- rownames(df)


  if(online){

    ensembl <- biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
    genes_position <- biomaRt::getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                                           'start_position', 'end_position'), filters = 'hgnc_symbol', values = gene_names, mart = ensembl)
  }else {
    genes_position <- eval(parse(text = paste0(genome, "_genes")))
    genes_position <-  genes_position %>% filter(hgnc_symbol %in% gene_names)
  }

  genes_position <- genes_position %>%  filter(chromosome_name %in% chrs)
  genes_position <- genes_position[with(genes_position,gtools::mixedorder(paste(chromosome_name, start_position, sep =  ":"))),]
  df <- df[genes_position$hgnc_symbol,]
  res <- split(as.data.frame(df), as.factor(genes_position$chromosome_name))
  res <- res[mixedorder(names(res))]
  genes_position <- split(genes_position, genes_position$chromosome_name)
  genes_position <- genes_position[mixedorder(names(genes_position))]

  return(list(res = res, pos = genes_position))

}




filterSC <- function(x, geneCut = 0.05, cellCut = 3000, capping = FALSE, cap_quantile = 0.95){

  filt <- ncol(x) * geneCut
  x <-  x[rowSums(x) > filt,]

  x <-  x[, which(colSums(x != 0) > cellCut)]

  if(capping)
    x <-  cap_genes(x, cap_quantile)

  return(x)

}

cap_genes <- function(x, quantile = 0.95){

  x <- H13001_heat[rowMeans(x) < quantile(rowMeans(x), probs = 0.95),]
  return(x)
}




get_data <- function(data, bindim = 100,chrs = c(1:22,"X"), filter = NULL, fun = sum,
                     type = c("binning", "smoothing","fixed_binning"), cnv_data = NULL, median_bindims = TRUE,
                     online = FALSE, genome = "hg38", startsWithchr = FALSE, correct_bins = TRUE)

{

  #data <- filterSC(data)


  if(!is.null(filter)){
    data <- data[rownames(data) %in% filter, ]
  }

  df_list <- getChromosomeDF(data, genome = genome, online = online, chrs = chrs)

  data_splitted <- df_list[[1]]
  chrs_cord <- df_list[[2]]



  if(type == "binning"){
    data_binned <- mapply(data_splitted, chrs_cord ,FUN = function(x,y) binned_mean(x,y,binSize=bindim ,fun))
    result <- t(do.call(rbind, data_binned))
    chrs_cord <- t(do.call(rbind, chrs_cord))
    return(list(counts = result,genes_cords  = chrs_cord))
  } else if(type == "smoothing"){
    data_binned <- mapply(data_splitted, chrs_cord ,FUN = function(x,y) binned_recalibration(x,y,binSize=bindim ,fun))
    result <- t(co.call(rbind, data_binned))
    chrs_cord <- t(do.call(rbind, chrs_cord))
    return(list(counts = result,genes_cords  = chrs_cord))
  } else if(type == "fixed_binning"){

# cnv_data <- CNAqc::smooth_segments(cnv_data, maximum_distance = 10^8)

    if(startsWithchr){
      cnv_data <- cnv_data %>% mutate(chr = gsub(chr, pattern = "chr", replacement = ""))
    }

    if(all(c("start", "end") %in% colnames(cnv_data))){
      cnv_data <- cnv_data %>% mutate(from = start, to = end)
    }

    cp <- cnv_data %>% mutate(dist = to - from,  to = as.integer(to), from = as.integer(from)) %>% filter(chr %in% chrs)

    if(correct_bins)
      cp <- correct_bins(cp, filt = 1.0e+07, hard = TRUE)
    bins_splitted <- split(cp , cp$chr)
    bins_splitted <- bins_splitted[mixedsort(names(bins_splitted))]
    #bins_splitted <- bins_splitted[names(bins_splitted) %in% names(chrs_cord)]

    chrs_cord <- chrs_cord[names(chrs_cord) %in% names(bins_splitted)]
    data_splitted <- data_splitted[names(data_splitted) %in% names(bins_splitted)]


    mask <- mapply(bins_splitted, chrs_cord, FUN = function(x,y){
      res <- vector(length = nrow(x))
      for(i in seq_len(nrow(x))){
        res[i] <- any(x$from[i] <  y$start_position & x$to[i] > y$end_position)
      }
      return(res)
    })


    bins_splitted <- mapply(bins_splitted, mask, FUN =  function(x,y) x[y,], SIMPLIFY = F)

    data_binned <-  mapply(data_splitted, chrs_cord, bins_splitted ,FUN = function(x,y,z) custom_fixed_binned_apply(x,y,z,fun))

    result <- data_binned[seq(1,length(data_binned),3)]
    result_bindim <- data_binned[seq(2,length(data_binned),3)]
    fixed_bindim <- data_binned[seq(3,length(data_binned),3)]

    result <- t(do.call(rbind, result))
    result_bindim <- t(do.call(rbind,result_bindim))
    fixed_bindim <- do.call(c, fixed_bindim )
    colnames(result_bindim) <- colnames(result)
    cp <- do.call(rbind, bins_splitted)

    if(median_bindims)
      cp$mu <- apply(result_bindim,2,median)
    else
      cp$mu <- fixed_bindim
    colnames(cp)[4] <- "ploidy_real"

    return(structure(list(counts = as.matrix(result), bindims = result_bindim, cnv = cp, gene_counts = as.matrix(do.call(rbind, data_splitted))),class = "rcongas" ))
  }
}





`[.rcongas` <- function(x, i, j) {
  x$cnv <- x$cnv[j,]
  x$counts <- x$counts[i,j]
  x$bindims <-  x$bindims[i,, drop = FALSE]

  return(x)
}