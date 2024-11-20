# Jul/12/2019
#
# Script to compute Beta-diversity (Bray/curtis) from counts table.
#
# Supplementary functions:
#     - compute tsne
#     - rotate tsne.
#
#
# Quick usage:
# It requires samplesid, taxonomy and count to initialize it (melted version of OTU table. Recommend taxonomy at genus level.)
# cbd = compute_beta_diversity_and_tsne(samplesid, taxonomy, count)
# cbd$compute_beta_diversity();
# cbd$compute_tsne()
# data_set_tsne = cbd$get_tsne();

#' This function helps you to compute beta-diversity, tsne and PCA plots
#'
#' @param sampleid vector of sampleid. sampleid, taxonomy and cout must have same length
#' @param taxonomy taxonomic classification. sampleid, taxonomy and cout must have same length
#' @param count count. sampleid, taxonomy and cout must have same length
#' @export
#' @name compute_beta_diversity_and_tsne
#' @example
#' We show how to compute a tsne map with 100 samples from
#' frozen_set_ag based on the beta-diversity at genus level
#' get_table_from_database("frozen_set_ag"); #getting samples from frozen set asv_counts_subset = get_counts_subset( frozen_set_ag$sampleid[1:100])
#' get_table_from_database("asv_annotation_blast_ag"); #getting annotation
#' m = merge(asv_counts_subset[,.(asv_key,sampleid,count)],asv_annotation_blast_ag[,.(asv_key,genus)])
#' cbd = compute_beta_diversity_and_tsne(m$sampleid, m$genus, m$count)
#' cbd$compute_beta_diversity();
#' cbd$compute_tsne()
#' data_set_tsne = cbd$get_tsne();
#'
# library("Rtsne");
# library("labdsv");
# library(ape); #For biplot function
compute_beta_diversity_and_tsne <- function(sampleid, taxonomy, count) {
  if (!requireNamespace("labdsv", quietly = TRUE)) {
    warning("The labdsv package must be installed to use this functionality")
    return(NULL)
  }
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    warning("The Rtsne package must be installed to use this functionality")
    return(NULL)
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    warning("The ape package must be installed to use this functionality")
    return(NULL)
  }
  # It will group counts by `taxonomy`

  data_root <- list() # list that is going to load the structure of this function.
  thisEnv <- environment()
  data_root$thisEnv <- thisEnv

  # Data structures that should be returned in the script.
  dt_count <- NULL
  composition_matrix <- NULL
  d_beta <- NULL
  tsne <- NULL
  pca <- NULL
  pca_raw <- NULL
  pcoa <- NULL
  pcoa_raw <- NULL

  data_root$get_composition <- function() {
    return(composition_matrix)
  }
  data_root$get_betadiversity <- function() {
    # 1
    return(d_beta)
  }
  data_root$get_tsne <- function() {
    return(tsne)
  }
  data_root$get_pca <- function() {
    stop("This function is depreciated; please consider PCoA of the composition matrix instead.")
    return(pca)
  }
  data_root$get_pcoa <- function() {
    return(pcoa)
  }

  data_root$get_pca_raw <- function() {
    stop("This function is depreciated; please consider PCoA of the composition matrix instead.")
    return(pca_raw)
  }


  dt_count <- data.table(
    sampleid = sampleid,
    taxonomy = taxonomy,
    count = count
  )

  dt_count[, .(count = sum(count)), by = .(sampleid, taxonomy)]
  dt_count[, count_total := sum(count), by = .(sampleid)]
  dt_count$count_relative <- dt_count$count / dt_count$count_total

  data_root$compute_composition_matrix <- function() {
    if (!is.null(composition_matrix)) {
      print("composition matrix already computed")
      return()
    }
    t1 <- Sys.time()
    composition_matrix <- dcast.data.table(
      data = dt_count,
      formula = sampleid ~ taxonomy,
      value.var = "count_relative",
      fun.aggregate = sum,
      fill = 0
    )
    composition_matrix <- data.frame(composition_matrix)
    row_names <- composition_matrix$sampleid

    # row.names(composition_matrix) = composition_matrix$sampleid;
    composition_matrix <- composition_matrix[, 2:dim(composition_matrix)[2]]
    rownames(composition_matrix) <- row_names
    t2 <- Sys.time()
    cat("Time:Composition_matrix:\n")
    print(t2 - t1)
    # print(row.names(composition_matrix));
    assign("composition_matrix", composition_matrix, data_root$thisEnv)
  }

  data_root$compute_beta_diversity <- function() {
    if (!is.null(d_beta)) {
      print("d_beta matrix already computed")
      return()
    }
    data_root$compute_composition_matrix()

    t1 <- Sys.time()
    method <- "bray/curtis"
    d_beta <- labdsv::dsvdis(composition_matrix, method)
    d_beta <- as.matrix(d_beta)
    t2 <- Sys.time()
    cat("Time:Bray-Curtis matrix:\n")
    print(t2 - t1)

    assign("d_beta", d_beta, data_root$thisEnv)
  }

  data_root$compute_pca <- function() {
    stop("This function is depreciated; please consider PCoA of the composition matrix instead.")
    if (!is.null(pca)) {
      print("pca is already computed")
      return()
    }

    t1 <- Sys.time()
    pca <- stats::prcomp(d_beta)
    t2 <- Sys.time()
    cat("Time:pca\n")
    print(t2 - t1)

    pca_raw <- pca
    pca <- data.table(
      sampleid = rownames(pca$x),
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      PC3 = pca$x[, 3],
      var1 = (pca$sdev[1]^2),
      var2 = (pca$sdev[2]^2),
      var3 = (pca$sdev[3]^2),
      var_total = sum(pca$sdev^2)
    )

    assign("pca", pca, data_root$thisEnv)
    assign("pca_raw", pca_raw, data_root$thisEnv)
  }

  data_root$compute_pcoa <- function() {
    if (!is.null(pcoa)) {
      print("pcoa is already computed")
      return()
    }

    t1 <- Sys.time()
    pcoa_raw <- ape::pcoa(d_beta, rn = rownames(d_beta))
    t2 <- Sys.time()
    cat("Time:pcoa\n")
    print(t2 - t1)

    pcoa <- pcoa_raw
    pcoa <- data.table(
      sampleid = rownames(pcoa_raw$vectors),
      PCoA1 = pcoa_raw$vectors[, 1],
      PCoA2 = pcoa_raw$vectors[, 2],
      PCoA3 = pcoa_raw$vectors[, 3],
      eigen1 = (pcoa_raw$values[1, 1]),
      eigen2 = (pcoa_raw$values[1, 2]),
      eigen3 = (pcoa_raw$values[1, 3]),
      eigen_sum = sum(pcoa_raw$values[1, ]),
      eigen_squared_sum = sum(pcoa_raw$values[1, ]^2)
    )

    assign("pcoa", pcoa, data_root$thisEnv)
    assign("pcoa_raw", pcoa_raw, data_root$thisEnv)
  }


  data_root$compute_tsne <- function(perplexity = 20, max_iter = 3000, theta = 0.1, seed = 1, re_compute = F) {
    if (!is.null(tsne) & !re_compute) {
      print("tsne is already computed")
      return()
    }
    set.seed(seed)

    t1 <- Sys.time()
    t <- Rtsne::Rtsne(d_beta,
      perplexity = perplexity,
      xpca = F,
      check_duplicates = F,
      theta = theta,
      max_iter = max_iter
    )

    t$sampleid <- row.names(d_beta)
    t2 <- Sys.time()
    cat("Time:tsne\n")
    print(t2 - t1)

    tsne <- data.frame(
      sampleid = t$sampleid,
      t1 = t$Y[, 1],
      t2 = t$Y[, 2],
      t1_scaled = f_scale_interval(t$Y[, 1]),
      t2_scaled = f_scale_interval(t$Y[, 2]),
      perplexity = perplexity,
      max_iter = max_iter,
      theta = theta
    )

    tsne$angle <- atan2(tsne$t2_scaled, tsne$t1_scaled)
    tsne$angle[tsne$angle < 0] <- tsne$angle[tsne$angle < 0] + 2 * pi

    assign("tsne", tsne, data_root$thisEnv)
  }

  data_root$compute_tsne_rotate <- function(rotation_angle_rad = 0, reflection_angle_rad = 0, perform_reflection = F) {
    # rotation_angle = (rotation_angle/360) * 2*pi;
    rotation_angle <- rotation_angle_rad
    rotation_matrix <- cbind(
      c(cos(rotation_angle), -sin(rotation_angle)),
      c(sin(rotation_angle), cos(rotation_angle))
    )
    # reflection_angle = (reflection_angle/360)* 2*pi; #I am going to rotate it to proper angle and just reflect it on the main axis.
    reflection_angle <- reflection_angle_rad
    reflection_matrix <- cbind(
      c(cos(2 * reflection_angle), sin(2 * reflection_angle)),
      c(sin(2 * reflection_angle), -cos(2 * reflection_angle))
    ) # Remove or add if needed.
    # reflection_matrix = cbind(c(1,0),c(0,1));
    # td$Y_rotated = td$Y %*%rotation_matrix %*%reflection_matrix; #I rotated the data in 120ª for visualization purpose.
    Y_raw <- as.matrix(tsne[, .(
      t1_scaled,
      t2_scaled
    )])

    Y1 <- Y_raw %*% rotation_matrix # I rotated the data in 120ª for visualization purpose.

    if (perform_reflection) {
      Y2 <- Y1 %*% reflection_matrix
    } else {
      Y2 <- Y1
    }

    tsne$t1_scaled <- f_scale_interval(Y2[, 1])
    tsne$t2_scaled <- f_scale_interval(Y2[, 2])

    tsne$angle <- atan2(tsne$t2_scaled, tsne$t1_scaled)
    tsne$angle[tsne$angle < 0] <- tsne$angle[tsne$angle < 0] + 2 * pi

    assign("tsne", tsne, data_root$thisEnv)
  }

  assign("this", "data_root", data_root$thisEnv)
  class(data_root) <- append(class(data_root), "beta_diversity_set")
  return(data_root)
}


f_scale_interval <- function(x, interval = c(-1, 1)) {
  x_min <- min(x)
  x_delta_range <- diff(range(x))
  x_scaled <- ((x - x_min) / x_delta_range)

  x_scaled <- interval[1] + x_scaled * (diff(interval))
  return(x_scaled)
}



#' This function helps you to compute beta-diversity, tsne and PCA plots
#'
#' @param dists distance matrix from the phyloseq/vegan::distance function
#' @export
#' @name make_tsne
#' @example

make_tsne <-  function(dists, perplexity = 20, max_iter = 3000, theta = 0.1, seed = 1){
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    warning("The Rtsne package must be installed to use this functionality")
    return(NULL)
  }
  tsne_pre <- Rtsne::Rtsne(dists,
                       perplexity = perplexity,
                       xpca = FALSE,
                       check_duplicates = FALSE,
                       theta = theta,
                       max_iter = max_iter
  ) 
  
  #tsne$sampleid <- labels(dists)
  tsne <- data.frame(
    sampleid = labels(dists),
    t1 = tsne_pre$Y[, 1],
    t2 = tsne_pre$Y[, 2],
    t1_scaled = f_scale_interval(tsne_pre$Y[, 1]),
    t2_scaled = f_scale_interval(tsne_pre$Y[, 2]),
    perplexity = perplexity,
    max_iter = max_iter,
    theta = theta
  )
  
  tsne$angle <- atan2(tsne$t2_scaled, tsne$t1_scaled)
  tsne$angle[tsne$angle < 0] <- tsne$angle[tsne$angle < 0] + 2 * pi
  tsne
}
  
