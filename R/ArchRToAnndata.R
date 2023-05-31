## This script used to covert ArchR object into anndata object
## 2023.5.10
require(ArchR)
require(jsonlite)
require(anndata)
require(Matrix)
require(parallel)
require(tidyverse)
require(SummarizedExperiment)



# Function used to build adata according to the assay type,when multi assay provided
#' @title Build_init_adata
#' @name Build_init_adata
#' @param assay.types vector of the type of assays, mainly peak, gene, motif
#' @param assays vector of assay names to export
#' @param idx.used name the assay
#' @param object ArchR object
#' @param groups vector of cell metadata colnames to export
#' @importFrom ArchR getCellColData getMatrixFromProject
#' @importFrom SummarizedExperiment rowRanges rowData colData assay seqnames start end
#' @importFrom anndata AnnData
#' @return an adata object
#'
Build_init_adata <- function(assay.types, assays, idx.used, object, groups) {
  ## fetech metadata
  if (!is.null(groups) && groups %in% colnames(getCellColData(object))) {
    cell_meta <- as.data.frame(getCellColData(ArchRProj = object))
    cell_meta <- cell_meta[, which(colnames(cell_meta) %in% groups)]
  } else {
    cell_meta <- as.data.frame(getCellColData(ArchRProj = object))
  }

  ###### fetech the exp data
  exp_sce <- getMatrixFromProject(ArchRProj = object, useMatrix = assays[idx.used])
  assay_names <- names(assays(exp_sce))
  # if (length(assay_names) > 1) {
  #   exp_counts <- assay(exp_sce, assay_names[1])
  # } else {
  #   exp_counts <- assay(exp_sce, assay_names)
  # }

  ### gain the feature names
  ## judge the sce object type, mainly:"RangedSummarizedExperiment" and "SummarizedExperiment"
  exp_sce_type <- class(exp_sce)[1]

  if (exp_sce_type == "RangedSummarizedExperiment") {
    rowdata <- rowRanges(exp_sce)
    feature_names <- paste(seqnames(rowdata), start(rowdata), end(rowdata), sep = "-")
    feature_meta <- as.data.frame(rowRanges(exp_sce))
    feature_meta$`feature` <- feature_names
    feature_meta <- feature_meta[, c("feature", "seqnames", "start", "end", "width", "strand", "idx")]
    rownames(feature_meta) <- feature_names
    # feature_names <- paste(seqnames(feature_meta), start(feature_meta), end(feature_meta), sep = "_")
  } else if (exp_sce_type == "SummarizedExperiment") {
    feature_names <- SummarizedExperiment::rowData(exp_sce)$`name`
    feature_meta <- data.frame(feature = feature_meta)
  } else {
    stop("the matrix sce object is not valided")
  }
  rownames(feature_meta) <- feature_meta$`feature`
  ## build the anndata object
  adata <- AnnData(
    X = t(assay(exp_sce, assays[idx.used])),
    obs = cell_meta,
    var = feature_meta
  )

  return(adata)
}




# Function used to add other exp and related data into adata uns
#' @title Add_other_exp
#' @name Add_other_exp
#' @param adata adata object
#' @param idx.used name the assay
#' @param assays vector of assay names to export
#' @param assay.types vector of the type of assays, mainly peak, gene, motif
#' @param object ArchR object
#' @param export_pwm wether to export motif pwm
#' @importFrom ArchR getMatrixFromProject getPeakAnnotation
#' @importFrom SummarizedExperiment rowRanges rowData colData assay
#' @importFrom jsonlite toJSON
#'
Add_other_exp <- function(adata, idx.used, assays, assay.types, object, export_pwm) {

  ## delete the used assays
  assays2 <- assays[-idx.used]
  assay.types2 <- assay.types[-idx.used]

  ## add other exp data into adata
  for (i in seq(length(assay.types2))) {
    exp_sce <- getMatrixFromProject(ArchRProj = object, useMatrix = assays2[i])
    cell_names <- rownames(colData(exp_sce))
    feature_names <- as.vector(rowData(exp_sce)$`name`)
    adata$uns[["assays"]][[assays2[i]]][["exp"]] <- t(assay(exp_sce))
    # adata$uns[["assays"]][[assays2[i]]][["exp"]] <- t(assay(exp_sce, assays2[i]))
    adata$uns[["assays"]][[assays2[i]]][["feature"]] <- feature_names
    adata$uns[["assays"]][[assays2[i]]][["cells"]] <- cell_names

    ### add motif data
    if (export_pwm && length(getPeakAnnotation(object, name = "Motif")) > 0 && assay.types2[i] == "motif") {
      motif_obj <- getPeakAnnotation(object, name = "Motif")
      motif_pwm_obj <- motif_obj[["motifs"]]
      motif_pwm_list <- lapply(names(motif_pwm_obj), function(x) {
        motif_pwm_obj[[x]]@`profileMatrix`
      })

      names(motif_pwm_list) <- names(motif_pwm_obj)
      motif_pwm_json <- jsonlite::toJSON(motif_pwm_list, auto_unbox = TRUE)

      adata$uns[["assays"]][[assays2[i]]][["motif_pwm"]] <- motif_pwm_json
    }
  }

  return(adata)
}





# Function used to covert archr object into anndata object
#' @title ArchRToAnndata
#' @name Function ArchRToAnndata
#' @param object ArchR object used to export to SGS cell browser
#' @param outpath outputpath of the anndata object
#' @param assays vector of assay names to export
#' @param assay.types vector of the type of assays, mainly peak, gene, motif
#' @param markersDF a named list of marker df vars
#' @param groups vector of cell metadata colnames to export
#' @param reductions vector of reduction names to export
#' @param export_links whether to expport coaccess link informs
#' @param export_pwm whether to export pwm informas
#' @importFrom ArchR getCellColData getMatrixFromProject getPeakAnnotation getCoAccessibility
#' @importFrom SummarizedExperiment rowRanges rowData colData assay seqnames start end
#' @importFrom anndata AnnData write_h5ad
#' @importFrom jsonlite toJSON
#' @export
#' @examples
#' @useDynLib SgsAnndata
#'
ArchRToAnndata <- function(object,
                           outpath,
                           assays,
                           assay.types,
                           markersDF = NULL,
                           groups = NULL,
                           reductions = NULL,
                           export_links = FALSE,
                           export_pwm = FALSE) {
  if ("peak" %in% assay.types) {
    ## delete the assay name
    idx <- match("peak", assay.types)
    if (length(idx) > 1) {
      idx.used <- c[1]
    } else {
      idx.used <- idx
    }

    adata <- Build_init_adata(assay.types, assays, idx.used, object, groups)

    ##### gain peak links
    if (!is.null(export_links) && length(getCoAccessibility(object)) > 0) {
      co_acc <- data.frame(getCoAccessibility(object))
      co_data <- data.frame(
        seqnames = co_acc$`seqnames`,
        start = co_acc$`start`,
        end = co_acc$`end`,
        width = co_acc$`width`,
        strand = co_acc$`strand`,
        score = co_acc$`value`,
        group = co_acc$`group`
      )

      adata$uns[["assays"]][[assays[idx.used]]][["coacc_link"]] <- co_data
    }

    ## add other exp data
    adata <- Add_other_exp(adata, idx.used, assays, assay.types, object, export_pwm)
  } else {
    message("no peak matrix provided")
    idx.used <- 1
    adata <- Build_init_adata(assay.types, assays, idx.used, object, groups)
    adata <- Add_other_exp(adata, idx.used, assays, assay.types, object, export_pwm)
  }

  ######## cell embedding
  dr <- object@`embeddings`
  if (!is.null(reductions) && all(reductions %in% names(dr))) {
    ReducNames <- reductions
  } else if (length(dr) >= 1) {
    ReducNames <- names(dr)
    message("Using all embeddings contained in the Seurat object: ", reducNames)
  } else {
    message("please check the reduction imformation and run again")
    # stop()
  }

  for (embedding in ReducNames) {
    adata$obsm[[embedding]] <- getEmbedding(ArchRProj = object, embedding = embedding, returnDF = TRUE)
  }

  ######## add marker data
  if (!is.null(markersDF)) {
    for (i in seq(length(markersDF))) {
      df_names <- names(markersDF)
      adata$uns[["assays"]][[df_names[i]]][["marker"]] <- markersDF[[i]]
    }
  } else {
    message("no marker provided")
  }


  ### write the adata into h5ad
  # adata$write_h5ad(file.path(outdir,sprintf("siganc.adata") ))
  adata$write_h5ad(outpath)
}
