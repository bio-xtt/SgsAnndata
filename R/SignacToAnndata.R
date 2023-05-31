## This script used to cover siganc object into anndata object
## 2023.4.10
require(Signac)
require(Seurat)
require(jsonlite)
require(anndata)
require(Matrix)
require(parallel)
require(SeuratObject)
options(Seurat.object.assay.version = "v5")


# Function used to build adata according to the assay type,when multi assay provided
#' @title Build_init_adata
#' @name Build_init_adata
#' @param assay.types vector of the type of assays, mainly peak, gene, motif
#' @param assays vector of assay names to export
#' @param idx.used assay name
#' @param object siganc object
#' @param groups vector of cell metadata colnames to export
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom SeuratObject Features
#' @importFrom anndata AnnData
#' @return adata
#'
Build_init_adata <- function(assay.types, assays, idx.used, object, groups) {
  ## fetech metadata
  if (!is.null(groups) && groups %in% colnames(object[[]])) {
    cell_meta <- object[[]]
    cell_meta <- cell_meta[, which(colnames(cell_meta) %in% groups)]
  } else {
    cell_meta <- object[[]]
  }

  ###### fetech data
  ## fetech feature x cell matrix and add related data
  DefaultAssay(object) <- assays[idx.used]
  feature_meta <- data.frame(feature = Features(object))
  rownames(feature_meta) <- feature_meta[["feature"]]

  ## attension the exp data fetch method between seurat v3 ~ seurat v5
  if (!is.null(object@assays[[assays[idx.used]]][["data"]])) {
    adata <- AnnData(
      X = t(GetAssayData(object = object, assay = assays[idx.used], layer = "data")),
      obs = cell_meta,
      var = feature_meta
    )
  } else {
    adata <- AnnData(
      X = t(GetAssayData(object = object, assay = assays[idx.used], layer = "counts")),
      obs = cell_meta,
      var = feature_meta
    )
  }

  return(adata)
}



# Function used to add other exp and related data into adata uns
#' @title Add_other_exp
#' @name Add_other_exp
#' @param adata adata object
#' @param idx.used assay name
#' @param assays vector of assay names to export
#' @param assay.types vector of the type of assays, mainly peak, gene, motif
#' @param object siganc object
#' @param export_pwm whether to export pwm informas
#' @importFrom Signac GetMotifData
#' @importFrom Seurat GetAssayData DefaultAssay Cells
#' @importFrom SeuratObject Features
#' @importFrom jsonlite toJSON
#' @return adata
#'
Add_other_exp <- function(adata, idx.used, assays, assay.types, object, export_pwm) {

  ## delete the used assays
  assays2 <- assays[-idx.used]
  assay.types2 <- assay.types[-idx.used]

  ## add other exp data into adata
  for (i in seq(length(assay.types2))) {
    DefaultAssay(object) <- assays2[i]
    if (!is.null(object@assays[[assays2[idx.used]]][["data"]])) {
      adata$uns[["assays"]][[assays2[i]]][["exp"]] <- t(GetAssayData(object = object, assay = assays2[idx.used], layer = "data"))
    } else {
      adata$uns[["assays"]][[assays2[i]]][["exp"]] <- t(GetAssayData(object = object, assay = assays2[idx.used], layer = "counts"))
    }

    adata$uns[["assays"]][[assays2[i]]][["feature"]] <- Features(object)
    adata$uns[["assays"]][[assays2[i]]][["cells"]] <- Cells(object)

    ### add motif data
    if (export_pwm && length(GetMotifData(object = object, slot = "pwm", assay = assays[idx.used])) > 0 && assay.types2[i] == "motif") {
      pwm_data <- GetMotifData(object = object, slot = "pwm", assay = assays[idx.used])
      motif_name <- GetMotifData(object = object, slot = "motif.names", assay = assays[idx.used])
      new_pwm_name <- as.vector(paste(names(pwm_data), motif_name, sep = "_"))
      names(pwm_data) <- new_pwm_name
      pwm_json <- jsonlite::toJSON(pwm_data, auto_unbox = TRUE)
      adata$uns[["assays"]][[assays2[i]]][["motif_pwm"]] <- pwm_json
    }
  }
  return(adata)
}





# Function used to covert signac object into anndata object
#' @title SignacToAnndata
#' @name Function SignacToAnndata
#' @param object siganc object
#' @param outpath anndata output path
#' @param assays vector of assay names to export
#' @param assay.types vector of the type of assays, mainly peak, gene, motif
#' @param markersDF a named list of marker df vars
#' @param outpath outputpath of the anndata object
#' @param groups vector of cell metadata colnames to export
#' @param reductions vector of reduction names to export
#' @param export_links whether to expport coaccess link informs
#' @param export_pwm whether to export pwm informas
#' @param bw_dir path of cluster bigwig dirs
#' @importFrom Signac GetMotifData Links
#' @importFrom Seurat GetAssayData DefaultAssay Cells Embeddings
#' @importFrom SeuratObject Features
#' @importFrom jsonlite toJSON
#' @importFrom anndata AnnData write_h5ad
#' @examples
#' @export
#' @useDynLib SgsAnndata
#'
SignacToAnndata <- function(object,
                            outpath,
                            assays,
                            assay.types,
                            markersDF = NULL,
                            groups = NULL,
                            reductions = NULL,
                            export_links = FALSE,
                            export_pwm = FALSE,
                            bw_dir = NULL) {

  ###### exp related data
  if ("peak" %in% assay.types) {
    ## delete the assay name
    idx <- match("peak", assay.types)
    if (length(idx) > 1) {
      idx.used <- c[1]
    } else {
      idx.used <- idx
    }
    adata <- Build_init_adata(assay.types, assays, idx.used, object, groups)

    ### gain peak links
    if (!is.null(export_links) && length(Links(object = object[[assays[idx.used]]])) > 0) {
      co_data <- as.data.frame(Links(object = object[[assays[idx.used]]]))
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
  dr <- object@`reductions`
  if (!is.null(reductions) && all(reductions %in% names(dr))) {
    ReducNames <- reductions
  } else if (length(dr) >= 1) {
    ReducNames <- names(dr)
    message("Using all embeddings contained in the Seurat object: ", reducNames)
  } else {
    message("please check the reduction imformation and run again")
  }

  for (embedding in ReducNames) {
    adata$obsm[[embedding]] <- Embeddings(object = object, embedding)
  }


  ######## add meta data
  if (!is.null(markersDF)) {
    for (i in seq(length(markersDF))) {
      df_names <- names(markersDF)
      adata$uns[["assays"]][[df_names[i]]][["marker"]] <- markersDF[[i]]
    }
  } else {
    message("no marker provided")
  }

  ### write the adata into h5ad
  adata$write_h5ad(outpath)
}
