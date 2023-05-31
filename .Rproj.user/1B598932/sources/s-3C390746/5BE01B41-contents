## This script used to cover seurat object into anndata object
## 2023.5.17
require(Seurat)
require(anndata)
require(Matrix)
require(SeuratObject)
options(Seurat.object.assay.version = "v5")


# Function used to build the init adata object
#' @title Build_init_adata
#' @name Build_init_adata
#' @param object seurat object
#' @param assay name of assay
#' @param groups vector of groups to export,if null export all
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom SeuratObject Features
#' @importFrom anndata AnnData
#' @return adata
#'
Build_init_adata <- function(object, assay, groups) {
  DefaultAssay(object) <- assay
  feature_meta <- data.frame(feature = Features(object))
  rownames(feature_meta) <- feature_meta$`feature`

  if (!is.null(groups) && groups %in% colnames(object[[]])) {
    cell_meta <- object[[]]
    cell_meta <- cell_meta[, which(colnames(cell_meta) %in% groups)]
  } else {
    cell_meta <- object[[]]
  }

  ### be attention
  if (class(object@`assays`[[assay]]) == "Assay") {
    if (!is.null(object@`assays`[[assay]][["data"]])) {
      adata <- AnnData(
        X = t(GetAssayData(object = object, assay = assay, layer = "data")),
        obs = cell_meta,
        var = feature_meta
      )
    } else {
      adata <- AnnData(
        X = t(GetAssayData(object = object, assay = assay, layer = "counts")),
        obs = cell_meta,
        var = feature_meta
      )
    }
  } else if (class(object@`assays`[[assay]]) == "Assay5") {
    options(Seurat.object.assay.version = "v5")

    if (!is.null(object@`assays`[[assay]]@layers[["data"]])) {
      adata <- AnnData(
        X = t(GetAssayData(object = object, assay = assay, layer = "data")),
        obs = cell_meta,
        var = feature_meta
      )
    } else {
      adata <- AnnData(
        X = t(GetAssayData(object = object, assay = assay, layer = "counts")),
        obs = cell_meta,
        var = feature_meta
      )
    }
  }
  return(adata)
}




###### Function used to parser spatial informs of different ST-seq
#' @title Gain_Spatial
#' @name Gain_Spatial
#' @param object seurat object
#' @param img.name name of the image
#' @param adata adata object
#' @return a list
#' @importFrom Seurat GetTissueCoordinates ScaleFactors
#'
Gain_Spatial <- function(object, img.name, adata) {

  ## "SlideSeq","VisiumV1","STARmap"(seqbased);"FOV"(image-based)
  ## class(object@`images`[[img.name]]) == "VisiumV1"
  if (class(object[[img.name]]) == "VisiumV1") {

    ## merge multi slice coords into one

    spatial_coords <- GetTissueCoordinates(object = object, image = img.name)[, c("imagecol", "imagerow")]


    ## gain the related library id
    library_id <- data.frame(library_id = rep(img.name, nrow(spatial_coords)), cell_id = rownames(spatial_coords))
    rownames(library_id) <- library_id$`cell`

    scale <- as.list(ScaleFactors(object[[img.name]]))
    scale_list <- list(
      "tissue_hires_scalef" = scale$`hires`,
      "tissue_lowres_scalef" = scale$`lowres`,
      "fiducial_diameter_fullres" = scale$`fiducial`,
      "spot_diameter_fullres" = scale$`spot`
    )
    img_arrary <- object@`images`[[img.name]]@`image`

    adata$uns[["spatial"]][[img.name]][["images"]][["hires"]] <- img_arrary * 255.0
    # adata$uns[["spatial"]][[img.name]][["images"]][["hires"]] <- apply(img_arrary*255.0, c(3),FUN = as.integer)
    adata$uns[["spatial"]][[img.name]][["scalefactors"]] <- scale_list
  } else if (class(object[[img.name]]) == "SlideSeq") {

    ## merge multi slice coords into one
    spatial_coords <- GetTissueCoordinates(object = object, image = img.name)[, c("y", "x")]

    ## gain the related library id
    library_id <- data.frame(library_id = rep(img.name, nrow(spatial_coords)), cell_id = rownames(spatial_coords))
    rownames(library_id) <- library_id$`cell_id`
  } else if (class(object[[img.name]]) == "STARmap") {

    ## merge multi slice coords into one
    spatial_coords <- GetTissueCoordinates(object = object, image = img.name)[, c("y", "x")]

    ## gain the related library id
    library_id <- data.frame(library_id = rep(img.name, nrow(spatial_coords)), cell_id = rownames(spatial_coords))
    rownames(library_id) <- library_id$`cell_id`
  } else if (class(object@`images`[[img.name]]) == "FOV") {

    ## merge multi slice coords into one
    spatial_coords <- GetTissueCoordinates(object = object[[img.name]], which = "centroids")
    rownames(spatial_coords) <- spatial_coords$`cell`
    spatial_coords <- spatial_coords[, c("y", "x")]

    # spatial_coords <- GetTissueCoordinates(object = object[[img.name]], which = "centroids")[ ,c("x","y")]
    ## gain the related library id
    library_id <- data.frame(library_id = rep(img.name, nrow(spatial_coords)), cell_id = rownames(spatial_coords))
    rownames(library_id) <- library_id$`cell_id`


    ## deal the segment informations
    if ("segmentation" %in% names(object[[img.name]])) {
      ## get segmentation coords, becareful, the segament data dim is diferent to tissue coords
      Segmentation_coords <- GetTissueCoordinates(object = object[[img.name]], which = "segmentation")
      Segmentation_coords <- as.matrix(Segmentation_coords[, c("x", "y")])
      ### add the segamnet data into uns
      adata$uns[["spatial"]][[img.name]][["images"]][["Segmentation_coords"]] <- Segmentation_coords
    }
  } else {
    message("spatial information fetech failed")
    stop()
  }

  result <- list(adata = adata, sp_coords = spatial_coords, library_id = library_id)
  return(result)
}







# Function used to ccovert seurat object into anndata(multi assay)
#' @title SeuratToAnndata
#' @name Function SeuratToAnndata
#' @param object seurat object
#' @param outpath outputpath of the anndata object
#' @param main.assays name of assay used to build adata.X
#' @param assays a vector of assay names to export in adata
#' @param groups vector of groups to export,if null export all
#' @param reductions vector of reduction names to export,if null export all
#' @param images vector of images names to export,if null export all
#' @param markersDF a dataframe of the marker files,need the colnames in c(gene, cluster, ...) order
#' @importFrom Seurat GetAssayData GetTissueCoordinates ScaleFactors  DefaultAssay GetAssayData Cells Images Embeddings
#' @importFrom SeuratObject Features
#' @importFrom anndata AnnData write_h5ad
#' @export
#' @examples
#' @useDynLib SgsAnndata
#'
SeuratToAnndata <- function(object,
                            outpath,
                            main.assays,
                            assays = NULL,
                            groups = NULL,
                            reductions = NULL,
                            images = NULL,
                            markersDF = NULL) {

  #### this ccovert function mainly for seurat V3-seurat V5
  message("Seurat Version installed: ", packageVersion("Seurat"))
  message("Object was created with Seurat version ", object@version)


  ##### build the init adata and add related assay
  adata <- Build_init_adata(object = object, assay = main.assays, groups = groups)
  if (!is.null(assays) && length(assays) >= 1) {
    for (i in seq(length(assays))) {
      DefaultAssay(object) <- assays[i]
      if (!is.null(object@assays[[assays[i]]][["data"]])) {
        adata$uns[["assays"]][[assays[i]]][["exp"]] <- t(GetAssayData(object = object, assay = assays[i], layer = "data"))
      } else {
        adata$uns[["assays"]][[assays[i]]][["exp"]] <- t(GetAssayData(object = object, assay = assays[i], layer = "counts"))
      }
      adata$uns[["assays"]][[assays[i]]][["feature"]] <- Features(object)
      adata$uns[["assays"]][[assays[i]]][["cells"]] <- Cells(object)
    }
  }


  ##### fetech the spatial information
  # judge the images classes:"SlideSeq","VisiumV1","STARmap","FOV"
  # class(slide.seq@images$image)
  DefaultAssay(object) <- main.assays
  if (!is.null(images)) {
    imageNames <- images
  } else {
    if (length(Images(object = object, assay = main.assays)) > 0) {
      images <- Images(object = object, assay = main.assays)
    } else {
      images <- Images(object = object)
    }
    # images <- Images(object = object, assay = assay)
    imageNames <- images
  }

  if (length(imageNames) >= 1 && class(object[[imageNames[1]]]) == "VisiumV1") {

    # for (img in imageNames) {
    #   adata <- Gain_Spatial(object = object, img.name = img, adata = adata)
    # }
    ## changed
    all_coords <- data.frame()
    library_ids <- data.frame()
    for (img in imageNames) {
      result <- Gain_Spatial(object = object, img.name = img, adata = adata)
      spatial_coords <- result[["sp_coords"]]
      all_coords <- rbind(spatial_coords, all_coords)

      library_id <- result[["library_id"]]
      library_ids <- rbind(library_id, library_ids)
    }

    # (adata=adata, sp_coords = spatial_coords, library_id = library_id)
    all_coords <- all_coords[rownames(object[[]]), ]
    library_ids <- library_ids[rownames(object[[]]), ]
    adata <- result[["adata"]]
    adata$obsm[["spatial"]] <- as.matrix(all_coords)
    adata$obs[["library_id"]] <- library_ids$`library_id`
  } else if (length(imageNames) >= 1) {
    for (img in imageNames) {
      result <- Gain_Spatial(object = object, img.name = img, adata = adata)
    }
    spatial_coords <- result[["sp_coords"]]
    library_id <- result[["library_id"]]
    adata <- result[["adata"]]
    adata$obsm[["spatial"]] <- as.matrix(spatial_coords)
    adata$obs[["library_id"]] <- library_id$`library_id`
  } else {
    message("Object has no image informations")
  }


  ###### cell embeddings
  dr <- object@`reductions`
  if (!is.null(reductions) && all(reductions %in% names(dr))) {
    ReducNames <- reductions
  } else if (length(dr) >= 1) {
    ReducNames <- names(dr)
    message("Using all embeddings contained in the Seurat object:")
  } else {
    message("please check the reduction imformation and run again")
    # stop()
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
