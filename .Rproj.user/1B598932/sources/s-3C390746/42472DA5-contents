## This script used to covert Giotto object into anndata object
## 2023.5.10
require(anndata)
require(Giotto)
require(magick)
require(imager)
require(png)



#' @title list_expression_names
#' @name list_expression_names
#' @description lists the available matrices names for a given spatial unit and feature type
#' @param gobject giotto object
#' @param spat_unit spatial unit ‘cell’
#' @param feat_type feature type ‘rna’
#' @return vector with names of available matrices
list_expression_names <- function(gobject,
                                  spat_unit = NULL,
                                  feat_type = NULL) {
  if (is.null(spat_unit)) stop("spat_unit must be given\n")
  if (is.null(feat_type)) stop("feat_type must be given\n")

  expression_names <- names(gobject@expression[[spat_unit]][[feat_type]])
  return(expression_names)
}



#' @title set_default_feat_type
#' @name set_default_feat_type
#' @param gobject giotto object
#' @param feat_type feature type ‘rna’
#' @param spat_unit spatial unit ‘cell’
#' @keywords internal
set_default_feat_type <- function(gobject,
                                  feat_type = NULL,
                                  spat_unit) {

  # set spatial unit
  if (is.null(feat_type)) {
    feat_type <- getOption("giotto.feat_type")

    if (is.null(feat_type)) {
      if (!is.null(gobject@`expression`) & length(gobject@`expression`) > 0) {
        feat_type <- names(gobject@`expression`[[spat_unit]])[[1]]
        if (is.null(feat_type)) stop("valid spat_unit input needed \n")
      } else if (!is.null(gobject@`feat_info`)) {
        feat_type <- names(gobject@`feat_info`)[[1]]
      } else {
        warning("No default for feat_type could be set \n")
      }
    }
  }

  return(feat_type)
}


#' @title set_default_spat_unit
#' @name set_default_spat_unit
#' @param gobject gobject
#' @param spat_unit spatial unit ‘cell’
#' @keywords internal
set_default_spat_unit <- function(gobject,
                                  spat_unit = NULL) {

  # set spatial unit
  if (is.null(spat_unit)) {
    spat_unit <- getOption("giotto.spat_unit")
    if (is.null(spat_unit)) {
      if (!is.null(gobject@`expression`) & length(gobject@`expression`) > 0) {
        spat_unit <- names(gobject@`expression`)[[1]]
      } else if (!is.null(gobject@`spatial_info`)) {
        spat_unit <- names(gobject@`spatial_info`)[[1]]
      } else {
        warning("No default for spat_unit could be set \n")
      }
    }
  }

  return(spat_unit)
}



### expression values slot ####
#' @title Get expression values
#' @name get_expression_values
#' @description Function to get expression values from giotto object
#' @param gobject giotto object used to make anndata
#' @param values expression values to extract (e.g. "raw", "normalized", "scaled")
#' @param spat_unit spatial unit ‘cell’
#' @param feat_type feature type ‘rna’
#' @param output what object type to retrieve the expression as. Currently either
#' 'matrix' for the matrix object contained in the exprObj or 'exprObj' (default) for
#' the exprObj itself are allowed.
#' @return expression matrix

get_expression_values <- function(gobject,
                                  values = NULL,
                                  spat_unit = NULL,
                                  feat_type = NULL,
                                  output = c("exprObj", "matrix")) {
  output <- match.arg(output, choices = c("exprObj", "matrix"))

  # 1. Set feat_type and spat_unit

  spat_unit <- set_default_spat_unit(
    gobject = gobject,
    spat_unit = spat_unit
  )
  feat_type <- set_default_feat_type(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type
  )

  # 2. Find object

  potential_values <- list_expression_names(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type
  )

  if (is.null(values)) values <- potential_values[[1]]

  ## special cases for giotto standard pipeline
  if (values == "scaled" & is.null(gobject@`expression`[[spat_unit]][[feat_type]][[values]])) {
    stop("run first scaling (& normalization) step")
  } else if (values == "normalized" & is.null(gobject@`expression`[[spat_unit]][[feat_type]][[values]])) {
    stop("run first normalization step")
  } else if (values == "custom" & is.null(gobject@`expression`[[spat_unit]][[feat_type]][[values]])) {
    stop("first add custom expression matrix")
  }

  if (!values %in% potential_values) stop("The spatial unit ", spat_unit, " for expression matrix ", feat_type, " and with name ", "'", values, "'", " can not be found \n")

  # 3. Get object in desired format

  expr_values <- gobject@`expression`[[spat_unit]][[feat_type]][[values]]

  if (output == "exprObj") {
    if (!inherits(expr_values, "exprObj")) {
      expr_values <- new("exprObj",
        name = values,
        exprMat = expr_values,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = spat_unit, # assumed
        misc = NULL
      )
    }

    if (!inherits(expr_values, "exprObj")) stop("Cannot convert to exprObj")

    # return exprObj
    return(expr_values)
  } else if (output == "matrix") {
    if (inherits(expr_values, "exprObj")) expr_values <- slot(expr_values, "exprMat")

    # return 'matrix'
    return(expr_values)
  }
}




#' @title GiottoToAnnData
#' @name Function GiottoToAnnData
#' @description Coveret Giotto To Anndata
#' @param object giotto object used to make anndata
#' @param outpath anndata output path
#' @param assays name of assays to export
#' @param reductions a vector of reduction names used to export
#' @param images name of image used to export
#' @param markersDF a named list of marker df vars
#' @importFrom anndata write_h5ad AnnData
#' @importFrom Giotto get_spatial_locations pDataDT get_dimReduction fDataDT
#' @importFrom png readPNG
#' @export
#' @useDynLib SgsAnndata
#'
GiottoToAnnData <- function(object = NULL,
                            outpath = NULL,
                            assays = NULL,
                            reductions = NULL,
                            images = NULL,
                            markersDF = NULL) {

  # object <- visium_brain
  # Check gobject
  invalid_obj <- !("giotto" %in% class(object))
  if (is.null(object) || invalid_obj) {
    stop(wrap_msg("Please provide a valid Giotto Object for conversion."))
  }

  # Check outpath directory, make it if it doesn't exist
  if (is.null(outpath)) outpath <- paste0(getwd(), "/GiottoToAnndata.h5ad")

  ## raw_x:a sparse matrix with rows named genes and columns named cells（"dgTMatrix"）；
  raw_x <- get_expression_values(
    gobject = object,
    values = "raw",
    spat_unit = "cell",
    feat_type = "rna",
    output = "matrix"
  )


  # fetech cell annotation meta informations
  # cell_meta_anno <- data.frame(visium_brain@cell_metadata$cell$rna@metaDT)
  cell_meta_anno <- pDataDT(object, spat_unit = NULL, feat_type = NULL)
  rownames(cell_meta_anno) <- cell_meta_anno$`cell_ID`

  # fetech feature annotation meta informations
  gene_metadata <- fDataDT(object, spat_unit = NULL, feat_type = NULL)
  rownames(gene_metadata) <- gene_metadata$`feat_ID`


  if (!is.null(raw_x)) {
    adata <- AnnData(
      X = t(raw_x),
      obs = cell_meta_anno,
      var = gene_metadata
    )
  } else {
    message("raw_x information fetech failed, please make sure it is not null! ")
  }



  ## fetach space location information
  ## visium_brain@spatial_locs$cell$raw@coordinates
  # spatial_coordinates <- visium_brain@spatial_locs$cell$raw@coordinates
  spatial_locs <- get_spatial_locations(
    gobject = object,
    output = "data.table",
    spat_unit = "cell"
  )


  if (!is.null(spatial_locs)) {
    # img_spatial <- sprintf("%s_spatial",visium_brain@images$image@name)
    # dim(obsm$spatial); obsm$spatial  [1]    2 2698
    adata$obsm[["spatial"]] <- as.matrix(spatial_locs[, c("sdimx", "sdimy")])
  } else {
    message("spatial_locs information fetech failed, please make sure it is not null! ")
  }


  ##### fetch dimension reduction information
  if (!is.null(object@dimension_reduction)) {
    # get a list of dimension reduction information
    dr <- object@dimension_reduction$`cells`$`cell`$`rna`
    if (length(dr) >= 1) {
      ReducNames <- list(names(dr))
      # ReducNames[[1]]
      for (embedding in ReducNames[[1]]) {
        adata$obsm[[embedding]] <- get_dimReduction(
          object,
          spat_unit = NULL,
          feat_type = NULL,
          reduction = c("cells", "feats"),
          reduction_method = embedding,
          name = embedding,
          output = "data.table",
          set_defaults = TRUE
        )
        message("Dimension reduction information ", embedding, " is successfully stored in anndata object.")
      }
    }
  } else {
    message("dimension_reduction is null!")
  }


  ############ Image storage processing
  ###### fetech image slice information
  if (!is.null(object@`images`)) {
    # Get a list of image information
    dr <- object@`images`
    if (length(dr) >= 1) {
      ImagesNames <- list(names(dr))

      for (image_name in ImagesNames[[1]]) {

        # obtaining image path
        image_path <- object@images[[image_name]]@file_path
        # read PNG or TIFF images
        img <- readPNG(image_path)
        # Convert the image to an RGB three-channel array
        img_arrary <- as.array(img)
        dim(img_arrary) # output something like [1] width x height x 3
        # Stored image information
        adata$uns[["spatial"]][[image_name]][["images"]][["hires"]] <- img_arrary * 255.0

        message("Images information ----- ", image_name, " is successfully stored in anndata object.")
      }
    }
  } else {
    message("images is null!")
  }

  ### write the adata into h5ad
  adata$write_h5ad(outpath)

  # adata
}
