#' Find alternate SingleCellExperiment containing desired reduction embeddings.
#' 
#' @description  Find which alternate experiment within a SingleCellExperiment object a given reduction can be found in. This is needed if a user wants to retrieve and plot the cell embeddings for a given reduction from a SingleCellExperiment object.
#'
#' @param alt_exp_name A string specifying the name of an alternate experiment from a SingleCellExperiment object in which to search if a given reduction is present. The names of alternate experiments can be obtained by looking through SingleCellExperiment::altExpNames(x), where x is a SingleCellExperiment object.
#' @param sce_object A SingleCellExperiment object containing data from a CITE-seq experiment processed with Seurat.
#' @param reduction_name A string specifying the name of a reduction (e.g., "UMAP", "PCA", etc.).
#'
#' @importFrom SingleCellExperiment altExp reducedDimNames
#'
#' @return A string specifying the name of an alternate experiment (e.g., "ADT", "RNA", "SCT", etc.) in which the cell embeddings for a given reduction can be found.
#' @noRd
#'
#' @examples \dontrun{
#' sce_object <- readRDS("path/to/RDS_file_containing_SCE_object.rds")
#' alt_exp_name <- SingleCellExperiment::altExpNames(sce_object)[1]
#' reduction_name <- "UMAP"
#' alt_exp_of_desired_reduction <- find_reduction_in_altSCE(alt_exp_name, sce_object, reduction_name)
#' }
find_reduction_in_altSCE <- function(alt_exp_name, sce_object, reduction_name) {
  alt_exp_object <- SingleCellExperiment::altExp(x = sce_object, 
                                                 e = alt_exp_name)
  alt_exp_reductions <- SingleCellExperiment::reducedDimNames(alt_exp_object)
  # match the reduction name exactly when searching if it's within a list
  reduction_finder_results <- grepl(paste0("^", reduction_name, "$"), 
                                    alt_exp_reductions)
  if (TRUE %in% reduction_finder_results) {
    return(alt_exp_name)
  }
}


#' Check if alternate experiments in a SingleCellExperiment (SCE) object inherit from a given class.
#'
#' @description By default, alternate experiments in a SingleCellExperiment object will inherit from the SummarizedExperiment class. However, sometimes SCE objects created from Seurat-processed data and converted to SCE within Seurat may have alternate experiments that inherit from both the SCE and SummarizedExperiment classes. This multiple inheritance affects the data slots present in an alternate experiment object as well as the methods that can be used to access data from alternate experiments in a SCE object. This function can be used to check if alternate experiments in a SingleCellExperiment (SCE) object also inherit from a class other than SummarizedExperiment, such as the SingleCellExperiment class. 
#'
#' @param sce_object A SingleCellExperiment object containing CITE-seq data and, by extension, at least one alternate experiment.
#' @param class_name A string containing the name of a class for which class inheritance of alternate experiments in a SingleCellExperiment object is to be tested. 
#'
#' @importFrom SingleCellExperiment altExps
#'
#' @return A Boolean value. If all alternate experiments in a SingleCellExperiment (SCE) object inherit from a given class, then this return value is TRUE. If one or more alternate experiments in a SingleCellExperiment (SCE) object do not inherit from a given class, then this return value is FALSE.
#' @noRd
#'
#' @examples \dontrun{
#' sce_object <- readRDS("path/to/RDS_file_containing_SCE_object.rds")
#' 
#' altExps_inherit_class(sce_object, class_name = "SingleCellExperiment")
#' }
altExps_inherit_class <- function(sce_object, class_name) {
  inherits_from_class <- NULL
  altExp_inheritance_results <- lapply(
    SingleCellExperiment::altExps(sce_object), 
    inherits, 
    what = class_name)
  # check if FALSE is anywhere in the inheritance results list
  # this accounts for cases where some alt experiments inherit from a given class but others don't 
  if (FALSE %in% altExp_inheritance_results) {
    inherits_from_class <- FALSE
  }
  else {
    inherits_from_class <- TRUE
  }
  return(inherits_from_class)
}


# helper function for get_choices
get_choices_from_seurat <- function(category, seurat_object, input_file_df, assay_name = NULL) {
  menu_choices <- NULL
  if (is.null(category) & !is.null(assay_name)) {
    menu_choices <- rownames(SeuratObject::GetAssayData(object = seurat_object, 
                                                        slot = "data",
                                                        assay = assay_name))
  }
  else if (category == "metadata") {
    menu_choices <- seurat_object[[]] %>%
      dplyr::select(where(is.factor) | where(is.character)) %>% 
      colnames()
  }
  else if (category == "reductions") {
    menu_choices <- SeuratObject::Reductions(seurat_object)
  }
  else if (category == "assays") {
    menu_choices <- SeuratObject::Assays(object = seurat_object)
  }
  return(menu_choices)
}


# helper function for get_choices
get_choices_from_sce <- function(category, sce_object, input_file_df, assay_name = NULL) {
  menu_choices <- NULL
  if (is.null(category) & !is.null(assay_name)) {
    if (SingleCellExperiment::mainExpName(sce_object) == assay_name) {
      menu_choices <- rownames(sce_object)
    }
    else {
      menu_choices <- rownames(SingleCellExperiment::altExp(x = sce_object,
                                                            e = assay_name))
    }
  }
  else if (category == "metadata") {
    menu_choices <- as.data.frame(SingleCellExperiment::colData(sce_object)) %>% 
      dplyr::select(where(is.factor) | where(is.character)) %>%
      colnames()
  }
  else if (category == "reductions") {
    if (altExps_inherit_class(sce_object, "SingleCellExperiment")) {
      menu_choices <- unlist(
        SingleCellExperiment::applySCE(sce_object, 
                                       SingleCellExperiment::reducedDimNames))
    }
    else {
      menu_choices <- SingleCellExperiment::reducedDimNames(sce_object)
    }
  }
  else if (category == "assays") {
    menu_choices <- c(SingleCellExperiment::mainExpName(sce_object), 
                      SingleCellExperiment::altExpNames(sce_object))
  }
  return(menu_choices)
}



#' Get dropdown menu options for selectInput elements in CITEViz.
#' 
#' @description  This function retrieves various categories of data (e.g., metadata, reductions, or assays) from user-uploaded input file(s) (e.g., an RDS file containing a Seurat object) and returns a character vector of items with which to populate a dropdown menu in the CITEViz user interface. This function can also get subchoices for a dropdown menu after a user has selected an assay they want to view. For example, if a user selects the "ADT" assay, then this function will return a vector of all the possible ADTs a user can choose to view from their input data. If a user selects the "RNA" assay, then this function will return a vector of all the genes a user can choose to view from their input data.
#'
#' @param category A string specifying the category of data to retrieve. Categories are based on slots in a Seurat object produced from a CITE-seq experiment Possible values are "metadata", "reductions", or "assays". Can be NULL if the user wants to retrieve data for a specific assay (e.g. if the user selects the "ADT" assay and wants a sub-menu of specific ADTs for which they can plot/view expression data)
#' @param input_data_type An integer value indicating if the user-uploaded input files are RDS files holding Seurat objects or SingleCellExperiment objects (1 = Seurat object, 2 = SingleCellExperiment object, etc.)
#' @param rds_object An RDS object containing metadata, assays, and reductions for a CITE-seq experiment. Can be NULL if the user uploaded an alternate supported file instead of an RDS file.
#' @param input_file_df A dataframe generated by retrieving file input data from a fileInput element in the UI (i.e., by running input$file_input in the app server). This dataframe contains information about each file uploaded by the user, such as the name of each file as it appears on the user's local filesystem (input_file_df$name) and the temporary datapath with which data from the file can be read in (input_file_df$datapath).
#' @param assay_name A string specifying the name of an assay (e.g., "ADT", "RNA", etc.) from a Seurat-processed CITE-seq experiment.
#' 
#' @import magrittr
#' @importFrom SeuratObject Assays Reductions GetAssayData
#' @importFrom dplyr select
#' @importFrom SingleCellExperiment altExp altExpNames applySCE colData mainExpName reducedDimNames
#'
#' @return A character vector of items with which to populate a dropdown menu in the CITEViz UI, sorted in ascending order.
#' @noRd
#'
#' @examples \dontrun{
#' # if RDS file is uploaded by user
#' myso <- readRDS("path/to/RDS_file.rds")
#' input_file_data <- input$file_input
#' 
#' get_choices("metadata", 
#' input_data_type = FALSE, 
#' rds_object = myso, 
#' input_file_df = input_file_data, 
#' assay_name = NULL)
#' }
#' 
get_choices <- function(category, input_data_type, rds_object, input_file_df, assay_name = NULL) {
  menu_choices <- NULL
  # if input data type is a Seurat object from an RDS file
  if (input_data_type == 1  & !is.null(rds_object)) {
    menu_choices <- get_choices_from_seurat(category, rds_object, input_file_df, assay_name)
  }
  # if input data type is a SingleCellExperiment object from an RDS file
  else if (input_data_type == 2  & !is.null(rds_object)) {
    menu_choices <- get_choices_from_sce(category, rds_object, input_file_df, assay_name)
  }
  return(sort(as.vector(menu_choices)))
}



# helper function for get_data
get_data_from_seurat <- function(category, seurat_object, input_file_df, assay_name = NULL, reduction_name = NULL, assay_data_to_get = NULL) {
  data <- NULL
 
  if (category == "metadata") {
    # get all metadata in seurat object
    data <- seurat_object[[]]
  }
  else if (category == "assays" & !is.null(assay_name)) {
    # get count data from seurat object for a given assay (ADT, RNA, etc) 
    # and transpose it so that colnames are RNAs/ADTs/etc and rownames are cell barcodes, 
    # for consistency with the way assay data is accessed in the app server
    data <- t(as.data.frame(
      SeuratObject::GetAssayData(object = seurat_object, 
                                 slot = "data",
                                 assay = assay_name))[assay_data_to_get, ])
  }
  else if (category == "reductions" & !is.null(reduction_name)) {
    # we need, at most, 3 reduction dimensions for plotting, so don't get more than 3 dimensions
    # if a reduction doesn't have 3 dimensions (i.e. if PCA data only contains PC1 and PC2),
    # then just get whatever is there
    all_data <- SeuratObject::Embeddings(object = seurat_object, 
                                         reduction = reduction_name)
    num_components <- ncol(all_data)
    if (num_components < 3) {
      data <- all_data
    }
    else {
      data <- all_data[, seq(1, 3)]
    }
  }
  return(data)
}


# helper function for get_data
get_data_from_sce <- function(category, sce_object, input_file_df, assay_name = NULL, reduction_name = NULL, assay_data_to_get = NULL) {
  
  data <- NULL
  if (category == "metadata") {
    data <- SingleCellExperiment::colData(sce_object)
  }
  else if (category == "assays" & !is.null(assay_name)) {
    count_datatype_to_get <- "logcounts"
    experiment_obj <- NULL
    if (SingleCellExperiment::mainExpName(sce_object) == assay_name) {
      experiment_obj <- sce_object
    }
    else {
      experiment_obj <- SingleCellExperiment::altExp(x = sce_object, 
                                                     e = assay_name)
    }
    count_datatypes_list <- SummarizedExperiment::assayNames(experiment_obj)
    if (!("logcounts" %in% count_datatypes_list) & ("normcounts" %in% count_datatypes_list)){
      count_datatype_to_get <- "normcounts"
    }
    data <- t(as.matrix(SummarizedExperiment::assay(x = experiment_obj,
                                                    i = count_datatype_to_get)))
  }
  else if (category == "reductions" & !is.null(reduction_name)) {
    main_exp_reductions <- SingleCellExperiment::reducedDimNames(sce_object)
    # match the reduction name exactly when searching if it's within a list
    reduction_finder_results <- grepl(pattern = paste0("^", reduction_name, "$"), 
                                      x = main_exp_reductions)
    if (TRUE %in% reduction_finder_results) {
      all_data <- SingleCellExperiment::reducedDim(sce_object, 
                                                   type = reduction_name, 
                                                   withDimnames = FALSE)
    }
    # else if the alt experiments inherit from SCE class and, thus, contain reduction embeddings data that is stored somewhere other than in the main experiment of the parent SCE object
    else if (altExps_inherit_class(sce_object, "SingleCellExperiment")){
      # find which alternate experiment a reduction can be found under
      alt_exp_names_list <- SingleCellExperiment::altExpNames(sce_object)
      assay_name <- lapply(alt_exp_names_list, 
                           FUN = find_reduction_in_altSCE, 
                           sce_object = sce_object, 
                           reduction_name = reduction_name) %>% unlist()
      alt_exp_obj <- SingleCellExperiment::altExp(x = sce_object,
                                                  e = assay_name)
      all_data <- SingleCellExperiment::reducedDim(alt_exp_obj,
                                                   type = reduction_name, 
                                                   withDimnames = FALSE)
    }
    num_components <- ncol(all_data)
    if (num_components < 3) {
      data <- all_data
    }
    else {
      data <- all_data[, seq(1, 3)]
    }
  }
  return(data)
}




#' Get CITE-seq data from an object read in form an RDS file. 
#' 
#' @description This function gets metadata, reduction embeddings (e.g. UMAP embeddings, PCA embeddings, etc.), or normalized assay expression levels (e.g., normalized ADT or RNA/gene expression, etc.) from an object read in from an RDS file.
#'
#' @param category A string specifying the category of data to retrieve. Categories are based on slots in a Seurat object produced from a CITE-seq experiment Possible values are "metadata", "reductions", or "assays".
#' @param input_data_type An integer value indicating if the user-uploaded input files are RDS files holding Seurat objects or SingleCellExperiment objects (1 = Seurat object, 2 = SingleCellExperiment object, etc.)
#' @param rds_object An RDS object containing metadata, assays, and reductions for a CITE-seq experiment. Can be NULL if the user uploaded an alternate supported file instead of an RDS file.
#' @param input_file_df A dataframe generated by retrieving file input data from a fileInput element in the UI (i.e., by running input$file_input in the app server). This dataframe contains information about each file uploaded by the user, such as the name of each file as it appears on the user's local filesystem (input_file_df$name) and the temporary datapath with which data from the file can be read in (input_file_df$datapath).
#' @param assay_name A string specifying the name of an assay (e.g., "ADT", "RNA", etc.) from a Seurat-processed CITE-seq experiment. Can be NULL if not retrieving any assay data.
#' @param reduction_name A string specifying the name of an reduction (e.g., "pca", "umap", etc.) from a Seurat-processed CITE-seq experiment. Can be NULL if not retrieving any reduction data.
#' @param assay_data_to_get Character vector specifying the column(s) from assay data to get (i.e. if the assay is "ADT", then a possible value for assay_data_to_get may be c("ab-CD10", "ab-CD38")). Can be NULL if not retrieving any assay data.
#' 
#' @import magrittr
#' @importFrom SeuratObject Embeddings GetAssayData
#' @importFrom SingleCellExperiment altExp colData logcounts mainExpName reducedDimNames reducedDim
#' @importFrom SummarizedExperiment assay assayNames
#' 
#' @return A dataframe containing metadata, assay count data, or reduction embeddings data that was generated in a Seurat-processed CITE-seq experiment. The rownames of the dataframe are cell barcodes, and the column names are metadata columns, assay columns (i.e. if assay is "RNA", then assay columns would be genes), or reduction embedding columns (i.e. if the reduction is "pca", the embedding data columns would be PC1, PC2, etc).
#' @noRd
#'
#' @examples \dontrun{
#' # if RDS file is uploaded by user
#' myso <- readRDS("path/to/RDS_file.rds")
#' input_file_data <- input$file_input
#' 
#' get_data("reductions", 
#' input_data_type = FALSE, 
#' rds_object = myso, 
#' input_file_df = input_file_data, 
#' reduction_name = "pca")
#' }
#' 
get_data <- function(category, input_data_type, rds_object, input_file_df, assay_name = NULL, reduction_name = NULL, assay_data_to_get = NULL) {
  data <- NULL
  # if input data type is a Seurat object from an RDS file
  if (input_data_type == 1 & !is.null(rds_object)) {
    data <- get_data_from_seurat(category, 
                                 rds_object, 
                                 input_file_df, 
                                 assay_name, 
                                 reduction_name, 
                                 assay_data_to_get)
  }
  # if input data type is a SingleCellExperiment object from an RDS file
  else if (input_data_type == 2  & !is.null(rds_object)) {
    data <- get_data_from_sce(category, 
                              rds_object, 
                              input_file_df, 
                              assay_name, 
                              reduction_name, 
                              assay_data_to_get)
  }
  return(as.data.frame(data))
}

