# define a Gate class for gate objects
#' Title
#'
#' @slot counter integer. 
#' @slot assay_name character. 
#' @slot input_cells list. 
#' @slot input_coords data.frame. 
#' @slot subset_cells list. 
#' @slot subset_coords data.frame. 
#' @slot x_axis character. 
#' @slot y_axis character. 
#' @slot gate_coords list. 
#' @slot name_subset_cells character. 
#' @slot num_input_cells integer. 
#' @slot num_subset_cells integer. 
#' @slot total_num_cells_in_sample integer. 
#' @slot pct_subset_from_previous numeric. 
#' @slot pct_subset_from_total numeric. 
#'
#' @return
#' @export
#'
#' @examples
setClass("Gate", slots = list(
    counter = "integer", #to keep track of gate number in the case of multiple gating steps
    assay_name = "character", #to keep track of which assay (RNA, ADT, etc) the gate data is coming from
    
    ## cells and coordinates, seurat data
    input_cells = "list",
    input_coords = "data.frame",
    subset_cells = "list",
    subset_coords = "data.frame",
    
    ## gate information
    x_axis = "character",
    y_axis = "character",
    gate_coords = "list",
    
    ## cell summary data
    name_subset_cells = "character", 
    num_input_cells = "integer",
    num_subset_cells = "integer",
    total_num_cells_in_sample = "integer",
    pct_subset_from_previous = "numeric", #percent of cells in the current gate that were subsetted from previous gate (100 * num_subset_cells/num_input_cells)
    pct_subset_from_total = "numeric" #percent of cells in the current gate that were subsetted from the original total number of cells in the sample
  ))
  
  
#constructor for Gate objects
#' Title
#'
#' @param counter 
#' @param assay_name 
#' @param input_cells 
#' @param input_coords 
#' @param subset_cells 
#' @param subset_coords 
#' @param x_axis 
#' @param y_axis 
#' @param gate_coords 
#' @param name_subset_cells 
#' @param num_input_cells 
#' @param num_subset_cells 
#' @param total_num_cells_in_sample 
#' @param pct_subset_from_previous 
#' @param pct_subset_from_total 
#'
#' @return
#' @export
#'
#' @examples
Gate <- function(counter=NA_integer_, assay_name = NA_character_, input_cells = list(), input_coords = data.frame(), 
                   subset_cells = list(), subset_coords = data.frame(), x_axis = NA_character_, y_axis = NA_character_,
                   gate_coords = list(), name_subset_cells = NA_character_, num_input_cells = NA_integer_, num_subset_cells = NA_integer_,
                   total_num_cells_in_sample = NA_integer_, pct_subset_from_previous = NA_real_, pct_subset_from_total = NA_real_) {
    new("Gate", 
        counter=counter, 
        assay_name = assay_name, 
        input_cells = input_cells, 
        input_coords = input_coords, 
        subset_cells = subset_cells, 
        subset_coords = subset_coords, 
        x_axis = x_axis, 
        y_axis = y_axis,
        gate_coords = gate_coords, 
        name_subset_cells = name_subset_cells, 
        num_input_cells = num_input_cells, 
        num_subset_cells = num_subset_cells,
        total_num_cells_in_sample = total_num_cells_in_sample, 
        pct_subset_from_previous = pct_subset_from_previous, 
        pct_subset_from_total = pct_subset_from_total)
}

#' Slot Accessor Method
#'
#' @param gate_obj 
#' @param slot_name 
#'
#' @return
#' @export
#'
#' @examples
GetData <- function(gate_obj, slot_name) {
  return(slot(gate_obj, slot_name))
}

#' Set Cell Subset Name
#'
#' @param gate_obj 
#' @param new_name 
#'
#' @return
#' @export
#'
#' @examples
SetSubsetName <- function(gate_obj, new_name) {
  slot(gate_obj, "name_subset_cells") <- new_name
  return(gate_obj)
}