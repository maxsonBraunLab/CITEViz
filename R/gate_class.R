#' Gate class for gate objects
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
#' @return a gate class object
#'
#' @examples
#' 
#' # See examples of Gate function
#' 
methods::setClass("Gate", slots = list(
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
  
  
#' Constructor for Gate objects
#' 
#' All parameters are slots within the gate object
#'
#' @param counter integer to keep track of gate number
#' @param assay_name character. Seurat object assay gates on
#' @param input_cells list. cells plotted
#' @param input_coords dataframe of plotted cell coordinates
#' @param subset_cells list. cells selected
#' @param subset_coords dataframe of selected cell coordinates
#' @param x_axis character. Feature on x-axis
#' @param y_axis character. Feature on y-axis
#' @param gate_coords list. coordinates of drawn gate
#' @param name_subset_cells character. user ipute of selected cell subset name
#' @param num_input_cells integer. number of total cells on plot
#' @param num_subset_cells integer. number of cells selected
#' @param total_num_cells_in_sample integer. number of total cells in seurat object
#' @param pct_subset_from_previous numeric. percentage of cells selected from cells plotted
#' @param pct_subset_from_total numeric. percentage of cells selected from cells in Seurat object
#'
#' @importFrom methods new
#'
#' @return Gate class object
#'
#' @examples
#' # create Gate object
#' example_gate <- Gate(counter = as.integer(1),
#'                    assay_name = "example",
#'                    input_cells = list(c("a", "b", "c")),
#'                    input_coords =  data.frame(x=c(1,2,3),y=c(4,5,6)),
#'                    subset_cells = list(c("a","b")), 
#'                    subset_coords = data.frame(x=c(1,2),y=c(4,5)), 
#'                    x_axis = "ADT-A", 
#'                    y_axis = "ADT-B", 
#'                    gate_coords = list(x=c(1,2,3,4),y=c(5,6,7,8)), 
#'                    name_subset_cells = "example_cells_A", 
#'                    num_input_cells = as.integer(1000), 
#'                    num_subset_cells = as.integer(500), 
#'                    total_num_cells_in_sample = as.integer(1000),
#'                    pct_subset_from_previous = 50, 
#'                    pct_subset_from_total = 50)
Gate <- function(counter=NA_integer_, assay_name = NA_character_, input_cells = list(), input_coords = data.frame(), 
                   subset_cells = list(), subset_coords = data.frame(), x_axis = NA_character_, y_axis = NA_character_,
                   gate_coords = list(), name_subset_cells = NA_character_, num_input_cells = NA_integer_, num_subset_cells = NA_integer_,
                   total_num_cells_in_sample = NA_integer_, pct_subset_from_previous = NA_real_, pct_subset_from_total = NA_real_) {
  methods::new("Gate", 
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
#' @param gate_obj a single Gate class object
#' @param slot_name name of slot to draw data from
#' 
#' @importFrom methods slot
#' 
#' @return data from slot in gate object
#'
#' @examples
#' 
#' # An internal function that operates as such
#' example_gate <- Gate(counter = as.integer(1), 
#'                    assay_name = "example", 
#'                    input_cells = list(c("a", "b", "c")),
#'                    input_coords =  data.frame(x=c(1,2,3),y=c(4,5,6)),
#'                    subset_cells = list(c("a","b")), 
#'                    subset_coords = data.frame(x=c(1,2),y=c(4,5)), 
#'                    x_axis = "ADT-A", 
#'                    y_axis = "ADT-B", 
#'                    gate_coords = list(x=c(1,2,3,4),y=c(5,6,7,8)), 
#'                    name_subset_cells = "example_cells_A", 
#'                    num_input_cells = as.integer(1000), 
#'                    num_subset_cells = as.integer(500), 
#'                    total_num_cells_in_sample = as.integer(1000),
#'                    pct_subset_from_previous = 50, 
#'                    pct_subset_from_total = 50)
#'                    
#' # Retrieve coordinates of gate
#' GetData(example_gate, "gate_coords")
#'
GetData <- function(gate_obj, slot_name) {
  return(methods::slot(gate_obj, slot_name))
}

#' Set Cell Subset Name
#'
#' @param gate_obj a single Gate class object
#' @param new_name character. New name of cell subset
#'
#' @importFrom methods slot
#'
#' @return Gate object with new name of subset cells in name_subset_cells slot
#'
#' @examples
#' 
#' # An internal function to change the name of the cell selection
#' example_gate <- Gate(counter = as.integer(1), 
#'                    assay_name = "example", 
#'                    input_cells = list(c("a", "b", "c")),
#'                    input_coords =  data.frame(x=c(1,2,3),y=c(4,5,6)),
#'                    subset_cells = list(c("a","b")), 
#'                    subset_coords = data.frame(x=c(1,2),y=c(4,5)), 
#'                    x_axis = "ADT-A", 
#'                    y_axis = "ADT-B", 
#'                    gate_coords = list(x=c(1,2,3,4),y=c(5,6,7,8)), 
#'                    name_subset_cells = "example_cells_A", 
#'                    num_input_cells = as.integer(1000), 
#'                    num_subset_cells = as.integer(500), 
#'                    total_num_cells_in_sample = as.integer(1000),
#'                    pct_subset_from_previous = 50, 
#'                    pct_subset_from_total = 50)
#'                    
#' # Change cell subset name name from "example_cells_A" to "custom_name"
#' example_gate <- SetSubsetName(example_gate, "custom_name")
SetSubsetName <- function(gate_obj, new_name) {
  methods::slot(gate_obj, "name_subset_cells") <- new_name
  return(gate_obj)
}