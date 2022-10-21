#' Gate class for gate objects
#'
#' @description The Gate class is the core data structure of CITEViz to store
#' gating information. Gates are initiated upon the first click of the "Gate"
#' button and turned into a reactive value to facilitate interactivity with 
#' R-Shiny. The Gate class holds the following information:
#' 
#' * counter <integer>: keep track of gate number
#' * assay_name <character>: keep track of which assay (RNA, ADT, etc) the gate data is coming from
#' * input_cells <list>: a list of cell barcodes going into the gate 
#' e.g. if a gate isn't initialized yet, this variable contains all the cells.
#' * input_coords <data.frame>: the coordinates of the input cells
#' * subset_cells <list>: a list of output cells
#' * subset_coords <data.frame>: the coordinates of the output cells
#' * x_axis <character>: the feature used to plot the x-axis
#' * y_axis <character>: the feature used to plot the y-axis
#' * gate_coords <list>: the coordinates of the gate boundaries drawn by plotly
#' * name_subset_cells <character>: the label of the set of output cells inputted by the user
#' * num_input_cells <integer>: number of input cells
#' * num_subset_cells <integer>: number of output cells
#' * total_num_cells_in_sample <integer>: total number of cell in a sample
#' * pct_subset_from_previous <numeric>: proportion of output cell WRT previous gate
#' * pct_subset_from_total <numeric>: proportion of output cell WRT total cells
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
#' @export
#'
methods::setClass("Gate", slots = list(
    counter = "integer",
    assay_name = "character",

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
    pct_subset_from_previous = "numeric",
    pct_subset_from_total = "numeric" # percent of cells in the current gate that were subsetted from the original total number of cells in the sample
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
#' @param name_subset_cells character. user input of selected cell subset name
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
Gate <- function(counter = NA_integer_, assay_name = NA_character_, input_cells = list(), input_coords = data.frame(),
    subset_cells = list(), subset_coords = data.frame(), x_axis = NA_character_, y_axis = NA_character_,
    gate_coords = list(), name_subset_cells = NA_character_, num_input_cells = NA_integer_, num_subset_cells = NA_integer_,
    total_num_cells_in_sample = NA_integer_, pct_subset_from_previous = NA_real_, pct_subset_from_total = NA_real_) {
    methods::new("Gate",
        counter = counter,
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
        pct_subset_from_total = pct_subset_from_total
    )
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
SetSubsetName <- function(gate_obj, new_name) {
    methods::slot(gate_obj, "name_subset_cells") <- new_name
    return(gate_obj)
}
