#' Gating datatable for display in UI
#'
#' @param temp_reactive_gating_df reactive gating dataframe to use as contents
#'
#' @import magrittr
#' 
#' @return an interactive data table
#'
#' @examples
#' 
#' # This is all done internally within the app
#' 
#' # create an empty gating data frame using create_gating_df()
#' example_df <- create_gating_df()
#' 
#' # create Gate object
#' example_gate1 <- Gate(counter = as.integer(1), 
#'                    assay_name = "example1", 
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
#' # update_gating_df() requires a list of gate objects to choose from 
#' 
#' example_gate2 <- Gate(counter = as.integer(2), 
#'                    assay_name = "example2", 
#'                    input_cells = list(c("a", "b")),
#'                    input_coords =  data.frame(x=c(1,2),y=c(4,5)),
#'                    subset_cells = list(c("a")), 
#'                    subset_coords = data.frame(x=c(1),y=c(4)), 
#'                    x_axis = "ADT-C", 
#'                    y_axis = "ADT-D", 
#'                    gate_coords = list(x=c(10,20,30,40),y=c(50,60,70,80)), 
#'                    name_subset_cells = "example_cells_B", 
#'                    num_input_cells = as.integer(500), 
#'                    num_subset_cells = as.integer(100), 
#'                    total_num_cells_in_sample = as.integer(1000),
#'                    pct_subset_from_previous = 20, 
#'                    pct_subset_from_total = 10)
#' 
#' example_gate_list <- list("gate1" = example_gate1, "gate2" = example_gate2)
#' 
#' example_df <- update_gating_df(gate_name_string = "gate1", 
#'     reactive_gate_list = example_gate_list, temp_gating_df = example_df)
#' example_df <- update_gating_df(gate_name_string = "gate2", 
#'     reactive_gate_list = example_gate_list, temp_gating_df = example_df)
#' 
#' create_gating_dt(example_df)
create_gating_dt <- function(temp_reactive_gating_df) {
    temp_gating_dt <- DT::datatable(temp_reactive_gating_df, 
                                    rownames = TRUE, 
                                    selection = "single", #only let user click 1 row at a time so only that one gate can be shown in plot. Otherwise, app will crash.
                                    editable = list(target = "cell", 
                                                    disable = list(columns = c(0:4, 6:ncol(temp_reactive_gating_df)))), # only cells in column 5(subset_name) can be edited by user
                                    extensions = c("Buttons", "Scroller"),
                                    options = list(deferRender = TRUE,
                                                scroller = TRUE,
                                                scrollY = 300,
                                                scrollX = TRUE,
                                                dom = "frtipB", 
                                                buttons = c("copy", "print"),
                                                columnDefs = list(list(visible = FALSE, targets = c(11:ncol(temp_reactive_gating_df))))
                                    )) %>%
    DT::formatRound(columns = c("Percent_Subsetted_From_Previous", 
                                "Percent_Subsetted_From_Total"), 
                    digits = 4)
    
    return(temp_gating_dt)
}