#' #create gate objects from user input selections
#'
#' @param input Reactive container for ui inputs
#' @param is_forward_gating Boolean for if the gating is forward or back gating
#' @param assay_count_data Feature expression data
#' @param gate_counter Integer of gate number
#' @param reactive_gate_list Reactive list of gates
#' @param reactive_selected_gate Reactive value of user chosen gate object
#' @param reactive_last_buttons_clicked Reactive list object tracking the last two buttons clicked of three possibilities: "gate_button","reset_button","clear_all_gates_button"
#'
#' @return a Gate object
#' @noRd
#' 
create_gate_from_input <- function(input, is_forward_gating = TRUE, assay_count_data, gate_counter, reactive_gate_list, reactive_selected_gate, reactive_last_buttons_clicked) {
    # "ui_input_suffix" refers to the suffix that is at the end of the UI input elements for forward and backgating
    # So for forward-gating, an example of an input from a UI element would be input$Assay. The corresponding input in the backgating page would be input$Assay_bg. the suffix for the backgating UI elements is "_bg"
    ui_input_suffix <- ""
    sel <- NA_character_
    brushed_coords <- NA_character_
    input_cells <- list()
    input_coords <- data.frame()

    # forward-gating logic
    if (is_forward_gating == TRUE) {
      # get plotly event data
      sel <- event_data("plotly_selected", source = "C")
      brushed_coords <- event_data("plotly_brushed", source = "C")
      # Check for Gate object. If one does not exists, create initial input for first gate.
      if (length(reactive_gate_list) == 0 | reactive_last_buttons_clicked$second_to_last == "reset_button" | reactive_last_buttons_clicked$second_to_last == "clear_all_gates_button") {
        input_cells <- list(rownames(assay_count_data))
        input_coords <- data.frame(x = assay_count_data[,1], y = assay_count_data[,2])
      }
      else {
        #set up inputs to gate object if there is already a gate object
        # if user has not selected a previous gate from datatable to re-gate from
        if (is.null(input$gating_pg_table_rows_selected)) {
          input_cells <- reactive_gate_list[[paste0("gate_", gate_counter - 1)]]@subset_cells
          input_coords <- data.frame(x = brushed_coords$x, y = brushed_coords$y)
        }
        # else if user has selected a previous gate from datatable that they want to re-gate from
        else {
          input_cells <- reactive_gate_list[[reactive_selected_gate]]@subset_cells
          input_coords <- reactive_gate_list[[reactive_selected_gate]]@subset_coords
        }
      }
    }
    # back-gating logic
    else {
      ui_input_suffix <- "_bg"
      # get plotly event data and input cell data
      sel <- event_data("plotly_selected", source = "D")
      brushed_coords <- event_data("plotly_brushed", source = "D")
      input_cells <- list(rownames(assay_count_data))
      input_coords <- data.frame(x = assay_count_data[,1], y = assay_count_data[,2])
    }
    #assign values to variables that will be referenced multiple times in Gate object creation
    subset_cells <- list(sel$customdata)
    num_input_cells <- length(unlist(input_cells))
    num_subset_cells <- length(unlist(subset_cells))
    total_num_cells_in_sample <- length(rownames(assay_count_data))

    # construct Gate object
    gate_obj <- Gate(counter = gate_counter,
                    assay_name = input[[paste0("Assay", ui_input_suffix)]],
                    input_cells = input_cells,
                    input_coords = input_coords,
                    subset_cells = subset_cells,
                    subset_coords = data.frame(x = assay_count_data[sel$customdata,1], 
                                               y = assay_count_data[sel$customdata,2]),
                    x_axis = input[[paste0("x_feature", ui_input_suffix)]],
                    y_axis = input[[paste0("y_feature", ui_input_suffix)]],
                    gate_coords = list(x = c(brushed_coords$x), y = c(brushed_coords$y)),
                    name_subset_cells = input[[paste0("cell_subset_name", ui_input_suffix)]],
                    num_input_cells = num_input_cells,
                    num_subset_cells = num_subset_cells,
                    total_num_cells_in_sample = total_num_cells_in_sample,
                    pct_subset_from_previous = 100 * num_subset_cells / num_input_cells,
                    pct_subset_from_total = 100 * num_subset_cells / total_num_cells_in_sample)
    return(gate_obj)
  }
  
  
#' Convert Shiny reactiveValues object to reactive gate list
#'
#' @param gating_reactiveValues Shiny reactiveValues object holding all reactive gate objects
#'
#' @return reactive gate list for app purposes
#' @export
#'
#' @noRd
#' 
get_reactive_gate_list <- function(gating_reactiveValues) {
    reactive_gate_list <- reactive({
      unordered_list <- reactiveValuesToList(gating_reactiveValues)
      ordered_list <- unordered_list[order(names(unordered_list))]
      ordered_list
    })
    return(reactive_gate_list)
  }
  
#' Create gating dataframe
#' 
#' @return Empty gate dataframe for app purposes
#' @export
#' 
#' @examples
#' 
#' # Takes no arguments
#' example_df <- create_gating_df()
#' 
create_gating_df <- function() {
    data.frame(
      Gate_ID = character(),
      Assay_Name = character(),
      X_Axis = character(),
      Y_Axis = character(),
      Subset_Name = character(),
      Num_Input_Cells = integer(),
      Num_Subset_Cells = integer(),
      Total_Cells_in_Sample = integer(),
      Percent_Subsetted_From_Previous = numeric(),
      Percent_Subsetted_From_Total = numeric(),
      
      Input_Cells = I(list()),
      Input_X_Coordinates = I(list()),
      Input_Y_Coordinates = I(list()),
      Subset_Cells = I(list()),
      Subset_X_Coordinates = I(list()),
      Subset_Y_Coordinates = I(list()),
      Gate_X_Coordinates = I(list()),
      Gate_Y_Coordinates = I(list())
    )
  }
  
#' Populate/update gating data frame
#'
#' @param gate_name_string character string of new gate name
#' @param reactive_gate_list shiny reactivelist of gates
#' @param temp_gating_df empty gating dataframe built from create_gating_df function
#'
#' @importFrom tibble add_row
#' @return updated gating dataframe 
#' @export
#'
#' @examples
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
#' # Inputs contents of Gate object into empty gating dataframe
#' example_df <- update_gating_df(gate_name_string = "gate1", 
#'     reactive_gate_list = example_gate_list, temp_gating_df = example_df)
#' example_df <- update_gating_df(gate_name_string = "gate2", 
#'     reactive_gate_list = example_gate_list, temp_gating_df = example_df)
#' 
update_gating_df <- function(gate_name_string, reactive_gate_list, temp_gating_df) {
    gate_obj <- reactive_gate_list[[gate_name_string]]
    if (!is.null(gate_obj)) {
      temp_gating_df <- tibble::add_row(temp_gating_df,
                                Gate_ID = gate_name_string,
                                Assay_Name = GetData(gate_obj,"assay_name"),
                                X_Axis = GetData(gate_obj,"x_axis"),
                                Y_Axis = GetData(gate_obj,"y_axis"),
                                Subset_Name = GetData(gate_obj,"name_subset_cells"),
                                Num_Input_Cells = GetData(gate_obj,"num_input_cells"),
                                Num_Subset_Cells = GetData(gate_obj,"num_subset_cells"),
                                
                                Total_Cells_in_Sample = GetData(gate_obj,"total_num_cells_in_sample"),
                                Percent_Subsetted_From_Previous = GetData(gate_obj,"pct_subset_from_previous"),
                                Percent_Subsetted_From_Total = GetData(gate_obj,"pct_subset_from_total"),
                                
                                Input_Cells = GetData(gate_obj,"input_cells"),
                                Input_X_Coordinates = list(GetData(gate_obj,"input_coords")$x),
                                Input_Y_Coordinates = list(GetData(gate_obj,"input_coords")$y),
                                Subset_Cells = GetData(gate_obj,"subset_cells"),
                                Subset_X_Coordinates = list(GetData(gate_obj,"subset_coords")$x),
                                Subset_Y_Coordinates = list(GetData(gate_obj,"subset_coords")$y),
                                Gate_X_Coordinates = list(GetData(gate_obj,"gate_coords")$x),
                                Gate_Y_Coordinates = list(GetData(gate_obj,"gate_coords")$y)
      )
    }
    return(temp_gating_df)
  }
  
  
#' Reset all values in gate_reactive_values to NULL
#'
#' @param gate_name_string name of gate to be deleted
#' @param local_gate_reactive_values Shiny reactiveValues
#'
#' @return empty reactiveValues list of Gate objects
#' @noRd
#' 
set_gates_to_null <- function(gate_name_string, local_gate_reactive_values) {
    local_gate_reactive_values[[gate_name_string]] <- NULL
    return(local_gate_reactive_values[[gate_name_string]])
  }