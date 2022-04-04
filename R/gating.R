#create gate objects from user input selections depending on if the user is forward-gating or back-gating
#' Title
#'
#' @param input
#' @param is_forward_gating 
#' @param assay_count_data 
#' @param gate_counter 
#' @param reactive_gate_list 
#' @param reactive_selected_gate 
#' @param reactive_last_buttons_clicked 
#'
#' @return
#' @export
#'
#' @examples
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

    #assign values to variables
    assay_name <- input[[paste0("Assay", ui_input_suffix)]]
    subset_cells <- list(sel$customdata)
    subset_coords <- data.frame(x = assay_count_data[sel$customdata,1], y = assay_count_data[sel$customdata,2])
    x_axis <- input[[paste0("x_feature", ui_input_suffix)]]
    y_axis <- input[[paste0("y_feature", ui_input_suffix)]]
    gate_coords <- list(x = c(brushed_coords$x), y = c(brushed_coords$y))
    name_subset_cells <- input[[paste0("cell_subset_name", ui_input_suffix)]]
    num_input_cells <- length(unlist(input_cells))
    num_subset_cells <- length(unlist(subset_cells))
    total_num_cells_in_sample <- length(rownames(assay_count_data))
    pct_subset_from_previous <- 100 * num_subset_cells / num_input_cells
    pct_subset_from_total <- 100 * num_subset_cells / total_num_cells_in_sample

    # construct Gate object
    gate_obj <- Gate(counter = gate_counter,
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
    return(gate_obj)
  }
  
  
#' Title
#'
#' @param gating_reactiveValues 
#'
#' @return
#' @export
#'
#' @examples
get_reactive_gate_list <- function(gating_reactiveValues) {
    reactive_gate_list <- reactive({
      unordered_list <- reactiveValuesToList(gating_reactiveValues)
      ordered_list <- unordered_list[order(names(unordered_list))]
      ordered_list
    })
    return(reactive_gate_list)
  }
  
#create gating dataframe
#' Title
#'
#' @return
#' @export
#'
#' @examples
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
  
#function to populate/update gating data frame
#' Title
#'
#' @param gate_name_string 
#' @param reactive_gate_list 
#' @param temp_gating_df 
#'
#' @return
#' @export
#'
#' @examples
update_gating_df <- function(gate_name_string, reactive_gate_list, temp_gating_df) {
    gate_obj <- reactive_gate_list[[gate_name_string]]
    if (!is.null(gate_obj)) {
      temp_gating_df <- add_row(temp_gating_df,
                                Gate_ID = gate_name_string,
                                Assay_Name = gate_obj@assay_name,
                                X_Axis = gate_obj@x_axis,
                                Y_Axis = gate_obj@y_axis,
                                Subset_Name = gate_obj@name_subset_cells,
                                Num_Input_Cells = gate_obj@num_input_cells,
                                Num_Subset_Cells = gate_obj@num_subset_cells,
                                
                                Total_Cells_in_Sample = gate_obj@total_num_cells_in_sample,
                                Percent_Subsetted_From_Previous = gate_obj@pct_subset_from_previous,
                                Percent_Subsetted_From_Total = gate_obj@pct_subset_from_total,
                                
                                Input_Cells = gate_obj@input_cells,
                                Input_X_Coordinates = list(gate_obj@input_coords$x),
                                Input_Y_Coordinates = list(gate_obj@input_coords$y),
                                Subset_Cells = gate_obj@subset_cells,
                                Subset_X_Coordinates = list(gate_obj@subset_coords$x),
                                Subset_Y_Coordinates = list(gate_obj@subset_coords$y),
                                Gate_X_Coordinates = list(gate_obj@gate_coords$x),
                                Gate_Y_Coordinates = list(gate_obj@gate_coords$y)
      )
    }
    return(temp_gating_df)
  }
  
  
# reset all values in gate_reactive_values to NULL
#' Title
#'
#' @param gate_name_string 
#' @param local_gate_reactive_values 
#'
#' @return
#' @export
#'
#' @examples
set_gates_to_null <- function(gate_name_string, local_gate_reactive_values) {
    local_gate_reactive_values[[gate_name_string]] <- NULL
    return(local_gate_reactive_values[[gate_name_string]])
  }