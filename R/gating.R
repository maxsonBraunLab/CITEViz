#' Create gate objects from user input selections
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
            input_coords <- data.frame(x = assay_count_data[, 1], y = assay_count_data[, 2])
        } else {
            # set up inputs to gate object if there is already a gate object
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
        input_coords <- data.frame(x = assay_count_data[, 1], y = assay_count_data[, 2])
    }
    # assign values to variables that will be referenced multiple times in Gate object creation
    subset_cells <- list(sel$customdata)
    num_input_cells <- length(unlist(input_cells))
    num_subset_cells <- length(unlist(subset_cells))
    total_num_cells_in_sample <- length(rownames(assay_count_data))

    # construct Gate object
    gate_obj <- Gate(
        counter = gate_counter,
        assay_name = input[[paste0("Assay", ui_input_suffix)]],
        input_cells = input_cells,
        input_coords = input_coords,
        subset_cells = subset_cells,
        subset_coords = data.frame(
            x = assay_count_data[sel$customdata, 1],
            y = assay_count_data[sel$customdata, 2]
        ),
        x_axis = input[[paste0("x_feature", ui_input_suffix)]],
        y_axis = input[[paste0("y_feature", ui_input_suffix)]],
        gate_coords = list(x = c(brushed_coords$x), y = c(brushed_coords$y)),
        name_subset_cells = input[[paste0("cell_subset_name", ui_input_suffix)]],
        num_input_cells = num_input_cells,
        num_subset_cells = num_subset_cells,
        total_num_cells_in_sample = total_num_cells_in_sample,
        pct_subset_from_previous = 100 * num_subset_cells / num_input_cells,
        pct_subset_from_total = 100 * num_subset_cells / total_num_cells_in_sample
    )
    return(gate_obj)
}


#' Convert Shiny reactiveValues object to reactive gate list
#'
#' @param gating_reactive_values Shiny reactiveValues object holding all reactive gate objects
#'
#' @return reactive gate list for app purposes
#'
#' @noRd
#'
get_reactive_gate_list <- function(gating_reactive_values) {
    reactive_gate_list <- reactive({
        unordered_list <- reactiveValuesToList(gating_reactive_values)
        ordered_list <- unordered_list[order(names(unordered_list))]
        ordered_list
    })
    return(reactive_gate_list)
}

#' Collect all gating information into a dataframe
#' @noRd
#' @return Empty gate dataframe for app purposes
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
#' @noRd
#' @return updated gating dataframe
#'
update_gating_df <- function(
    gate_name_string,
    reactive_gate_list,
    temp_gating_df
) {
    gate_obj <- reactive_gate_list[[gate_name_string]]
    if (!is.null(gate_obj)) {
        temp_gating_df <- tibble::add_row(temp_gating_df,
            Gate_ID = gate_name_string,
            Assay_Name = get_gate_data(gate_obj, "assay_name"),
            X_Axis = get_gate_data(gate_obj, "x_axis"),
            Y_Axis = get_gate_data(gate_obj, "y_axis"),
            Subset_Name = get_gate_data(gate_obj, "name_subset_cells"),
            Num_Input_Cells = get_gate_data(gate_obj, "num_input_cells"),
            Num_Subset_Cells = get_gate_data(gate_obj, "num_subset_cells"),
            Total_Cells_in_Sample = get_gate_data(gate_obj, "total_num_cells_in_sample"),
            Percent_Subsetted_From_Previous = get_gate_data(gate_obj, "pct_subset_from_previous"),
            Percent_Subsetted_From_Total = get_gate_data(gate_obj, "pct_subset_from_total"),
            Input_Cells = get_gate_data(gate_obj, "input_cells"),
            Input_X_Coordinates = list(get_gate_data(gate_obj, "input_coords")$x),
            Input_Y_Coordinates = list(get_gate_data(gate_obj, "input_coords")$y),
            Subset_Cells = get_gate_data(gate_obj, "subset_cells"),
            Subset_X_Coordinates = list(get_gate_data(gate_obj, "subset_coords")$x),
            Subset_Y_Coordinates = list(get_gate_data(gate_obj, "subset_coords")$y),
            Gate_X_Coordinates = list(get_gate_data(gate_obj, "gate_coords")$x),
            Gate_Y_Coordinates = list(get_gate_data(gate_obj, "gate_coords")$y)
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

#' Normal Gate Scatterplot
#' 
#' @param input
#' @param input_data_type
#' @param rds_object
#' @param last_buttons_clicked
#' @param gate_list
#' 
#' @importFrom plotly plot_ly event_register layout add_histogram2dcontour add_markers config event_data
#' 
#' @return scatterplot for normal gating
#' 
#' @noRd
gate_scatterplot <- function(
    input,
    input_data_type,
    rds_object,
    last_buttons_clicked,
    gate_list,
    selected_gate
) {

    # Give some documentation here on the last buttons clicked method.

    # code to execute when one of the above input events occurs
    req(
        input$x_feature,
        input$y_feature
    )

    count_data <- get_data(
        category = "assays",
        input_data_type = input_data_type,
        rds_object = rds_object,
        input_file_df = input_file_df,
        assay_name = input$Assay,
        reduction_name = NULL,
        assay_data_to_get = c(input$x_feature, input$y_feature)
    )

    # generate dataframe for custom colorscale for contour plot, where each hex color code is mapped to a specific z-value between 0 and 1 (inclusive)
    # colorscale needs to be in this format for Plotly's add_histogram2dcontour(colorscale = ...) parameter
    gating_color_scale <- data.frame(
        z = c(0.0, 0.20, 0.40, 0.60, 0.80, 1.0),
        col = c("#FFFFFF", "#4564FE", "#76EFFF", "#FFF900", "#FFA300", "#FF1818")
    )

    # initialize a base_scatterplot variable before if/else statements below so that the plot object can be accessed outside of the if/else statements
    base_scatterplot <- NULL

    # generate a base scatterplot based on user input
    # needs some more documentation on its behavior here
    if ((last_buttons_clicked$last == "NA" | last_buttons_clicked$last == "reset_button" | last_buttons_clicked$last == "clear_all_gates_button") & is.null(input$gating_pg_table_rows_selected)) {

        base_scatterplot <- plot_ly(
            data = count_data,
            x = ~ count_data[, input$x_feature],
            y = ~ count_data[, input$y_feature],
            customdata = rownames(count_data),
            mode = "markers",
            source = "C") %>%
        add_histogram2dcontour(
            showscale = FALSE,
            ncontours = 10,
            colorscale = gating_color_scale,
            contours = list(coloring = "heatmap")) %>%
        add_markers(
            x = count_data[, input$x_feature],
            y = count_data[, input$y_feature],
            marker = list(size = 2),
            color = I("black"),
            alpha = 0.6
        )

    } else {

        selected_cell_barcodes <- NULL
        count_data_subset <- NULL

        if (last_buttons_clicked$last == "gate_button" & is.null(input$gating_pg_table_rows_selected)) {
            selected_cell_barcodes <- event_data("plotly_selected", source = "C")$customdata
            count_data_subset <- count_data[rownames(count_data) %in% selected_cell_barcodes, ]

        } else if (!is.null(input$gating_pg_table_rows_selected)) {
            selected_cell_barcodes <- get_gate_data(gate_list[[selected_gate]], "subset_cells")[[1]]
            count_data_subset <- count_data[rownames(count_data) %in% selected_cell_barcodes, ]
        }

        base_scatterplot <- plot_ly(
            data = count_data_subset,
            x = ~ count_data_subset[, input$x_feature],
            y = ~ count_data_subset[, input$y_feature],
            customdata = rownames(count_data_subset),
            mode = "markers",
            source = "C") %>%
        add_histogram2dcontour(
            showscale = FALSE, ncontours = 10, colorscale = gating_color_scale,
            contours = list(coloring = "heatmap")) %>%
        add_markers(
            x = count_data_subset[, input$x_feature],
            y = count_data_subset[, input$y_feature],
            marker = list(size = 2.5),
            color = I("black"),
            alpha = 0.6
        )
    }

    # add configuration and layout options to base scatterplot, register selection events
    p <- base_scatterplot %>%
        config(
            toImageButtonOptions = list(
                format = "png",
                scale = 10
            ) # scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
        ) %>%
        # Layout changes the aesthetic of the plot
        layout(
            title = "Normalized Feature Scatter Plot",
            xaxis = list(title = input$x_feature),
            yaxis = list(title = input$y_feature),
            showlegend = FALSE,
            dragmode = "select"
        ) %>% # Determines the mode of drag interactions. "select" and "lasso" apply only to scatter traces with markers or text. "orbit" and "turntable" apply only to 3D scenes.
        event_register("plotly_selected")

    return(p)

}

#' Normal Gate Scatterplot
#' 
#' @param input
#' @param input_data_type
#' @param rds_object
#' @param last_buttons_clicked
#' @param gate_list
#' @param selected_gate
#' 
#' @importFrom plotly plot_ly event_register layout add_histogram2dcontour
#' @importFrom plotly add_markers config event_data
#' 
#' @return scatterplot for normal gating
#' 
#' @noRd
gate_reduction <- function(
    input,
    input_data_type,
    rds_object,
    last_buttons_clicked,
    gate_list,
    selected_gate
) {

    req(
        input$file_input,
        input$gating_color_dimred,
        input$gating_reduction,
    )

    # create string for reduction to plot
    reduc <- input$gating_reduction

    # selected metadata to color clusters by
    color <- input$gating_color_dimred

    # gather metadata
    metadata_df <- get_data(
        category = "metadata",
        input_data_type = input_data_type,
        rds_object = rds_object,
        input_file_df = input_file_df,
        assay_name = NULL,
        reduction_name = NULL
    )

    # interpolate the base color palette
    custom_palette <- get_palette(length(unique(metadata_df[[color]])))

    # create dataframe from reduction selected
    cell_data <- get_data(
        category = "reductions",
        input_data_type = input_data_type,
        rds_object = rds_object,
        input_file_df = input_file_df,
        assay_name = NULL,
        reduction_name = reduc
    )

    # create list containing all column names of cell_data
    cell_col <- colnames(cell_data)

    # initialize selected_cells
    selected_cells <- NULL

    # define selected cells by either user
    # input or based on a previous gate accessed by a click
    if (!is.null(input$gating_pg_table_rows_selected)) {
        selected_cells <- get_gate_data(gate_list[[selected_gate]], "subset_cells")[[1]]
    } else {
        selected_cells <- event_data("plotly_selected", source = "C")$customdata
    }

    # color the reduction plot by user metadata (prior to any gates)
    # or black-and-white to highlight selected cells after gates are made
    if (is.null(selected_cells)) {
        plotly_color_list <- c(paste0("metadata_df$", color), "custom_palette")
    } else {
        plotly_color_list <- c("rownames(cell_data) %in% selected_cells", 'c("grey", "black")')
    }

    # plot gating dimension reduction plot
    p <- plot_ly(
        data = cell_data,
        x = ~ cell_data[, 1],
        y = ~ cell_data[, 2],
        customdata = rownames(cell_data),
        color = stats::as.formula(paste0("~", plotly_color_list[1])), # color by selected metadata in object
        colors = stats::as.formula(paste0("~", plotly_color_list[2])),
        type = "scatter",
        mode = "markers",
        marker = list(size = 3, width = 2)) %>%
    config(
        toImageButtonOptions = list(
            format = "png",
            scale = 10)) %>%
    layout(
        title = toupper(reduc),
        xaxis = list(title = cell_col[1]),
        yaxis = list(title = cell_col[2]),
        legend = list(itemsizing = "constant")
    )

    return(p)

}