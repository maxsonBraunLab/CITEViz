#' Create Feature Expression Plot
#'
#'
#' @param input Shiny internal parameter object containing UI user input values
#' @param input_data_type An integer value indicating if the user-uploaded input files are RDS files holding Seurat objects or SingleCellExperiment objects (1 = Seurat object, 2 = SingleCellExperiment object, etc.)
#' @param rds_object An RDS object containing metadata, assays, and reductions for a CITE-seq experiment. Can be NULL if the user uploaded an alternate supported file instead of an RDS file.
#'
#' @import magrittr
#' @importFrom ggplot2 aes ggplot geom_tile labs scale_x_continuous scale_y_continuous
#' @importFrom shiny req
#' @importFrom dplyr left_join
#' @importFrom plotly plot_ly config layout event_register
#' @return dimension reduction plot with the correct colors for a 2-feature co-expression plot
#' @noRd
expression_plot <- function(input, input_data_type, rds_object) {

    req(
        input$file_input, input$reduction_expr_1d,
        input$Assay_1d, input$feature_1d
    )

    # create string for reduction to plot
    reduc <- input$reduction_expr_1d

    # selected feature to color clusters by
    color_x <- input$feature_1d

    count_data <- get_data(
        category = "assays",
        input_data_type = input_data_type,
        rds_object = rds_object,
        input_file_df = input_file_df,
        assay_name = input$Assay_1d,
        reduction_name = NULL,
        assay_data_to_get = color_x
    )

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

    p <- plot_ly(cell_data,
        x = ~ cell_data[, 1],
        y = ~ cell_data[, 2],
        customdata = rownames(cell_data),
        type = "scatter",
        mode = "markers",
        marker = list(
            size = 3,
            color = ~ count_data[, color_x],
            colorbar = list(
                title = color_x,
                len = 0.5
            ),
            colorscale = "Viridis",
            reversescale = TRUE
        ),
        source = "expression_1d_plot"
    ) %>%
    plotly::config(
        toImageButtonOptions = list(
            format = "png",
            scale = 10
        )
    ) %>%
    plotly::layout(
        title = toupper(reduc),
        xaxis = list(title = cell_col[1]),
        yaxis = list(title = cell_col[2]),
        dragmode = "select"
    ) %>%
    event_register("plotly_selected")
    
    return(p)

}