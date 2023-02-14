#' Quality Control: Distribution Plot
#'
#' Visualize common QC metrics in CITE-Seq using density plots.
#'
#' @importFrom shiny req
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes labs theme ggtitle xlab ylab geom_vline scale_color_manual scale_fill_manual scale_x_log10 geom_bar geom_density
#' @importFrom ggplotly config layout
#' @noRd
#'
#' @return numeric vector for color values used in coexpression plot
qc_dist_plot <- function(input, input_data_type, rds_object) {

    # required info
    req(
        input$file_input,
        input$QC,
        input$color_qc
    )

    color <- input$color_qc

    # This assigns the params variable a list of strings that act a varying parameters depending on input of QC
    params <- switch(input$QC,
        "RNA Count Per Cell" = c("nCount_RNA", "Distribution of Counts per Cell", "Number of Counts"),
        "Gene Count Per Cell" = c("nFeature_RNA", "Distribution of Genes Detected per Cell", "Number of Unique Genes"),
        "ADT Count Per Cell" = c("nCount_ADT", "ADT Counts per Cell", "Number of Counts"),
        "Unique ADTs Per Cell" = c("nFeature_ADT", "Distribution of CITE-seq Antibodies per Cell", "Number of Unique Antibodies")
    )

    metadata_df <- get_data(
        category = "metadata",
        input_data_type = input_data_type,
        rds_object = rds_object,
        input_file_df = input_file_df,
        assay_name = NULL,
        reduction_name = NULL
    )

    quant <- quantile(
        x = metadata_df[, params[1]],
        probs = c(0.5, 0.75, 0.95),
        na.rm = TRUE
    )

    # interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
    custom_palette <- get_palette(length(unique(metadata_df[[color]])))

    # generate base plot template with features that all QC distribution plots will have
    base_distrib_plot <- metadata_df %>%
        ggplot(aes(x = !!as.name(params[1]), fill = !!as.name(color), color = !!as.name(color))) +
        labs(fill = color, color = color) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(params[2]) +
        xlab(params[3]) +
        geom_vline(xintercept = quant, size = 0.5, alpha = 0.5, linetype = "dashed", color = "grey30") +
        scale_color_manual(values = custom_palette) +
        scale_fill_manual(values = custom_palette)

    # initialize QC distribution plot before if/else statements below so that the plot object can be accessed outside of the if/else statements
    final_distrib_plot <- base_distrib_plot

    # create density/bar plot for selected input. If integrated object is uploaded, then the original identity of the cells will separate into graphs per sample
    if (input$QC %in% "ADT Count Per Cell") {
        final_distrib_plot <- base_distrib_plot + scale_x_log10() + geom_density(alpha = 0.25) + ylab("Density")
    } else if (input$QC %in% "Unique ADTs Per Cell") {
        final_distrib_plot <- base_distrib_plot + geom_bar(alpha = 0.5, position = "dodge") + ylab("Frequency")
    } else {
        final_distrib_plot <- base_distrib_plot + geom_density(alpha = 0.25) + ylab("Density")
    }

    # show distribution plot
    final_distrib_plot <- final_distrib_plot %>%
        ggplotly() %>%
        config(
            toImageButtonOptions = list(format = "png", scale = 10) # scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
        ) %>%
        layout(title = list(font = list(size = 14)), hovermode = FALSE)

    return(final_distrib_plot)
}

#' Quality Control: Box Plot
#'
#' Visualize common QC metrics in CITE-Seq using box plots.
#'
#' @importFrom shiny req
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes labs theme ggtitle xlab ylab geom_vline scale_color_manual scale_fill_manual scale_x_log10 scale_y_log10 geom_bar geom_density
#' @importFrom ggplotly config layout
#' @noRd
#'
#' @return numeric vector for color values used in coexpression plot
qc_box_plot <- function(input, input_data_type, rds_object) {

    # required info
    req(
        input$file_input,
        input$QC,
        input$color_qc
    )

    color <- input$color_qc

    # This assigns the params variable a list of strings that act a varying parameters depending on input of QC
    params <- switch(input$QC,
        "RNA Count Per Cell" = c("nCount_RNA", "Distribution of Counts per Cell", "Number of Counts"),
        "Gene Count Per Cell" = c("nFeature_RNA", "Distribution of Genes Detected per Cell", "Number of Unique Genes"),
        "ADT Count Per Cell" = c("nCount_ADT", "ADT Counts per Cell", "Number of Counts"),
        "Unique ADTs Per Cell" = c("nFeature_ADT", "Distribution of CITE-seq Antibodies per Cell", "Number of Unique Antibodies")
    )

    metadata_df <- get_data(
        category = "metadata",
        input_data_type = input_data_type,
        rds_object = rds_object,
        input_file_df = input_file_df,
        assay_name = NULL,
        reduction_name = NULL
    )

    quant <- quantile(
        x = metadata_df[, params[1]],
        probs = c(0.5, 0.75, 0.95),
        na.rm = TRUE
    )

    # interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
    custom_palette <- get_palette(length(unique(metadata_df[[color]])))

    # generate base plot template with features that all QC distribution plots will have
    base_box_plot <- metadata_df %>%
        ggplot(aes(x = !!as.name(color), y = !!as.name(params[1]), fill = !!as.name(color), color = !!as.name(color))) +
        labs(fill = color, color = color) +
        geom_boxplot(alpha = 0.5, width = 0.5) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(params[2]) +
        xlab("Sample") +
        ylab(params[3]) +
        geom_violin(alpha = 0.2) +
        geom_hline(yintercept = quant, size = 0.5, alpha = 0.5, linetype = "dashed", color = "grey30") +
        scale_color_manual(values = custom_palette) +
        scale_fill_manual(values = custom_palette)

    # initialize QC box plot before if/else statements below so that the plot object can be accessed outside of the if/else statements
    final_box_plot <- base_box_plot

    # create box plot for selected input. If integrated object is uploaded, then the original identity of the cells will separate into boxes per sample
    if (input$QC %in% "ADT Count Per Cell") {
        final_box_plot <- base_box_plot + scale_y_log10()
    }

    # show box plot
    final_box_plot <- final_box_plot %>%
        ggplotly() %>%
        config(
            toImageButtonOptions = list(
                format = "png",
                scale = 10
            ) # scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
        ) %>%
        layout(title = list(font = list(size = 14)), hovermode = FALSE)

    return(final_box_plot)
}