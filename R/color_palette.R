#' Custom color palette with color interpolation to use for QA plots
#'
#' create custom color palette with color interpolation to use for QA plots
#' interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
#' takes in an integer argument for number of colors to generate in palette (default is 3)
#' returns a character vector of color hex codes for the desired number of colors
#' 
#' @param num_colors number of colors
#'
#' @importFrom grDevices colorRampPalette
#'
#' @return A color palette in the form of a character vector of hex color codes.
#' @export
#'
#' @examples 
#' # passing an integer directly into function
#' get_palette(8)
#' 
#' # passing a variable holding an integer value into function
#' number_of_colors <- 10
#' get_palette(number_of_colors)
get_palette <- function(num_colors = 3) {
    base_palette <- grDevices::colorRampPalette(colors = c("turquoise3",
                                                "darkgreen",
                                                "red",
                                                "royalblue1",
                                                "orange2",
                                                "lightseagreen",
                                                "yellow",
                                                "royalblue4",
                                                "yellow3",
                                                "darkorchid1",
                                                "lawngreen",
                                                "magenta3"))
    return(base_palette(num_colors))
  }