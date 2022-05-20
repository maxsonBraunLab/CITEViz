# ----- Bilinear interpolation function -----
# Takes in a given x and y coordinate, the dimension of a square grid, and four values representing the base red, green, or blue (RGB) color values (0-255) of the four quadrants of the square grid
# Returns the bilinear interpolated red, green, or blue value (0-255) of the given input coordinates (x,y)
#' Title
#'
#' @param x 
#' @param y 
#' @param ngrid 
#' @param quad11 
#' @param quad21 
#' @param quad12 
#' @param quad22 
#'
#' @return
#' @export
#'
#' @examples
get_bilinear_val <- function(x,y,ngrid,quad11,quad21,quad12,quad22){
  temp_val <- quad11*(ngrid-x)*(ngrid-y) + quad21*x*(ngrid-y) + quad12*(ngrid-x)*y + quad22*x*y
  bilinear_val <- temp_val / (ngrid*ngrid)
  return(bilinear_val)
}


# ----- Create 2D color dataframe for gene/ADT coexpression -----
#' Title
#'
#' @param ngrid 
#'
#' @importFrom grDevices rgb
#'
#' @return
#' @export
#'
#' @examples
get_color_matrix_df <- function(ngrid = 16) {
  color_matrix_df <- expand.grid(x_value = 0:ngrid, y_value = 0:ngrid)
  color10 = c(255,0,0) # numeric vector of RGB values for red quadrant of 2d color matrix
  color01 = c(0,0,255) # numeric vector of RGB values for blue quadrant of 2d color matrix
  color00 = c(217,217,217) # numeric vector of RGB values for light gray quadrant of 2d color matrix
  color11 = c(255,0,255) # numeric vector of RGB values for pink/violet quadrant of 2d color matrix
  color_matrix_df$R <- get_bilinear_val(color_matrix_df$x_value, color_matrix_df$y_value, ngrid, color00[1], color10[1], color01[1], color11[1])
  color_matrix_df$G <- get_bilinear_val(color_matrix_df$x_value, color_matrix_df$y_value, ngrid, color00[2], color10[2], color01[2], color11[2])
  color_matrix_df$B <- get_bilinear_val(color_matrix_df$x_value, color_matrix_df$y_value, ngrid, color00[3], color10[3], color01[3], color11[3])
  color_matrix_df$hex_color_mix <- rgb(color_matrix_df$R, color_matrix_df$G, color_matrix_df$B, maxColorValue = 255)
  
  return(color_matrix_df)
}

# ----- Create 2D color legend plot for gene/ADT coexpression -----
#' Title
#'
#' @param input 
#'
#' @return
#' @export
#'
#' @examples
create_2d_color_legend <- function(input, myso) {
  #selected metadata to color clusters by
  color_x <- input$x_axis_feature
  color_y <- input$y_axis_feature
  
  SeuratObject::DefaultAssay(myso) <- input$Assay_x_axis
  count_data_x <- SeuratObject::FetchData(object = myso, vars = color_x, slot = "data")
  # extract only the count values as a vector from the original count data dataframe
  count_data_x <- count_data_x[[color_x]]
  
  SeuratObject::DefaultAssay(myso) <- input$Assay_y_axis
  count_data_y <- SeuratObject::FetchData(object = myso, vars = color_y, slot = "data")
  # extract only the count values as a vector from the original count data dataframe
  count_data_y <- count_data_y[[color_y]]
  
  ngrid <- 16
  color_matrix_df <- get_color_matrix_df(ngrid)
  
  #show plot of 2D color legend
  color_matrix_df %>%
    ggplot2::ggplot(aes(x = x_value, y = y_value)) + 
    geom_tile(fill = color_matrix_df$hex_color_mix) +
    labs(x = input$x_axis_feature, y = input$y_axis_feature) +
    scale_x_continuous(breaks = c(0, ngrid), label = c(paste0("low\n", round(min(count_data_x), digits=2)), 
                                                       paste0("high\n", round(max(count_data_x), digits=2)))) + 
    scale_y_continuous(breaks = c(0, ngrid), label = c(paste0("low\n", round(min(count_data_y), digits=2)), 
                                                       paste0("high\n", round(max(count_data_y), digits=2)))) 
}
