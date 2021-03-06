% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coexpression.R
\name{get_bilinear_val}
\alias{get_bilinear_val}
\title{Bilinear Interpolation}
\usage{
get_bilinear_val(x, y, ngrid, quad11, quad21, quad12, quad22)
}
\arguments{
\item{x}{integer vector}

\item{y}{integer vector}

\item{ngrid}{integer taken from get_color_matrix_df}

\item{quad11}{integer}

\item{quad21}{integer}

\item{quad12}{integer}

\item{quad22}{integer}
}
\value{
numeric vector for color values used in coexpression plot
}
\description{
Takes in a given x and y coordinate, the dimension of a square grid, and four values representing the base red, green, or blue (RGB) color values (0-255) of the four quadrants of the square grid
Returns the bilinear interpolated red, green, or blue value (0-255) of the given input coordinates (x,y)
}
\examples{
ngrid <- 16
color_matrix_df <- expand.grid(x_value = 0:ngrid, y_value = 0:ngrid)
color10 <- c(255,0,0) # numeric vector of RGB values for red quadrant of 2d color matrix
color01 <- c(0,0,255) # numeric vector of RGB values for blue quadrant of 2d color matrix
color00 <- c(217,217,217) # numeric vector of RGB values for light gray quadrant of 2d color matrix
color11 <- c(255,0,255) # numeric vector of RGB values for pink/violet quadrant of 2d color matrix
color_matrix_df$red_values <- get_bilinear_val(color_matrix_df$x_value, color_matrix_df$y_value, ngrid, color00[1], color10[1], color01[1], color11[1])
color_matrix_df$green_values <- get_bilinear_val(color_matrix_df$x_value, color_matrix_df$y_value, ngrid, color00[2], color10[2], color01[2], color11[2])
color_matrix_df$blue_values <- get_bilinear_val(color_matrix_df$x_value, color_matrix_df$y_value, ngrid, color00[3], color10[3], color01[3], color11[3])
}
