# --------------------------------------------------
# UTILS
#
# pspline_smooth: Smooths a 3D curve using P-splines
# reflect_point: Reflects a point across a plane
# --------------------------------------------------

# --------------------------------------------------
#  pspline_smooth
# --------------------------------------------------
#' Smooths a 3D curve using P-splines
#'
#' This function fits a cubic P-spline to each coordinate (X, Y, Z) of a 3D curve
#' and generates a smooth version of the curve with a specified number of points.
#' Use for smoothing noisy or/and resampling semi-landmark curves.
#'
#' @param points A numeric matrix of size n × 3 representing the 3D coordinates of the curve.
#' @param n_out Integer. Number of output points along the smoothed curve. Default is 100.
#' @param k Integer. Number of basis functions for the spline (controls smoothness). Default is `7`.
#'
#' @return A numeric matrix of size n_out × 3 with smoothed X, Y, Z coordinates.
#' @export
pspline_smooth <- function(points, n_out = 100, k = 7) {
  t <- seq(0, 1, length.out = nrow(points))
  t_out <- seq(0, 1, length.out = n_out)
  smoothed <- sapply(1:3, function(d) {
    fit <- mgcv::gam(points[, d] ~ s(t, bs = "ps", k = k))
    predict(fit, newdata = data.frame(t = t_out))
  })
  matrix(smoothed, ncol = 3, dimnames = list(NULL, c("X", "Y", "Z")))
}

# --------------------------------------------------
#  reflect_point
# --------------------------------------------------
#' Reflects a point across a plane
#'
#' This function computes the reflection of a 3D point/landmark across a plane.
#' The plane is defined by a point `origin` and a normal vector `pc3`.
#'
#' @param X Numeric vector of length 3. The point to reflect.
#' @param origin Numeric vector of length 3. A point on the plane.
#' @param pc3 Numeric vector of length 3. Normal vector of the plane (should be normalized).
#'
#' @return Numeric vector of length 3 representing the reflected point.
#' @export
#'
#' @examples
#' origin <- c(0, 0, 0)
#' normal <- c(0, 0, 1)   # plane XY
#' point <- c(1, 2, 3)
#' reflected <- reflect_point(point, origin, normal)
#' # reflected should be (1, 2, -3)

reflect_point <- function(X, origin, pc3) {
  d <- sum((X - origin) * pc3)       # signed distance from point to plane
  X - 2 * d * pc3                    # reflected point across plane
}
