# --------------------------------------------------
# CLASSIC DESCRIPTORS
#
# Projection_2D_gpa : Computes a Generalized Procrustes Analysis of 2D projected peristomes
# project_on_plane : Projects 3D aperture curves onto a symmetry plane
# --------------------------------------------------


# --------------------------------------------------
#  Function: Projection_2D_gpa
# --------------------------------------------------

#' Computes a Generalized Procrustes Analysis of 2D projected peristomes
#'
#' Performs a Generalized Procrustes Analysis (GPA) on a list of 2D curves
#' (e.g., projected apertures or ribs) to align them for morphometric analysis.
#' Optionally, semi-landmarks can be allowed to slide along curves to minimize bending energy.
#'
#' @param projection A nested list of projected curves.
#'   Each element corresponds to a specimen and contains a list of apertures/curves
#'   as 2-column matrices (X, Y coordinates), typically from \link{project_on_plane} outputs.
#' @param slide Logical, default \code{FALSE}. If \code{TRUE}, semi-landmarks are allowed to slide along curves.
#' @param fixed_ldm Numeric vector, indices of landmarks to keep fixed during sliding (used only if \code{slide = TRUE}).
#'
#' @return A list containing:
#' \describe{
#'   \item{gpa}{The \code{geomorph::gpagen} object with aligned coordinates.}
#'   \item{coords}{The aligned coordinates array (p × 2 × n) for further analysis.}
#'   \item{meta}{A data frame linking each aligned curve to its specimen and aperture ID.}
#'   \item{meanshape}{The consensus (mean) shape of all aligned curves.}
#' }
#'
#' @details
#' This function is designed for 2D morphometric projections of 3D structures.
#' Each curve should be a 2-column matrix representing coordinates in a chosen plane
#' (e.g., PC1 vs PC2 from a prior 3D projection). If \code{slide = TRUE}, the function
#' constructs semi-landmark curves between fixed landmarks to allow sliding along tangents,
#' improving alignment and reducing artificial bending.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' projected_ribs <- project_on_plane(data, mode = "2D_lateral")
#' res_gpa <- Projection_2D_gpa(projected_ribs, slide = TRUE, fixed_ldm = c(1, 25))
#' plot(res_gpa$gpa)
#' }
#'
#' @export

Projection_2D_gpa <- function(projection, slide = FALSE, fixed_ldm = c()) {

  # =============================================
  # 1) Validate input and extract 2D matrices
  #    Input structure: projection$specimen_id$aperture_id = 2-column matrix (X,Y)
  # =============================================
  proj_list <- lapply(projection, function(spec){
    lapply(spec, function(aperture_matrix){
      # Each aperture must be a 2-column matrix
      if(!is.matrix(aperture_matrix) || ncol(aperture_matrix) != 2)
        stop("Each aperture must be a 2-column matrix (X, Y coordinates).")
      return(aperture_matrix)
    })
  })

  # Flatten the nested list to a single list of matrices
  proj_list_flat <- unlist(proj_list, recursive = FALSE)

  # Check that there are at least two curves for GPA
  if(length(proj_list_flat) < 2)
    stop("At least 2 curves are required for GPA alignment.")

  # Get number of landmarks (p), dimension (k = 2), and number of curves (n)
  p <- nrow(proj_list_flat[[1]])
  k <- ncol(proj_list_flat[[1]])
  n <- length(proj_list_flat)

  # Stack the matrices into a 3D array: landmarks × dimensions × curves
  projections_array <- array(unlist(proj_list_flat),
                             dim = c(p, k, n))

  # =============================================
  # 2) Perform GPA
  #    Optionally slide semi-landmarks along curves to minimize bending energy
  # =============================================
  if(slide){
    # Sliding is applied to landmarks that are NOT fixed
    sliding <- setdiff(seq_len(p), fixed_ldm)

    # Define curve connections for sliding semi-landmarks
    # Each row connects three consecutive points (start, end, and midpoint) along the curve
    curve_matrix <- cbind(sliding[-length(sliding)],
                          sliding[-1],
                          c(sliding[-c(1,2)], NA))
    curve_matrix <- curve_matrix[complete.cases(curve_matrix),]

    # Perform GPA with sliding semi-landmarks
    gpa <- geomorph::gpagen(projections_array,
                            PrinAxes = TRUE,
                            curves = curve_matrix)
  } else {
    # Standard GPA without sliding
    gpa <- geomorph::gpagen(projections_array, PrinAxes = TRUE)
  }

  # =============================================
  # 3) Reconstruct metadata for each curve
  # =============================================
  # Recreate specimen and aperture identifiers corresponding to each curve
  specimen_names <- rep(names(projection),
                        times = sapply(projection, length))
  aperture_names <- unlist(lapply(projection, names))

  meta <- data.frame(
    index = seq_len(n),
    specimen = specimen_names,
    aperture = aperture_names
  )

  # =============================================
  # 4) Return results
  # =============================================
  return(list(
    gpa = gpa,            # GPA object with aligned shapes
    coords = gpa$coords,  # aligned coordinates array (p × 2 × n)
    meta = meta,          # metadata linking curves to specimens
    meanshape = gpa$consensus  # mean shape of all curves
  ))
}


# --------------------------------------------------
#  Function: project_on_plane
# --------------------------------------------------
#' Projects 3D aperture curves onto a symmetry plane
#'
#' @description
#' For each specimen, this function projects every non-ventral aperture curve
#' into a local geometric reference frame defined from ventral landmarks (\link{compute_local_plane}).
#' Projection may be returned in full 3D coordinates or reduced to 2D views.
#' In 2D modes, curves are automatically resampled so that points follow
#' **equal arc-length spacing**, making curves comparable independent of size
#' or deformation.
#'
#' @param data A hierarchical list containing landmark coordinates, typically from \link{reconstruct_missingldm}.
#' It must include a \code{$landmarks} element, where each specimen contains:
#' \itemize{
#' \item a ventral reference curve named \code{"V"},
#' \item several aperture curves named with numeric identifiers.
#' }
#'
#' @param k Integer. Size of the local tangent window (number of points) used by
#' \link{compute_local_plane} to estimate the ventral direction.
#'
#' @param mode Character string specifying the projection output:
#' \describe{
#' \item{\code{"3D"}}{Returns full 3D coordinates in the local frame
#' (PC1, PC2, PC3).}
#' \item{\code{"2D_lateral"}}{Returns a lateral view (PC1–PC2) of the
#' ventral half of the aperture, resampled to equal arc-length spacing.}
#' \item{\code{"2D_frontal"}}{Returns a frontal view (PC1–PC3) of the
#' full aperture, resampled to equal arc-length spacing.}
#' }
#'
#' @param weighting Character. Passed to \link{compute_local_plane}, controls
#' whether ventral landmarks are weighted during local plane estimation.
#'
#' @param plot Logical. If \code{TRUE}, displays the aperture curve, ventral
#' landmarks, and local reference axes in an interactive \code{rgl} window
#' for diagnostic purposes.
#'
#' @details
#' For each specimen and each aperture curve, the function performs the following steps:
#' \enumerate{
#' \item Identifies the median landmark of the aperture curve.
#' \item Computes a local orthonormal basis (PC1, PC2, PC3) from ventral landmarks
#' using \code{compute_local_plane()}.
#' \item Centers the aperture on the local origin (median landmark).
#' \item Projects the curve onto the local axes.
#' \item Optionally trims the curve to its ventral half in lateral mode.
#' \item Resamples 2D projections to uniform arc-length spacing using
#' linear interpolation.
#' }
#'
#' No global alignment or Procrustes superimposition is performed; all projections
#' are purely local and specimen-specific.
#'
#' @return
#' A nested list structured as:
#' \itemize{
#' \item \code{result[[specimen_id]][[curve_id]]} → matrix of projected coordinates:
#' \itemize{
#' \item \code{n x 3} matrix for \code{"3D"} mode,
#' \item \code{n x 2} matrix for \code{"2D_lateral"} and \code{"2D_frontal"} modes.
#' }
#' }
#'
#' @examples
#' \dontrun{
#' # 3D projection
#' proj_3d   <- project_on_plane(data, mode="3D")
#'
#' # Lateral PC1–PC2 view, equalized spacing
#' proj_lat <- project_on_plane(data, mode="2D_lateral")
#'
#' # Frontal PC1–PC3 view, equalized spacing
#' proj_frn <- project_on_plane(data, mode="2D_frontal")
#' }
#'
#' @seealso \link{compute_local_plane}
#' @export

project_on_plane <- function(data,
                             k = 5,
                             mode = c("2D_lateral", "2D_frontal", "3D"),
                             weighting = "none",
                             plot = FALSE) {

  mode <- match.arg(mode)
  results <- list()

  for (specimen_id in names(data$landmarks)) {

    specimen <- data$landmarks[[specimen_id]]
    ventral_points <- as.matrix(specimen$V)

    projection_coords <- list()

    # We loop through openings in numeric order
    opening_names <- names(specimen)
    opening_names <- opening_names[grepl("^[0-9]+$", opening_names)]
    opening_names <- as.character(sort(as.numeric(opening_names)))

    for (curve_name in opening_names) {

      curve <- specimen[[curve_name]]
      n_ldm <- nrow(curve)

      # landmark median
      median_ldm <- as.numeric(curve[(n_ldm + 1) / 2, ])

      # local plane
      plane <- compute_local_plane(
        V = ventral_points,
        median_ldm = median_ldm,
        origin = median_ldm,        # projection origin
        weighting = weighting,
        tangent_window = k,
        stabilized_direction = TRUE
      )

      origin <- plane$origin
      pc1 <- plane$pc1
      pc2 <- plane$pc2
      pc3 <- plane$pc3

      # ------------------------------------
      # 2D mode adjustments
      # ------------------------------------
      curve_use <- curve

      if (mode == "2D_lateral") {
        # keep ventral half
        n_ldm <- nrow(curve)
        curve_use <- curve_use[1:round((n_ldm + 1)/2), ]
      }

      # ------------------------------------
      # Projection
      # ------------------------------------
      curve_centered <- curve_use -
        matrix(origin, nrow(curve_use), 3, byrow = TRUE)

      # Full PCA projection matrix
      proj_mat <- cbind(pc1, pc2, pc3)
      curve_proj <- as.matrix(curve_centered) %*% proj_mat

      if (plot == TRUE){

        points3d(curve_use)
        points3d(ventral_points)

        L <- 5
        pc1_end <- origin + as.numeric(pc1)*L
        pc2_end <- origin + as.numeric(pc2)*L
        pc3_end <- origin + as.numeric(pc3)*L

        points3d(origin, col="red",  lwd=20)  # PC1

        points3d(pc1_end, col="green",  lwd=10)  # PC1
        points3d(pc2_end, col="yellow", lwd=10)  # PC2
        points3d(pc3_end, col="purple", lwd=10)  # PC3

      }

      # ------------------------------------
      # Resampling helper
      # ------------------------------------
      resampling_2D <- function(coords){
        d <- sqrt(rowSums(diff(coords)^2))
        s <- c(0, cumsum(d))
        s_norm <- s / max(s)

        target_s <- seq(0,1, length.out = nrow(coords))
        coords <- cbind(
          approx(s_norm, coords[,1], xout = target_s)$y,
          approx(s_norm, coords[,2], xout = target_s)$y
        )
      }

      # ------------------------------------
      # Format output by mode
      # ------------------------------------
      if (mode == "2D_lateral") {
        coords <- resampling_2D(curve_proj[, 1:2])  # PC1 vs PC2
      } else if (mode == "2D_frontal") {
        coords <- resampling_2D(curve_proj[, c(2, 3)])  # PC2 vs PC3
      } else {
        coords <- curve_proj  # full 3D
      }

      projection_coords[[curve_name]] <- coords
    }

    results[[specimen_id]] <- projection_coords

    if (plot == TRUE){
      Sys.sleep(10)
      close3d()
    }
  }

  return(results)
}
