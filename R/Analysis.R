# --------------------------------------------------
# ANALYSIS
#
# Proc_peristome : Perform a Procrustes alignment on peristome curves
# find_localmaxima : Detect local curvature maxima along 2D or 3D curves
# find_i2 : Detect i2 landmark position along each curve
# --------------------------------------------------


# --------------------------------------------------
#  Function: Proc_peristome
# --------------------------------------------------
#' Perform a Procrustes alignment on peristome curves
#'
#' @description
#' This function performs a Generalized Procrustes Analysis (GPA) on a set of
#' apertural curves (e.g., ammonite peristomes) extracted from multiple specimens.
#' It aligns all curves to a common coordinate system, optionally allowing
#' sliding landmarks to reduce tangential variation along curves.
#'
#' @param data A list-like object containing 3D coordinates of landmarks for each specimen.
#' It must have the structure:
#' `data$landmarks[[specimen_id]][[aperture_number]]`, where each aperture
#' is a matrix of size `(n_landmarks x n_dimensions)`.
#'
#' @param slide Logical, whether to perform sliding of semi-landmarks along curves
#' during GPA (default = `FALSE`).
#'
#' @param fixed_ldm Numeric vector indicating which landmarks are fixed
#' (i.e., not allowed to slide). Used only when `slide = TRUE`.
#'
#' @details
#' The function first extracts all apertural curves from each specimen, stores them in
#' a unified 3D array `(landmarks x dimensions x apertures)`, and then performs
#' a Procrustes superimposition using `geomorph::gpagen()`.
#'
#' If `slide = TRUE`, the function defines the semi-landmarks and their curve
#' connections using `curve_matrix`, allowing landmark sliding to minimize bending energy.
#'
#' After alignment, the function also reconstructs a list grouping the aligned apertures
#' by specimen (useful for downstream visualization or morphometric analysis).
#'
#' @return A `gpagen` object containing:
#' \itemize{
#'   \item `coords` - 3D array of aligned landmark coordinates.
#'   \item `Csize` - centroid sizes of each aperture.
#'   \item `consensus` - mean shape of all apertures.
#'   \item and other standard outputs from `geomorph::gpagen`.
#' }
#'
#' @examples
#' \dontrun{
#' gpa <- Proc_peristome(reconstructed_data, slide = TRUE, fixed_ldm = c(1, 200, 399))
#' summary(gpa)
#' plot(gpa)
#' }
#'
#' @importFrom geomorph gpagen
#' @importFrom stats complete.cases
#' @export

Proc_peristome <- function(data, slide = FALSE, fixed_ldm = c()) {

  all_apertures <- list()
  specimen_map <- c()  # Keeps track of which specimen each aperture belongs to

  # ---- STEP 1: Collect all apertural curves from all specimens ----
  for (specimen_id in names(data$landmarks)) {
    specimen_data <- data$landmarks[[specimen_id]]

    # Extract aperture names that are purely numeric (e.g., "1", "2", "3", ...)
    aperture_names <- grep("^\\d+$", names(specimen_data), value = TRUE)
    apertures <- lapply(aperture_names, function(name) as.matrix(specimen_data[[name]]))

    # If the specimen has apertures, append them to the global list
    if (length(apertures) > 0) {
      all_apertures <- c(all_apertures, apertures)
      specimen_map <- c(specimen_map, rep(specimen_id, length(apertures)))
    }
  }

  # ---- STEP 2: Check that we have enough curves to align ----
  if (length(all_apertures) < 2) {
    stop("At least 2 curves are needed to perform a Procrustes analysis.")
  }

  num_landmarks <- nrow(all_apertures[[1]])
  num_dimensions <- ncol(all_apertures[[1]])
  num_apertures <- length(all_apertures)

  # ---- STEP 3: Build a 3D array (landmarks x dimensions x apertures) ----
  apertures_array <- array(NA, dim = c(num_landmarks, num_dimensions, num_apertures))
  for (i in seq_along(all_apertures)) {
    apertures_array[, , i] <- all_apertures[[i]]
  }

  # ---- STEP 4: Perform GPA (with or without sliding landmarks) ----
  if (slide) {
    # Define fixed and sliding landmarks
    fixed_lm <- fixed_ldm
    sliding_lm <- setdiff(seq_len(num_landmarks), fixed_lm)

    # Build curve connections between adjacent semi-landmarks
    curve_matrix <- cbind(
      sliding_lm[-length(sliding_lm)],
      sliding_lm[-1],
      c(sliding_lm[-c(1,2)], NA)
    )
    curve_matrix <- curve_matrix[complete.cases(curve_matrix), ]

    # Perform GPA with sliding landmarks
    gpa_result <- gpagen(apertures_array, PrinAxes = TRUE, curves = curve_matrix)
  } else {
    # Perform standard GPA without sliding
    gpa_result <- gpagen(apertures_array, PrinAxes = TRUE)
  }

  # ---- STEP 5: Reorganize aligned data by specimen (optional, for downstream use) ----
  aligned_data <- list()
  for (specimen_id in unique(specimen_map)) {
    indices <- which(specimen_map == specimen_id)
    aligned_data[[specimen_id]] <- gpa_result$coords[, , indices]
  }

  # Return the complete GPA object (the aligned data can be extracted if needed)
  return(gpa_result)
}


# --------------------------------------------------
#  Function: find_localmaxima
# --------------------------------------------------

#' Detect local curvature maxima along 2D or 3D curves
#'
#' @description
#' This function computes local curvature profiles along left/right aperture curves,
#' optionally projected into 2D planes (lateral or frontal). Curvature is estimated
#' from smoothed spline derivatives, and local maxima are identified as points of
#' strongest bending along the curve.
#'
#' @param data A nested list containing landmark data for multiple specimens.
#'        Each specimen must include a list of curves (including "V" for ventral reference points).
#' @param k Integer. Number of ventral landmarks used for defining the local PCA orientation
#'        of each curve (default: 5).
#' @param spar Numeric. Smoothing parameter for spline fitting of coordinate trajectories,
#'        between (0, 1]; see [stats::smooth.spline()] (default: 0.8). Decreasing in spar value
#'        results to decrease of spotted maxima.
#' @param plot Logical. If `TRUE`, the function visualizes each curve and its detected maxima
#'        in a 3D interactive window using **rgl** (default: FALSE).
#' @param mode Character. Defines how curves are analyzed:
#'        \itemize{
#'          \item `"3D"` - curvature computed directly in 3D space from smoothed x,y,z splines.
#'          \item `"2D_lateral"` - curvature computed in the local PCA plane (PC1-PC2).
#'          \item `"2D_frontal"` - curvature computed in an orthogonal projection (PC1-PC3).
#'        }
#'
#' @details
#' For each specimen and curve:
#' \itemize{
#'   \item A local PCA is computed on the `k` ventral landmarks nearest to the curve's median point.
#'   \item The curve is projected into the PCA reference frame.
#'   \item Smooth splines are fitted independently to each coordinate dimension
#'         (`x(t)`, `y(t)`, and `z(t)` if applicable).
#'   \item First and second derivatives are computed from the spline fits.
#'   \item Curvature is derived from differential geometry formulas:
#'         \itemize{
#'           \item In 2D: kappa = (x' y'' - y' x'') / (x'^2 + y'^2)^(3/2)
#'           \item In 3D: kappa = ||(r'(t) * r''(t))|| / ||r'(t)||^3
#'         }
#'   \item The curvature **sign** is preserved by projecting the local normal
#'         vector (`r'(t) * r''(t)`) onto a mean reference direction.
#'   \item Local maxima are detected as curvature peaks (points where the
#'         first derivative of curvature changes from positive to negative).
#' }
#'
#' @return A nested list of results, structured as:
#' \itemize{
#'   \item `curvature` - vector of curvature magnitudes along the curve.
#'   \item `curvature_sign` - signed curvature preserving concavity information.
#'   \item `maxima_indices` - integer indices of detected local curvature maxima.
#'   \item `projection` - the 2D or 3D projection of the analyzed curve.
#' }
#'
#' @examples
#' \dontrun{
#' # Detect curvature maxima in 3D mode
#' results_3d <- find_localmaxima(reconstructed_data, plot = TRUE, mode = "3D", spar = 0.5)
#'
#' # Or in 2D lateral projection
#' results_lat <- find_localmaxima(reconstructed_data, mode = "2D_lateral")
#'
#' # Or in 2D frontal projection
#' results_front <- find_localmaxima(reconstructed_data, mode = "2D_frontal")
#' }
#'
#' @importFrom stats smooth.spline predict
#' @importFrom pracma cross
#' @importFrom rgl lines3d points3d clear3d
#' @export

find_localmaxima <- function(data, k = 5, spar = 0.8,  plot = FALSE, mode = c("2D_lateral", "2D_frontal", "3D")) {
  mode <- match.arg(mode)
  results <- list()

  for (specimen_id in names(data$landmarks)) {
    specimen <- data$landmarks[[specimen_id]]
    ventral_points <- as.matrix(specimen$V)
    origin <- colMeans(ventral_points)

    specimen_results <- list()

    for (curve_name in names(specimen)) {
      # Skip ventral and umbilical reference curves
      if (grepl("^[UV]", curve_name)) next

      curve <- specimen[[curve_name]]

      # Get a central landmark for neighborhood selection
      if (grepl("^L", curve_name)) {
        if (grepl("^R", curve_name)) {
          median_landmark <- curve[nrow(curve), ]
        } else {
          median_landmark <- curve[nrow(curve) / 2, ]
        }
      } else {
        median_landmark <- curve[1, ]
      }

      median_landmark <- as.numeric(median_landmark)

      # Find the k closest ventral landmarks to the median landmark
      distances <- sqrt(rowSums((t(t(ventral_points) - median_landmark))^2))
      nearest_indices <- order(distances)[1:k]
      selected_ventral <- ventral_points[nearest_indices, , drop = FALSE]

      # Perform local PCA on ventral points (always in 3D)
      selected_ventral_centered <- selected_ventral - matrix(origin, nrow(selected_ventral), 3, byrow = TRUE)
      ventral_pca <- prcomp(selected_ventral_centered, center = FALSE, scale. = FALSE)

      # Project the curve into the PCA space
      curve_centered <- curve - matrix(origin, nrow(curve), 3, byrow = TRUE)
      curve_proj <- predict(ventral_pca, newdata = curve_centered)

      # Extract the projection coordinates based on the selected mode
      if (mode == "2D_lateral") {
        coords <- curve_proj[, 1:2]  # PC1 vs PC2
      } else if (mode == "2D_frontal") {
        coords <- curve_proj[, c(1,3)]  # PC1 vs PC3
      } else {
        coords <- curve_proj  # PC1, PC2, PC3
      }

      # Compute curvature
      t <- seq(0, 1, length.out = nrow(coords))

      if (mode != "3D") {
        # --- Mode 2D lateral and 2D frontal ---
        # Fit spline curve on 2D points
        spline_x <- smooth.spline(t, coords[,1], spar = spar)
        spline_y <- smooth.spline(t, coords[,2], spar = spar)

        # First derivative
        dx <- predict(spline_x, t, deriv = 1)$y
        dy <- predict(spline_y, t, deriv = 1)$y

        # Second derivative
        d2x <- predict(spline_x, t, deriv = 2)$y
        d2y <- predict(spline_y, t, deriv = 2)$y

        # 2D curvature formula
        curvature <- abs(dx * d2y - dy * d2x) / (dx^2 + dy^2)^(3/2)
        curvature_sign <- (dx * d2y - dy * d2x) / (dx^2 + dy^2)^(3/2)

      } else {
        # --- Mode 3D ---
        spline_x <- smooth.spline(t, coords[,1], spar = spar)
        spline_y <- smooth.spline(t, coords[,2], spar = spar)
        spline_z <- smooth.spline(t, coords[,3], spar = spar)

        dx  <- predict(spline_x, t, deriv = 1)$y
        dy  <- predict(spline_y, t, deriv = 1)$y
        dz  <- predict(spline_z, t, deriv = 1)$y

        d2x <- predict(spline_x, t, deriv = 2)$y
        d2y <- predict(spline_y, t, deriv = 2)$y
        d2z <- predict(spline_z, t, deriv = 2)$y

        T1 <- cbind(dx, dy, dz)
        T2 <- cbind(d2x, d2y, d2z)

        # Classical 3D curvature
        cross_prod <- pracma::cross(T1, T2)
        num <- sqrt(rowSums(cross_prod^2))
        denom <- rowSums(T1^2)^(3/2)
        curvature <- num / denom

        # Signed curvature: project cross(T1,T2) on PCA normal vector
        normal_vec <- pracma::cross(T1, T2)
        normal_vec <- normal_vec / sqrt(rowSums(normal_vec^2))
        ref_vec <- colMeans(normal_vec, na.rm = TRUE) # reference orientation
        sign_curv <- sign(rowSums(normal_vec * matrix(ref_vec, nrow = nrow(normal_vec), ncol = 3, byrow = TRUE)))
        curvature_sign <- curvature * sign_curv
      }

      # Detect local maxima of curvature
      maxima_indices <- sort(which(diff(sign(diff(curvature))) == -2) + 1)

      if (plot == TRUE) {
        rgl::lines3d(curve, col = "black", lwd = 2,)
        if (length(maxima_indices) > 0) {
          rgl::points3d(curve[maxima_indices, , drop = FALSE], col = "red", size = 8)
        }
      }

      # Store results
      specimen_results[[curve_name]] <- list(
        curvature = curvature,
        curvature_sign = curvature_sign,
        maxima_indices = maxima_indices,
        projection = coords
      )
    }
    if (plot == TRUE) {
      readline(prompt = paste0("\u2705 Check the specimen ", specimen_id, " and press [Enter] to continue..."))
      clear3d()
    }

    results[[specimen_id]] <- specimen_results
  }
  close3d()

  return(results)
}


# --------------------------------------------------
#  Function: find_i2
# --------------------------------------------------

#' Detect i2 landmark position along each curve
#'
#' @description
#' Identify the index `i2` corresponding to the main convex curvature maximum
#' along each left/right aperture curve. If no suitable maximum is found in the
#' expected interval, the midpoint of the interval is used as a fallback.
#'
#' @param results Output list from [find_localmaxima()], containing curvature,
#'        maxima indices, and curve coordinates for each specimen and curve.
#'
#' @return A nested list of i2 indices for each specimen and curve.
#'
#' @details
#' For each curve:
#' \itemize{
#'   \item Expected search interval for i2:
#'         \itemize{
#'           \item `"R"` curves: indices 45-110
#'           \item `"L"` curves: indices 90-160
#'         }
#'   \item If no local maxima are found in this interval, the midpoint
#'         of the interval is used instead.
#'   \item Among convex maxima (positive curvature), the one closest
#'         to the ventral reference point is selected.
#' }
#' @examples
#' \dontrun{
#'
#' # Examples with maxima in lateral projection
#' results_lat <- find_localmaxima(data = data_missingsym, mode = "2D_lateral")
#'
#' i2_list <- find_i2(results_lateral)
#' }

#' @export

find_i2 <- function(results) {
  i2_points <- list()

  for (specimen_id in names(results)) {
    specimen_res <- results[[specimen_id]]
    specimen_i2 <- list()

    for (curve_name in names(specimen_res)) {
      curvature <- specimen_res[[curve_name]]$curvature
      maxima_indices <- specimen_res[[curve_name]]$maxima_indices
      coords <- specimen_res[[curve_name]]$projection  # projected coordinates

      # n_ldm <- length(coords)

      # Define search interval depending on side
      if (grepl("R", curve_name)) {
        min_idx <- 40
        max_idx <- 110
        ref_point <- coords[1, ]
      } else if (grepl("L", curve_name)) {
        min_idx <- 90
        max_idx <- 160
        ref_point <- coords[nrow(coords), ]
      } else {
        specimen_i2[[curve_name]] <- NA
        warning(sprintf("Specimen %s curve %s: side not recognized, i2 set to NA", specimen_id, curve_name))
        next
      }

      # Filter maxima within interval
      valid_maxima <- maxima_indices[maxima_indices >= min_idx & maxima_indices <= max_idx]

      if (length(valid_maxima) == 0) {
        # No maxima: use midpoint fallback
        i2_idx <- round((min_idx + max_idx) / 2)
        message(sprintf("\u26A0 Specimen %s curve %s: no maxima found, using midpoint index %d as fallback",
                        specimen_id, curve_name, i2_idx))
        specimen_i2[[curve_name]] <- i2_idx
        next
      }

      # Keep convex maxima only
      convex_maxima <- valid_maxima[curvature[valid_maxima] > 0]

      if (length(convex_maxima) == 0) {
        # No convex maxima: use midpoint fallback
        i2_idx <- round((min_idx + max_idx) / 2)
        message(sprintf("\u26A0 Specimen %s curve %s: no convex maxima found, using midpoint index %d as fallback",
                        specimen_id, curve_name, i2_idx))
        specimen_i2[[curve_name]] <- i2_idx
        next
      }

      # If multiple convex maxima: pick the one closest to ventral reference
      if (length(convex_maxima) == 1) {
        i2_idx <- convex_maxima
      } else {
        distances <- sapply(convex_maxima, function(i) sum((coords[i, ] - ref_point)^2))
        i2_idx <- convex_maxima[which.min(distances)]
      }

      specimen_i2[[curve_name]] <- i2_idx
    }

    i2_points[[specimen_id]] <- specimen_i2
  }

  return(i2_points)
}



# --------------------------------------------------
#  Function: view_i2
# --------------------------------------------------

#' Plot i2 landmarks on 3D specimen curves
#'
#' @description
#' Visualize the position of i2 landmarks on 3D ammonite apertures.
#' Each specimen is shown sequentially, displaying all apertural curves in black
#' and highlighting the i2 points in red.
#'
#' @param results List returned by [find_localmaxima()], containing
#'   per-specimen curve coordinates or other analysis outputs.
#' @param i2_points List returned by [find_i2()], containing i2 landmark indices
#'   for each aperture of each specimen.
#' @param data Original landmark dataset, typically of the form
#'   `list(landmarks = list(SpecimenID = list(CurveName = matrix(X, Y, Z))))`.
#'
#' @details
#' The function opens a 3D `rgl` window for each specimen.
#' All aperture curves are plotted in black (`lines3d`), and i2 landmarks are
#' displayed as large red points (`points3d`).
#' Press **Enter** in the R console to move to the next specimen.
#'
#' @importFrom rgl title3d lines3d points3d clear3d
#'
#' @examples
#' \dontrun{
#' results_lat <- find_localmaxima(data, mode = "2D_lateral")
#' i2_list <- find_i2(results_lat)
#' view_i2(results_lat, i2_list, data)
#' }
#'
#' @export

view_i2 <- function(results,
                    i2_points,
                    data) {
  for (specimen_id in names(results)) {
    specimen <- results[[specimen_id]]
    specimen_landmarks <- data$landmarks[[specimen_id]]

    # Initialize 3D window
    title3d(main = specimen_id)

    # Draw all curves
    for (curve_name in names(specimen_landmarks)) {
      curve <- specimen_landmarks[[curve_name]]
      if (!is.null(curve)) {
        lines3d(curve, col = "black", lwd = 2)
      }
    }

    # Draw i2 points
    specimen_i2 <- i2_points[[specimen_id]]
    for (curve_name in names(specimen_i2)) {
      i2_idx <- specimen_i2[[curve_name]]

      if (!is.na(i2_idx)) {
        curve <- specimen_landmarks[[curve_name]]
        if (!is.null(curve) && i2_idx >= 1 && i2_idx <= nrow(curve)) {
          points3d(curve[i2_idx, , drop = FALSE], col = "red", size = 10)
        }
      }
    }

    # Wait for user input before continuing
    readline(prompt = paste0("\U1F50D Check specimen ", specimen_id, ". Press [Enter] to continue..."))

    clear3d()
  }
}


