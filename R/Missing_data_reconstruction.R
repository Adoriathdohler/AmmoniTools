
# --------------------------------------------------
#  MISSING DATA RECONSTRUCTION
#
# compute_local_plane : Computes a local symmetry plane from ammonite ventral landmarks
# symmetrize_sides : Symmetrizes and reconstructs missing portions of aperture and umbilical curves by exploiting bilateral symmetry of ammonite shells
# merge_sides : Merges left and right aperture sides into complete 3D peristome curves
# impute_weighted_TPS : Reconstructs missing landmarks coordinates using weighted TPS interpolation
# reconstruct_missingldm : Estimates missing landmark coordinates using distance-weighted TPS interpolation, prioritizing geometrically, ontogenetically and/or taxonomically closest reference apertures
# --------------------------------------------------

utils::globalVariables(c("Specimen", "name", "spec_tab"))


# --------------------------------------------------
#  Function: compute_local_plane
# --------------------------------------------------

#' Computes a local symmetry plane from ammonite ventral landmarks
#'
##' \strong{Internal function – not intended to be called directly by users.}
##'
#' This function computes a local orthonormal reference frame associated with
#' a single ammonite aperture, based on the geometry of the ventral curve.
#' The resulting frame is used internally to project aperture shapes
#' into consistent local 2D or 3D coordinate systems (e.g. frontal or lateral
#' projections).
#'
#' The local reference frame is defined by three orthogonal unit vectors:
#'
#' \itemize{
#'   \item \strong{PC1}: local growth direction, tangent to the ventral curve.
#'   \item \strong{PC2}: in-plane direction orthogonal to PC1, oriented away from
#'   the ammonite umbilicus to ensure consistent left–right orientation.
#'   \item \strong{PC3}: normal vector of the local plane, corresponding to the
#'   local bilateral symmetry plane.
#' }
#'
#' A weighted PCA of the ventral curve is used to estimate the plane geometry,
#' allowing either a global or locally constrained fit depending on the chosen
#' weighting scheme.
#'
#' @section Intended use:
#' This function is a low-level geometric utility used internally by higher-level
#' functions such as \link{symmetrize_sides} or \link{project_on_plane}.
#' Users should normally not call this function directly unless developing
#' or debugging methodological extensions.
#'
#' @param V Numeric matrix (\eqn{n \times 3}) containing the 3D coordinates
#' of the ventral curve.
#' @param median_ldm Numeric vector of length 3 giving the coordinates of the
#' median landmark of the aperture.
#' @param origin Numeric vector of length 3 used as the centering origin for PCA.
#' @param weighting Character string specifying the weighting scheme applied to
#' the ventral curve during PCA. One of:
#' \describe{
#'   \item{"none"}{Uniform weights (global plane estimation).}
#'   \item{"linear"}{Linear decay of weights with curvilinear distance from
#'   the median landmark.}
#'   \item{"exp"}{Exponential decay of weights, emphasizing local geometry.}
#' }
#' @param weight_power Numeric. Controls the strength of decay for
#' \code{"linear"} and \code{"exp"} weighting schemes.
#' @param tangent_window Integer. Number of ventral points on each side of the
#' median landmark used to estimate the local tangent direction when applicable.
#' @param global_center Numeric vector of length 3. Optional.
#' Approximate center of the ammonite used to orient PC2 consistently.
#' If \code{NULL}, the center is computed as the mean position of \code{V}.
#' @param stabilized_direction Logical.
#' If \code{TRUE}, PC1 is explicitly defined from the local ventral tangent,
#' enforcing a consistent biological growth direction.
#' If \code{FALSE}, PC1 is derived from the first eigenvector of the weighted PCA.
#'
#' @return
#' A list with the following components:
#' \describe{
#'   \item{\code{origin}}{Reference origin used for centering.}
#'   \item{\code{pc1}}{Unit vector corresponding to the local growth direction.}
#'   \item{\code{pc2}}{Unit vector lying in the local plane and oriented
#'   opposite to the ammonite umbilicus.}
#'   \item{\code{pc3}}{Unit normal vector of the local plane
#'   (local bilateral symmetry plane).}
#'   \item{\code{distances}}{Curvilinear distances of ventral points from
#'   the median landmark.}
#'   \item{\code{weights}}{Normalized weights applied to ventral points
#'   during the weighted PCA.}
#' }
#'
#' @details
#' When \code{stabilized_direction = TRUE}, the growth direction (PC1) is
#' explicitly derived from the ventral curve geometry, preventing residual
#' axis flips caused by eigenvector sign indeterminacy.
#'
#' PC3 is first estimated as the third eigenvector of the weighted PCA and then
#' orthogonalized relative to PC1. PC2 is computed as the cross product
#' of PC3 and PC1 and oriented toward the ammonite center to ensure a consistent
#' coordinate system across apertures.
#'
#' @seealso \link{symmetrize_sides} \link{project_on_plane}
#' @importFrom pracma cross
#' @keywords internal
#' @export
compute_local_plane <- function(V,
                                median_ldm,
                                origin,
                                weighting = c("none", "linear", "exp"),
                                weight_power = 1,
                                tangent_window = 5,
                                global_center = NULL,
                                stabilized_direction = FALSE) {

  weighting <- match.arg(weighting)
  median_ldm <- as.numeric(median_ldm)
  origin <- as.numeric(origin)

  # ---------------------------------------------------
  # 0. Compute global center (if not provided)
  #    → used to orient PC2 consistently
  # ---------------------------------------------------
  if (is.null(global_center)) {
    global_center <- colMeans(V)
  }

  # ---------------------------------------------------
  # 1. Compute curvilinear distances
  # ---------------------------------------------------
  dists_segment <- sqrt(rowSums(diff(V)^2))
  curv_dist <- c(0, cumsum(dists_segment))

  # Anchor = point in V closest to median landmark
  dist_to_median <- sqrt(rowSums((V - matrix(as.numeric(median_ldm), nrow(V), 3, byrow = TRUE))^2))
  anchor_idx <- which.min(dist_to_median)
  anchor_idx <- max(1, min(anchor_idx, nrow(V)))

  # ---------------------------------------------------
  # 2. PC1 = direction de croissance locale (toujours)
  # ---------------------------------------------------
  if(stabilized_direction == TRUE){
    # ---- Compute PC1 = local growth direction ---
    if (anchor_idx < nrow(V)) {
      next_idx <- anchor_idx + 1
      pc1 <- V[next_idx, ] - V[anchor_idx, ]
    } else {
      next_idx <- anchor_idx - 1
      pc1 <- -(V[next_idx, ] - V[anchor_idx, ])
    }
    pc1 <- pc1 / sqrt(sum(pc1^2))
  }


  # ---------------------------------------------------
  # 3. Weighted PCA on ventral curve → provisional PC3
  # ---------------------------------------------------
  distances <- abs(curv_dist - curv_dist[anchor_idx])

  w <- switch(weighting,
              none   = rep(1, length(distances)),
              linear = pmax(1 - (weight_power * distances / max(distances)), 0),
              exp    = exp(-weight_power * distances / max(distances)))

  W <- diag(w / sum(w))
  w_norm <- (w - min(w)) / (max(w) - min(w))

  centered <- V - matrix(origin, nrow(V), 3, byrow = TRUE)
  cov_w <- t(centered) %*% W %*% centered
  eig <- eigen(cov_w, symmetric = TRUE)

  if(stabilized_direction == FALSE){
    # PC1 of PCA = main direction of ventral curve (tangent)
    pc1 <- eig$vectors[, 1]
    pc1 <- pc1 / sqrt(sum(pc1^2))

    # --- Enforce biological growth direction ---
    if (anchor_idx < nrow(V)) {
      growth_vec <- V[anchor_idx + 1, ] - V[anchor_idx, ]
    } else {
      growth_vec <- V[anchor_idx, ] - V[anchor_idx - 1, ]
    }

    if (sum(pc1 * growth_vec) < 0) {
      pc1 <- -pc1
    }
  }

  # PCA chooses arbitrary sign → to be fixed
  pc3 <- eig$vectors[,3]

  # ---------------------------------------------------
  # 4. Make PC3 ⟂ PC1
  # ---------------------------------------------------
  pc3 <- pc3 - sum(pc3 * pc1) * pc1
  pc3 <- pc3 / sqrt(sum(pc3^2))

  # ---------------------------------------------------
  # 5. Compute PC2 from the orthonormal cross product
  # ---------------------------------------------------
  pc2 <- pracma::cross(pc3, pc1)
  pc2 <- pc2 / sqrt(sum(pc2^2))

  # ---------------------------------------------------
  # 6. ORIENT PC2 toward the center of ammonite
  # ---------------------------------------------------
  to_center <- global_center - median_ldm
  if (sum(pc2 * to_center) < 0) {
    pc2 <- -pc2
    pc3 <- -pc3
  }

  # Return
  list(
    origin    = origin,
    pc1       = as.numeric(pc1),
    pc2       = as.numeric(pc2),
    pc3       = as.numeric(pc3),
    distances = distances,
    weights   = w_norm
  )
}


# --------------------------------------------------
#  Function: symmetrize_sides
# --------------------------------------------------

#' Symmetrizes and reconstructs missing portions of aperture and umbilical curves by exploiting bilateral symmetry of ammonite shells
#'
#' This function reconstructs missing left/right (L/R) or umbilical (U/Ubis) curves
#' of an ammonite aperture by reflecting existing curves across a local symmetry plane.
#' When both sides exist, it performs a full geometric symmetrization using thin-plate
#' spline interpolation between corresponding landmarks. The symmetry plane is
#' estimated from ventral landmarks of the same specimen, optionally weighted by their
#' proximity to the ventral edge.
#'
#' @param data A list containing:
#'   \describe{
#'     \item{\code{landmarks}}{A nested list of 3D landmarks for each curve of each specimen.}
#'     \item{\code{curve_info}}{A data.frame describing curves (used to detect missing sides).}
#'     Needed data structure is typically from \link{reorder_landmarks} outputs.
#'   }
#' @param weighting Character, passed to \link{compute_local_plane}.
#'   Method used to weight ventral landmarks when computing the
#'   local symmetry plane. Options are:
#'   \describe{
#'     \item{\code{"none"}}{Uniform weights (all landmarks contribute equally), i.e. global symmetry plane.}
#'     \item{\code{"linear"}}{Weight decreases linearly with distance from the ventral edge.}
#'     \item{\code{"exp"}}{Weight decreases exponentially with distance from the ventral edge.}
#'   }
#'   Default is \code{"none"}.
#' @param weight_power  Numeric. Decay parameter controlling the effect of weighting
#'   (used for \code{"linear"} and \code{"exp"}). Default: \code{1}.
#' @param k_spline Integer. Number of knots used for P-spline smoothing of the ventral curve.
#'   Default: \code{7}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{curve_info}}{Updated curve information table including any new symmetrized curves.}
#'     \item{\code{landmarks}}{Updated landmarks list with reconstructed or symmetrized curves (L/R or U/Ubis).}
#'   }
#'
#' @details
#' For each specimen, the algorithm:
#' \enumerate{
#'   \item Smooths the ventral curve (\code{V}) with P-splines to obtain a stable geometric reference.
#'   \item Computes a local symmetry plane from weighted PCA of ventral points (weights depend on the chosen method).
#'   \item For each L/R pair:
#'         \itemize{
#'           \item If both sides exist -> performs full symmetric alignment using local reflection and TPS interpolation.
#'           \item If only one side exists -> reflects it across the local plane to reconstruct the opposite side.
#'         }
#'   \item For umbilical curves (\code{U}/\code{Ubis}), a global plane is estimated and reflection applied symmetrically.
#' }
#'
#' This ensures that local shell symmetry is respected and missing aperture halves are reconstructed consistently.
#'
#' @examples
#' \dontrun{
#' # Load dataset
#' data("example_spec_tab", package = "AmmoniTools")
#' data("example_data", package = "AmmoniTools")
#'
#' # Resample curve with 25 landmarks for each curve type
#' resampled_data <- resample_curves(example_data, n_landmarks = list(R = 25, L = 25, U = 25, V = 25))
#'
#' # Reorder curves
#' reordered_data <- reorder_landmarks(resampled_data)
#'
#' #Symmetrize sides
#' data_sym <- symmetrize_sides(reordered_data,
#'                              weighting = "exp",
#'                              weight_power  = 5,
#'                              k_spline = 7)
#' }
#'
#' @seealso \link{compute_local_plane}
#' @importFrom stats prcomp
#' @importFrom Morpho tps3d
#' @export

symmetrize_sides <- function(data, weighting = c("none", "linear", "exp"),
                             weight_power  = 5, k_spline = 7) {

  weighting <- match.arg(weighting)
  curve_info <- data$curve_info
  landmarks <- data$landmarks


  # ------------------------------------------------------------------
  # Helper: identify missing landmark indices for a given curve
  # ------------------------------------------------------------------
  get_missing_indices <- function(specimen_id, curve_id, n_ldm) {
    row <- which(curve_info$Specimen == specimen_id & curve_info$Curve == curve_id)
    if (length(row) == 0) return(integer(0))
    missing_vals <- curve_info$missing_ldm[[row]]
    # Cases with one missing region
    if (length(missing_vals) == 1) {
      if (missing_vals > (n_ldm/2)) return(unique(missing_vals[1]:n_ldm))
      else return(unique(1:missing_vals[1]))
    }
    # Cases with multiple missing landmarks: take one near each end
    if (length(missing_vals) >= 2) {
      mid <- n_ldm / 2
      first_half_vals <- missing_vals[missing_vals <= mid]
      second_half_vals <- missing_vals[missing_vals > mid]
      first_sel <- if (length(first_half_vals) > 0) min(first_half_vals) else 1
      second_sel <- if (length(second_half_vals) > 0) max(second_half_vals) else n_ldm
      return(unique(c(1:first_sel, second_sel:n_ldm)))
    }
    integer(0)
  }

  if (!"missing_ldm_updated" %in% names(curve_info)) {
    curve_info$missing_ldm_updated <- vector("list", nrow(curve_info))
  }

  # ==================================================================
  # Main loop: process all specimens
  # ==================================================================
  for (specimen_id in names(landmarks)) {
    specimen_data <- landmarks[[specimen_id]]
    if (!("V" %in% names(specimen_data))) next

    # Smooth ventral points
    V_raw <- as.matrix(specimen_data$V)
    ventral_points <- pspline_smooth(V_raw, n_out = nrow(V_raw), k = k_spline)
    landmarks[[specimen_id]][["V"]] <- ventral_points

    # ------------------------------------------------------------
    # 2. Process each pair of L/R curves
    # ------------------------------------------------------------
    lr_names <- grep("^[LR]", names(specimen_data), value = TRUE)
    nums <- unique(gsub("^[LR]", "", lr_names))

    for (num in nums) {

      L_name <- paste0("L", num)
      R_name <- paste0("R", num)
      L_exists <- L_name %in% names(specimen_data)
      R_exists <- R_name %in% names(specimen_data)

      # ============================================================
      # Case A: both sides exist -> perform full symmetric alignment
      # ============================================================
      if (L_exists && R_exists) {

        L_pts <- as.matrix(specimen_data[[L_name]])
        R_pts <- as.matrix(specimen_data[[R_name]])

        # Find median ldm
        median_ldm <- as.matrix(colMeans(rbind(L_pts[1,], R_pts[nrow(R_pts),])))

        # Find local plane
        plane <- compute_local_plane(V = ventral_points, median_ldm, origin = median_ldm,
                                     weighting, weight_power)
        origin <- plane$origin
        pc3 <- plane$pc3

        # --- proceed with symmetry operations using pc3 & origin_local  ---
        R_pts <- R_pts[nrow(R_pts):1, , drop = FALSE]  # reverse to match orientation
        n_ldm <- nrow(L_pts)
        if (nrow(R_pts) != n_ldm) next

        # Identify missing landmark regions from metadata
        missing_L <- get_missing_indices(specimen_id, L_name, n_ldm)
        missing_R <- get_missing_indices(specimen_id, R_name, n_ldm)

        # Adjust right-side indices (since it's reversed)
        invert_index <- function(idx, n) if (length(idx) > 0) n - idx + 1 else integer(0)
        missing_R <- invert_index(missing_R, n_ldm)

        # Determine indices for common and one-sided landmarks
        common_valid_idx <- setdiff(1:n_ldm, union(missing_L, missing_R))
        only_L_idx <- setdiff(setdiff(1:n_ldm, missing_L), common_valid_idx)
        only_R_idx <- setdiff(setdiff(1:n_ldm, missing_R), common_valid_idx)
        if (length(common_valid_idx) < 5) next  # not enough anchors for TPS

        # --- Reflect both sides across the local plane and average them
        L_common <- L_pts[common_valid_idx, , drop = FALSE]
        R_common <- R_pts[common_valid_idx, , drop = FALSE]
        L_reflected_common <- t(apply(L_common, 1, reflect_point, origin = origin, pc3 = pc3))
        R_reflected_common <- t(apply(R_common, 1, reflect_point, origin = origin, pc3 = pc3))
        R_after_common <- (R_common + L_reflected_common) / 2
        L_after_common <- (L_common + R_reflected_common) / 2

        # --- Initialize final coordinates
        L_final <- matrix(NA, nrow = n_ldm, ncol = 3)
        R_final <- matrix(NA, nrow = n_ldm, ncol = 3)
        colnames(L_final) <- c("X", "Y", "Z")
        colnames(R_final) <- c("X", "Y", "Z")

        L_final[common_valid_idx, ] <- L_after_common
        R_final[common_valid_idx, ] <- R_after_common

        # --- Interpolate missing landmarks using TPS
        if (length(only_L_idx) > 0) {
          tps_L_missing <- Morpho::tps3d(refmat = L_common, tarmat = L_after_common,
                                         x = L_pts[only_L_idx, , drop = FALSE])
          L_final[only_L_idx, ] <- tps_L_missing
          R_final[only_L_idx, ] <- t(apply(tps_L_missing, 1, reflect_point, origin = origin, pc3 = pc3))
        }
        if (length(only_R_idx) > 0) {
          tps_R_missing <- Morpho::tps3d(refmat = R_common, tarmat = R_after_common,
                                         x = R_pts[only_R_idx, , drop = FALSE])
          R_final[only_R_idx, ] <- tps_R_missing
          L_final[only_R_idx, ] <- t(apply(tps_R_missing, 1, reflect_point, origin = origin, pc3 = pc3))
        }

        # --- Save results (restore R order)
        landmarks[[specimen_id]][[L_name]] <- L_final
        landmarks[[specimen_id]][[R_name]] <- R_final[nrow(R_final):1, , drop = FALSE]

        landmarks[[specimen_id]][["V"]] <- ventral_points # save smooth V curve

        # --- Update metadata
        update_missing_info <- function(curve_name, new_curve) {
          row_idx <- which(curve_info$Specimen == specimen_id & curve_info$Curve == curve_name)
          if (length(row_idx) == 0) return(NULL)
          missing_idx <- which(apply(new_curve, 1, function(x) any(is.na(x))))
          # curve_info$missing_ldm_updated[[row_idx]] <<- list(missing_idx)
          curve_info$missing_ldm_updated[[row_idx]] <<- missing_idx
        }
        update_missing_info(L_name, L_final)
        update_missing_info(R_name, R_final)
      }

      # ============================================================
      # Case B: only one side exists -> simple reflection
      # ============================================================
      if (xor(L_exists, R_exists)) {

        existing_name <- if (L_exists) L_name else R_name
        new_name <- if (L_exists) R_name else L_name
        existing_curve <- as.matrix(specimen_data[[existing_name]])
        n_ldm <- nrow(existing_curve)

        # Find median ldm
        if (L_exists) median_ldm <- specimen_data[[L_name]][1,] else median_ldm <- specimen_data[[R_name]][n_ldm,]

        # Find local plane
        plane <- compute_local_plane(V = ventral_points, median_ldm, origin = median_ldm,
                                     weighting, weight_power)
        origin <- plane$origin
        pc3 <- plane$pc3

        # --- proceed with symmetry operations using pc3 & origin_local  ---

        reflected_curve <- t(apply(existing_curve, 1, reflect_point, origin = origin, pc3 = pc3))
        reflected_curve <- reflected_curve[nrow(reflected_curve):1, ]
        landmarks[[specimen_id]][[new_name]] <- as.matrix(reflected_curve)
        landmarks[[specimen_id]][[existing_name]] <- as.matrix(existing_curve)

        # Update curve_info with new entry
        new_row <- data.frame(
          Specimen = specimen_id,
          Curve = new_name,
          Type = ifelse(L_exists, "R", "L"),
          missing_ldm = I(list(integer(0))),
          Symmetrized = TRUE,
          stringsAsFactors = FALSE
        )
        missing_cols <- setdiff(names(curve_info), names(new_row))
        for (col in missing_cols) new_row[[col]] <- NA
        new_row <- new_row[, names(curve_info)]
        curve_info <- rbind(curve_info, new_row)
      }
    }

    # ============================================================
    # Case C: umbilical curves (U / Ubis)
    # ============================================================
    U_curves <- names(specimen_data)[grepl("^U", names(specimen_data))]
    if (length(U_curves) == 1) {
      curve_name <- U_curves
      n_ldm <- nrow(specimen_data[[curve_name]])
      global_plan <- prcomp(ventral_points)
      pc3_global <- global_plan$rotation[, 3]
      origin_global <- colMeans(ventral_points)
      reflected_curve <- t(apply(specimen_data[[curve_name]], 1, reflect_point, origin = origin_global, pc3 = pc3_global))
      # reflected_curve <- reflected_curve[nrow(reflected_curve):1, ]
      reflected_name <- if (curve_name == "U") "Ubis" else "U"
      landmarks[[specimen_id]][[reflected_name]] <- as.matrix(reflected_curve)
      landmarks[[specimen_id]][[curve_name]] <- as.matrix(specimen_data[[curve_name]])
    }else{
      curve_name1 <- U_curves[1]
      curve_name2 <- U_curves[2]
      landmarks[[specimen_id]][[curve_name1]] <- as.matrix(specimen_data[[curve_name1]])
      landmarks[[specimen_id]][[curve_name2]] <- as.matrix(specimen_data[[curve_name2]])
    }
  }

  # Return updated structure
  return(list(curve_info = curve_info, landmarks = landmarks))
}


# --------------------------------------------------
#  Function: merge_sides
# --------------------------------------------------

#' Merges left and right aperture sides into complete 3D peristome curves
#'
#' This function merges the left (`L*`) and right (`R*`) aperture curves of each specimen
#' into single, continuous apertures. The right curve is placed first, followed by the left
#' curve (with its first landmark removed to avoid duplication at the junction).
#'
#' Neutral curves (such as "U", "Ubis", or "V") that are not labeled as left or right are
#' preserved unchanged.
#'
#' @param data_sym A list containing:
#'   \describe{
#'     \item{\code{landmarks}}{A 3-level nested list of 3D landmark coordinates:
#'     \code{[[specimen]][[curve]][[landmark, XYZ]]}.}
#'     \item{\code{curve_info}}{A data.frame containing metadata for each curve.}
#'     Needed data structure is typically from \link{symmetrize_sides} outputs.
#'   }
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{landmarks}}{A 3-level list where left/right aperture curves have been merged.}
#'     \item{\code{curve_info}}{The original curve metadata, unchanged.}
#'   }
#'
#' @details
#' The function identifies matching pairs of left and right aperture curves
#' (e.g., "L3" and "R3") and concatenates their 3D coordinates into one continuous
#' aperture curve. This merging ensures each aperture is represented by a single,
#' complete set of landmarks, simplifying subsequent morphometric analyses.
#'
#' The left curve's first landmark is removed before merging to prevent duplication
#' at the ventral junction, where left and right sides meet.
#'
#' Neutral curves (those not starting with "L" or "R") are copied as-is.
#'
#' @examples
#' \dontrun{
#' # Load dataset
#' data("example_spec_tab", package = "AmmoniTools")
#' data("example_data", package = "AmmoniTools")
#'
#' # Resample curve with 25 landmarks for each curve type
#' resampled_data <- resample_curves(example_data, n_landmarks = list(R = 25, L = 25, U = 25, V = 25))
#'
#' # Reorder curves
#' reordered_data <- reorder_landmarks(resampled_data)
#'
#' # Symmetrize sides
#' data_sym <- symmetrize_sides(reordered_data,
#'                              weighting = "exp",
#'                              weight_power  = 5,
#'                              k_spline = 7)
#'
#' # Merge LR sides
#' merged_data <- merge_sides(data_sym)
#' }
#'
#' @export

merge_sides <- function(data_sym) {

  # --------------------------------------------------
  # Step 1 - Extract components
  # --------------------------------------------------
  landmarks <- data_sym$landmarks
  curve_info <- data_sym$curve_info

  combined_landmarks <- list()

  # --------------------------------------------------
  # Step 2 - Iterate over each specimen
  # --------------------------------------------------
  for (specimen_id in names(landmarks)) {

    specimen_data <- landmarks[[specimen_id]]

    # Identify left and right curves by prefix (e.g., "L3", "R3")
    left_landmarks <- grep("^L\\d+$", names(specimen_data), value = TRUE)
    right_landmarks <- grep("^R\\d+$", names(specimen_data), value = TRUE)

    # Identify neutral (non-left/right) curves
    all_landmarks <- names(specimen_data)
    neutral_landmarks <- setdiff(all_landmarks, c(left_landmarks, right_landmarks))

    # Extract numeric IDs shared by both left and right curves
    left_numbers <- sub("^L", "", left_landmarks)
    right_numbers <- sub("^R", "", right_landmarks)
    common_numbers <- intersect(left_numbers, right_numbers)

    combined_specimen <- list()

    # --------------------------------------------------
    # Step 3 - Merge matching left/right curves
    # --------------------------------------------------
    for (num in common_numbers) {

      left_name <- paste0("L", num)
      right_name <- paste0("R", num)

      # Retrieve the coordinate matrices
      left_data <- specimen_data[[left_name]]
      right_data <- specimen_data[[right_name]]

      # Remove the first point of the left side to avoid duplication at the junction
      if (nrow(left_data) > 1) {
        left_data <- left_data[-1, , drop = FALSE]
      }

      # Concatenate right + left to form a single aperture curve
      combined_curve <- rbind(right_data, left_data)

      # Store under numeric aperture ID
      combined_specimen[[num]] <- combined_curve
    }

    # --------------------------------------------------
    # Step 4 - Add neutral curves unchanged
    # --------------------------------------------------
    for (neutral_name in neutral_landmarks) {
      combined_specimen[[neutral_name]] <- specimen_data[[neutral_name]]
    }

    # --------------------------------------------------
    # Step 5 - Store results
    # --------------------------------------------------
    combined_landmarks[[specimen_id]] <- combined_specimen
  }

  # --------------------------------------------------
  # Step 6 - Return updated structure
  # --------------------------------------------------
  return(list(
    landmarks = combined_landmarks,
    curve_info = curve_info
  ))
}


# --------------------------------------------------
#  Function: impute_weighted_TPS
# --------------------------------------------------

#' Estimates missing landmark coordinates using distance-weighted TPS interpolation, prioritizing geometrically, ontogenetically and/or taxonomically closest reference apertures
#'
#' \strong{Internal function – not intended to be called directly by users.}
#'
#' This function reconstructs missing 3D landmarks on incomplete ammonite
#' aperture curves by borrowing geometric information from complete curves
#' (donors) using Thin Plate Spline (TPS) interpolation.
#'
#' Donor curves are selected following a hierarchical similarity scheme
#' (e.g. same specimen, same species, same genus) and are weighted according
#' to their Procrustes distance to the incomplete curve. The final reconstruction
#' is obtained as a weighted average of TPS predictions from the selected donors.
#'
#' @section Intended use:
#' This function is a low-level internal routine used exclusively by
#' \link{reconstruct_missingldm} as part of the missing landmark reconstruction
#' pipeline. Users should normally not call this function directly, as it assumes
#' a specific data structure and preprocessing workflow handled upstream.
#'
#' @param complete Named list of complete curves (numeric matrices:
#' \eqn{n \times 3}). Names must follow the format
#' \code{"SpecimenID:CurveName"}.
#' @param incomplete Named list of incomplete curves with missing landmarks
#' encoded as \code{NA}. Naming format must match \code{complete}.
#' @param metadata A data.frame or tibble containing specimen-level metadata.
#' Must include a \code{Specimen} column and may include higher-level grouping
#' variables such as \code{Species} or \code{Genus}.
#' @param hierarchy List defining the order of donor prioritization
#' (e.g. \code{list("Specimen", "Species", "Genus")}).
#' Donors are first searched at the most specific level, then progressively
#' relaxed if needed.
#' @param max_donors Integer. Maximum number of donor curves used to compute
#' the weighted reconstruction.
#' @param weighting_function Character string specifying how donor weights
#' are derived from Procrustes distances:
#' \describe{
#'   \item{"inverse"}{Inverse-distance weighting
#'   (\eqn{w = 1 / d^{p}}).}
#'   \item{"gaussian"}{Gaussian kernel weighting
#'   (\eqn{w = \exp(-d^2 / (2\sigma^2))}).}
#' }
#' @param weight_power Numeric. Power exponent used for inverse-distance
#' weighting.
#' @param sigma Numeric. Standard deviation of the Gaussian kernel.
#' If \code{NULL}, it is estimated from the donor distance distribution.
#' @param verbose Logical. If \code{TRUE}, prints progress messages during
#' imputation.
#' @param allow_fallback Logical. If TRUE (default), when no donor curves
#' are found at any hierarchical level, all remaining complete curves are used
#' as fallback donors. If FALSE, no fallback is applied and the incomplete
#' curve is skipped (not reconstructed).
#'
#' @return
#' A list with two components:
#' \describe{
#'   \item{\code{coords}}{Named list of reconstructed curves, including all
#'   original complete curves and successfully imputed incomplete ones.
#'   Incomplete curves for which no valid donors were found (when
#'   \code{allow_fallback = FALSE}) are omitted.}
#'   \item{\code{provenance}}{List of data frames documenting, for each imputed
#'   curve, the donor curves used and their associated weights. Curves that
#'   were skipped do not appear in this list.}
#' }
#'
#' @details
#' For each incomplete curve, the algorithm proceeds as follows:
#' \enumerate{
#'   \item Identify landmarks that are present (non-missing).
#'   \item Select donor curves following the specified hierarchy.
#'   \item If no donors are found:
#'     \itemize{
#'       \item If \code{allow_fallback = TRUE}, all remaining complete curves
#'       are used as fallback donors.
#'       \item If \code{allow_fallback = FALSE}, the curve is skipped and not reconstructed.
#'     }
#'   \item Align incomplete and donor curves using Generalized Procrustes Analysis.
#'   \item Compute Procrustes distances between the incomplete curve and each donor.
#'   \item Retain the closest donors (up to \code{max_donors}).
#'   \item Estimate missing landmarks using TPS transformations from each donor.
#'   \item Combine donor predictions using distance-based weights.
#' }
#'
#' Curves with fewer than three observed landmarks are skipped, as TPS
#' estimation is not geometrically defined in that case.
#'
#' @seealso \link{reconstruct_missingldm}
#' @importFrom Morpho tps3d computeTransform applyTransform
#' @importFrom geomorph gpagen
#' @importFrom dplyr left_join filter bind_rows
#' @importFrom tibble tibble
#' @keywords internal
#' @export

impute_weighted_TPS <- function(complete, incomplete, metadata, hierarchy = list(),
                                max_donors = 3,
                                weighting_function = c("inverse", "gaussian"),
                                weight_power = 3, sigma = NULL,
                                allow_fallback = TRUE,
                                verbose = FALSE) {

  # ---------------------------
  # Step 0 - Setup
  # ---------------------------
  weighting_function <- match.arg(weighting_function)
  res <- list()
  provenance <- list()

  # Exclude curves unrelated to lateral sides
  complete <- complete[!grepl(":[UVR]|:Ubis", names(complete))]
  incomplete <- incomplete[!grepl(":[UVR]|:Ubis", names(incomplete))]

  # Helper to extract specimen/curve IDs
  extract_ids <- function(nm) {
    parts <- strsplit(nm, ":", fixed = TRUE)[[1]]
    list(specimen = as.integer(parts[1]), curve = parts[2])
  }

  # Annotate metadata for complete curves
  complete_meta <- tibble(
    name = names(complete),
    Specimen = vapply(names(complete), function(nm) extract_ids(nm)$specimen, integer(1)),
    curve = vapply(names(complete), function(nm) extract_ids(nm)$curve, character(1))
  ) %>%
    left_join(metadata, by = c("Specimen" = "Specimen"))

  # ---------------------------
  # Step 1 - Copy complete curves directly
  # ---------------------------
  for (nm in names(complete)) {
    res[[nm]] <- complete[[nm]]
    provenance[[nm]] <- tibble(donor = NA_character_, weight = NA_real_)
  }

  # ---------------------------
  # Step 2 - Impute each incomplete curve
  # ---------------------------
  for (nm in names(incomplete)) {
    ids <- extract_ids(nm)
    target_specimen <- ids$specimen
    curve_name <- ids$curve
    incomplete_curve <- incomplete[[nm]]

    if (is.null(incomplete_curve)) next

    # Retrieve metadata for target specimen
    target_meta <- metadata %>% filter(Specimen == target_specimen)
    if (nrow(target_meta) == 0) next

    priority_set <- tibble()
    used_names <- character(0)

    # ---------------------------
    # Step 2a - Hierarchical donor selection
    # ---------------------------
    for (level in hierarchy) {
      if (is.character(level)) level <- list(level)
      filters <- rep(TRUE, nrow(complete_meta))

      for (var in level) {
        # Skip invalid metadata columns
        if (!var %in% names(complete_meta) || !var %in% names(target_meta)) next

        value <- target_meta[[var]][1]
        if (is.na(value)) {
          filters <- rep(FALSE, nrow(complete_meta))
        } else {
          filters <- filters & (complete_meta[[var]] == value)
        }

        # Select donors not already used
        if (any(filters)) {
          mask <- filters & !(complete_meta$name %in% used_names)
          selected <- complete_meta[mask, ]
          priority_set <- bind_rows(priority_set, selected)
          used_names <- c(used_names, selected$name)
          break
        }
      }
    }

    # Fallback if no donors found
    if (nrow(priority_set) == 0) {

      if (!allow_fallback) {
        if (verbose)
          message("No donors found for ", nm, " — curve removed (no fallback allowed).")
        next
      }

      fallback <- complete_meta %>% filter(!(name %in% used_names))
      priority_set <- bind_rows(priority_set, fallback)

      if (verbose)
        message("Fallback used for ", nm)
    }

    # ---------------------------
    # Step 2b - Identify present landmarks
    # ---------------------------
    missing_rows <- apply(incomplete_curve, 1, function(x) any(is.na(x)))
    present_idx <- which(!missing_rows)
    if (length(present_idx) < 3) {
      warning(paste("Too few landmarks to estimate for", nm))
      next
    }

    inc_partial <- as.matrix(incomplete_curve[present_idx, , drop = FALSE])
    donor_names <- priority_set$name
    donor_partials <- lapply(donor_names, function(dn) complete[[dn]][present_idx, , drop = FALSE])

    # ---------------------------
    # Step 2c - GPA alignment + distance computation
    # ---------------------------
    all_curves <- abind::abind(c(donor_partials, list(inc_partial)), along = 3)
    gpa <- geomorph::gpagen(all_curves, print.progress = FALSE)

    n_donors <- length(donor_names)
    distances <- vapply(1:n_donors, function(i) {
      sqrt(sum((gpa$coords[, , i] - gpa$coords[, , n_donors + 1])^2))
    }, numeric(1))

    # ---------------------------
    # Step 2d - Select top donors
    # ---------------------------
    top_idx <- order(distances)[1:min(max_donors, length(distances))]
    donor_names <- donor_names[top_idx]
    distances <- distances[top_idx]

    # ---------------------------
    # Step 2e - Weighted TPS interpolation
    # ---------------------------
    weighted_estimates <- array(0, dim = c(nrow(incomplete_curve), 3, length(top_idx)))
    weights <- numeric(length(top_idx))

    for (j in seq_along(top_idx)) {
      donor_curve <- as.matrix(complete[[donor_names[j]]])
      donor_partial <- as.matrix(donor_curve[present_idx, , drop = FALSE])
      donor_full <- donor_curve

      dist <- distances[j]
      weight <- if (weighting_function == "inverse") {
        1 / (dist^weight_power + 1e-6)
      } else {
        if (is.null(sigma)) sigma <- mean(distances)
        exp(-(dist^2) / (2 * sigma^2))
      }
      colnames(inc_partial) <- colnames(donor_partial)
      #rownames(inc_partial) <- NULL
      inc_partial <- as.matrix(inc_partial)

      coeff <- Morpho::computeTransform(x = inc_partial, y = donor_partial, type = "tps")
      tps_coords <- Morpho::applyTransform(donor_full, coeff)

      weighted_estimates[, , j] <- tps_coords
      weights[j] <- weight
    }

    weights <- weights / sum(weights)
    estimated_curve <- apply(weighted_estimates, c(1, 2), function(x) sum(x * weights))

    # ---------------------------
    # Step 2f - Store results
    # ---------------------------
    res[[nm]] <- estimated_curve
    provenance[[nm]] <- tibble(donor = donor_names, weight = weights)

    if (verbose)
      message("Imputed: ", nm, " with ", paste(donor_names, collapse = ", "))
  }

  # ---------------------------
  # Step 3 - Return results
  # ---------------------------
  list(
    coords = res,
    provenance = provenance
  )
}


# --------------------------------------------------
#  Function: reconstruct_missingldm
# --------------------------------------------------
#' Reconstructs incomplete apertures curves using TPS-based interperistome reconstruction
#'
#' This function detects incomplete aperture curves in a dataset of ammonite landmarks,
#' imputes missing landmarks using hierarchical weighted Thin Plate Spline (TPS)
#' interpolation, and replaces the missing curves in the global data structure.
#' The imputation is performed using the {\link{impute_weighted_TPS}} function,
#' which borrows information from complete curves of other specimens according to
#' a hierarchical donor selection scheme (e.g., same specimen -> same species -> same genus).
#'
#' @param merged_data A list containing:
#'   \describe{
#'     \item{\code{landmarks}}{A nested list of 3D landmark coordinates
#'     (one matrix per curve per specimen).}
#'     \item{\code{curve_info}}{A data.frame describing each curve
#'     (used for downstream processing).}
#'     Needed data structure is typically from \link{merge_sides} outputs.
#'   }
#' @param metadata A data.frame or tibble containing at least one column \code{Specimen},
#'   and optionally higher-level attributes such as \code{Species} or \code{Genus}
#'   to guide hierarchical donor selection.
#' @param hierarchy A list defining the hierarchical levels for donor prioritization.
#'   For example: \code{list("Specimen", "Species", "Genus")}.
#' @param weight_power Numeric. Exponent used in inverse distance weighting (default = 3).
#' @param max_donors Integer. Maximum number of donor curves to use for each imputation (default = 3).
#' @param weighting_function Character. Type of weighting function for donor distances:
#'   \describe{
#'     \item{\code{"inverse"}}{Weights = 1 / (distance^weight_power)}
#'     \item{\code{"gaussian"}}{Weights = exp(-(distance^2)/(2*sigma^2))}
#'   }
#'   Default: \code{"gaussian"}.
#' @param verbose Logical. If TRUE, prints detailed messages during the imputation process.
#' @param sigma Numeric. Standard deviation parameter for Gaussian weighting.
#'   If NULL, estimated automatically from donor distances.
#' @param allow_fallback Logical. Passed to \link{impute_weighted_TPS}.
#' If TRUE (default), fallback donors are used when no hierarchical match is found.
#' If FALSE, incomplete curves without valid donors are removed from the dataset.
#'
#' @return A list with three components:
#'   \describe{
#'     \item{\code{curve_info}}{
#'       Updated curve information (identical to input).
#'     }
#'     \item{\code{landmarks}}{
#'       Updated nested landmark list including reconstructed curves.
#'       Curves that could not be reconstructed (e.g., no donors available when
#'       \code{allow_fallback = FALSE}, or still containing \code{NA} after TPS
#'       interpolation) are removed from the dataset.
#'     }
#'     \item{\code{prov_reconstructed_ldm}}{
#'       List of data frames documenting, for each successfully imputed
#'       curve, the donor curves used and their associated weights.
#'     }
#'   }
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Separate complete and incomplete curves within each specimen
#'         (ventral and umbilical curves are ignored).
#'   \item Impute missing landmarks on incomplete curves using
#'         \link{impute_weighted_TPS} with hierarchical donor selection.
#'   \item Remove all originally incomplete curves from the dataset.
#'   \item Reinsert only successfully reconstructed curves.
#'   \item Remove any reconstructed curves that still contain missing
#'         coordinates after TPS interpolation.
#'   \item Return a fully reconstructed dataset ready for geometric analysis.
#' }
#'
#' This guarantees that the returned dataset contains only geometrically
#' complete aperture curves.
#'
#'
#' @examples
#' \dontrun{
#' # Load dataset
#' data("example_spec_tab", package = "AmmoniTools")
#' data("example_data", package = "AmmoniTools")
#'
#' # Resample curve with 25 landmarks for each curve type
#' resampled_data <- resample_curves(example_data, n_landmarks = list(R = 25, L = 25, U = 25, V = 25))
#'
#' # Reorder curves
#' reordered_data <- reorder_landmarks(resampled_data)
#'
#' # Symmetrize sides
#' data_sym <- symmetrize_sides(reordered_data,
#'                              weighting = "exp",
#'                              weight_power  = 5,
#'                              k_spline = 7)
#'
#' # Merge LR sides
#' merged_data <- merge_sides(data_sym)
#'
#' # Reconstruct missing landmarks using TPS-base interperistome reconstruction
#' reconstructed_data <- reconstruct_missingldm(
#'   merged_data,
#'   metadata = spec_tab,
#'   hierarchy = list("Specimen", "Species", "Genus"),
#'   weighting_function = "gaussian",
#'   max_donors = 3,
#'   verbose = TRUE,
#'   allow_fallback = TRUE
#' )
#'
#' }
#'
#' @seealso \link{impute_weighted_TPS}
#' @export
#'
#' @importFrom dplyr left_join filter bind_rows
#' @importFrom tibble tibble
#' @importFrom Morpho computeTransform applyTransform
#' @importFrom geomorph gpagen

reconstruct_missingldm <- function(merged_data,
                                   metadata,
                                   hierarchy = list(),
                                   weight_power = 1,
                                   max_donors = 5,
                                   weighting_function = "gaussian",
                                   verbose = TRUE,
                                   sigma = NULL,
                                   allow_fallback = TRUE) {

  # --------------------------------------------------
  # Step 1 - Separate complete and incomplete curves
  # --------------------------------------------------
  # This internal helper identifies which aperture curves contain missing landmarks (NA)
  # and splits them into two lists: complete and incomplete.
  separate_complete_incomplete <- function(data) {
    complete <- list()
    incomplete <- list()

    for (specimen_id in names(data)) {
      for (curve_name in names(data[[specimen_id]])) {

        # Skip auxiliary curves (U, Ubis, V) that should not be imputed
        if (curve_name %in% c("U", "Ubis", "V")) next

        mat <- data[[specimen_id]][[curve_name]]

        # Remove column names / dimnames
        if (!is.null(dimnames(mat))) attr(mat, "dimnames") <- NULL

        key <- paste0(specimen_id, ":", curve_name)

        # --- Detect NA and classify ---
        if (any(is.na(mat))) {
          incomplete[[key]] <- mat
        } else {
          complete[[key]] <- mat
        }
      }
    }

    list(complete = complete, incomplete = incomplete)
  }

  # --- Apply separation function
  split_data <- separate_complete_incomplete(merged_data$landmarks)

  complete <- split_data$complete
  incomplete <- split_data$incomplete

  # --------------------------------------------------
  # Step 2 - Impute missing landmarks using weighted TPS
  # --------------------------------------------------
  data_withoutmissing <- impute_weighted_TPS(
    complete,
    incomplete,
    metadata = metadata,
    hierarchy = hierarchy,
    weight_power = weight_power,
    max_donors = max_donors,
    weighting_function = weighting_function,
    sigma = sigma,
    allow_fallback = allow_fallback,
    verbose = verbose
  )

  # Make a copy of the data object to modify safely
  data_sym_original <- merged_data

  # --------------------------------------------------
  # Step 3 - Clean and replace reconstructed curves
  # --------------------------------------------------

  # 1) Remove all originally incomplete curves
  for (nm in names(incomplete)) {

    specimen_curve <- strsplit(nm, ":")[[1]]
    specimen_id <- specimen_curve[1]
    curve_name <- specimen_curve[2]

    if (!is.null(data_sym_original$landmarks[[specimen_id]][[curve_name]])) {
      data_sym_original$landmarks[[specimen_id]][[curve_name]] <- NULL

      if (verbose)
        message("Removed original incomplete curve: ", nm)
    }
  }

  # 2) Reinsert only successfully reconstructed curves
  for (nm in names(data_withoutmissing$coords)) {

    specimen_curve <- strsplit(nm, ":")[[1]]
    specimen_id <- specimen_curve[1]
    curve_name <- specimen_curve[2]

    coords <- data_withoutmissing$coords[[nm]]

    # Skip if still contains NA (extra safety)
    if (any(is.na(coords))) {
      if (verbose)
        message("Skipped (still NA after imputation): ", nm)
      next
    }

    full_df <- as.data.frame(coords)
    colnames(full_df) <- c("X", "Y", "Z")

    if (is.null(data_sym_original$landmarks[[specimen_id]])) {
      data_sym_original$landmarks[[specimen_id]] <- list()
    }

    data_sym_original$landmarks[[specimen_id]][[curve_name]] <- full_df

    if (verbose)
      message("Inserted reconstructed curve: ", nm)
  }

  # --------------------------------------------------
  # Step 4 - Return reconstructed dataset
  # --------------------------------------------------
  reconstructed_data <- list(
    curve_info = data_sym_original$curve_info,
    landmarks  = data_sym_original$landmarks,
    prov_reconstructed_ldm = data_withoutmissing$provenance
  )

  return(reconstructed_data)
}



