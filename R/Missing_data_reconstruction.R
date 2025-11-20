
# --------------------------------------------------
#  MISSING DATA RECONSTRUCTION
#
# symmetrize_sides : Symmetrize and reconstruct left/right sides and umbilical aperture curves
# merge_sides : Merge 3D coordinates of left and right sides of each aperture
# test_symdist : Compute mean distances for symmetry reconstruction (used as an error metric)
# impute_weighted_TPS : Impute missing landmarks using weighted TPS interpolation
# reconstruct_missingldm : Reconstruct missing aperture landmarks by hierarchical weighted TPS
# --------------------------------------------------

utils::globalVariables(c("Specimen", "name", "spec_tab"))

# --------------------------------------------------
#  Function: symmetrize_sides
# --------------------------------------------------

#' Symmetrize and reconstruct left/right sides and umbilical aperture curves
#'
#' This function reconstructs missing left/right (L/R) or umbilical (U/Ubis) curves
#' of an ammonite aperture by reflecting existing curves across a symmetry plane.
#' When both sides exist, it performs a full geometric symmetrization using thin-plate
#' spline interpolation between corresponding landmarks. The symmetry plane is
#' estimated from ventral landmarks of the same specimen, optionally weighted by their
#' proximity to the ventral edge.
#'
#' @param data A list containing:
#'   \describe{
#'     \item{\code{landmarks}}{A nested list of 3D landmarks for each curve of each specimen.}
#'     \item{\code{curve_info}}{A data.frame describing curves (used to detect missing sides).}
#'   }
#' @param weighting Character. Method used to weight ventral landmarks when computing the
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
#' data_sym <- symmetrize_sides(data,
#'                              weighting = "exp",
#'                              weight_power  = 1,
#'                              k_spline = 7)
#' }
#'
#' @importFrom stats prcomp
#' @importFrom Morpho tps3d
#' @export


symmetrize_sides <- function(data, weighting = c("none", "linear", "exp"),
                             weight_power  = 1, k_spline = 7) {

  weighting <- match.arg(weighting)
  curve_info <- data$curve_info
  landmarks <- data$landmarks

  # ------------------------------------------------------------------
  # Helper function: reflect a point across a plane
  # ------------------------------------------------------------------
  # The plane is defined by a normal vector `pc3` and passes through `origin`.
  reflect_point <- function(X, origin, pc3) {
    d <- sum((X - origin) * pc3)       # signed distance from point to plane
    X - 2 * d * pc3                    # reflected point across plane
  }

  # ------------------------------------------------------------------
  # Helper: smooth a 3D curve using P-splines (mgcv::gam)
  # ------------------------------------------------------------------
  pspline_smooth <- function(points, n_out = 100, k = k_spline) {
    t <- seq(0, 1, length.out = nrow(points))
    t_out <- seq(0, 1, length.out = n_out)
    smoothed <- sapply(1:3, function(d) {
      fit <- mgcv::gam(points[, d] ~ s(t, bs = "ps", k = k))
      predict(fit, newdata = data.frame(t = t_out))
    })
    matrix(smoothed, ncol = 3, dimnames = list(NULL, c("X", "Y", "Z")))
  }

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

  # ==================================================================
  # Main loop: process all specimens
  # ==================================================================
  for (specimen_id in names(landmarks)) {
    specimen_data <- landmarks[[specimen_id]]
    if (!("V" %in% names(specimen_data))) next

    # ------------------------------------------------------------
    # 1. Compute local ventral plane
    # ------------------------------------------------------------
    ventral_points <- pspline_smooth(as.matrix(specimen_data$V), n_out = length(specimen_data$V), k = k_spline)
    landmarks[[specimen_id]][["V"]] <- as.matrix(specimen_data$V)

    origin <- colMeans(ventral_points)
    centered <- ventral_points - matrix(origin, nrow(ventral_points), 3, byrow = TRUE)
    dist_curve <- c(0, cumsum(sqrt(rowSums(diff(ventral_points)^2))))

    # Weighting scheme along the ventral curve
    w <- switch(weighting,
                none   = rep(1, length(dist_curve)),
                linear = pmax(1 - (weight_power  * dist_curve / max(dist_curve)), 0),
                exp    = exp(-weight_power  * dist_curve / max(dist_curve))
    )
    W <- diag(w / sum(w))
    cov_w <- t(centered) %*% W %*% centered
    pc3 <- eigen(cov_w)$vectors[, 3]  # third PC defines the plane's normal

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
      reflected_curve <- reflected_curve[nrow(reflected_curve):1, ]
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

#' Merge 3D coordinates of left and right sides of each aperture
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
#' merged_data <- merge_sides(data_sym)
#' str(merged_data$landmarks)
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

#' Impute missing landmarks using weighted TPS interpolation
#'
#' This function estimates missing 3D landmarks of incomplete ammonite aperture curves
#' by borrowing information from complete curves (donors) using a Thin Plate Spline (TPS)
#' interpolation. Donor curves are selected based on an optionally hierarchical similarity scheme
#' (e.g., same specimen -> same species -> same genus) and weighted according to their
#' Procrustes distance from the incomplete curve.
#'
#' @param complete A named list of complete curves (3D matrices: n_landmarks * 3).
#'   Names must follow the format `"SpecimenID:CurveName"`.
#' @param incomplete A named list of incomplete curves with missing landmarks (NA).
#'   Names must follow the same format as \code{complete}.
#' @param metadata A data.frame or tibble containing at least a column `Specimen`
#'   and optionally higher-level descriptors (e.g., `Species`, `Genus`).
#' @param hierarchy A list defining the order of donor prioritization, e.g.:
#'   \code{list("Specimen", "Species", "Genus")}.
#'   The function first looks for donors with the same specimen, then same species, etc.
#' @param max_donors Integer. Maximum number of donor curves to use for weighted averaging.
#'   Default: \code{3}.
#' @param weighting_function Character. Defines how weights are computed from Procrustes distances:
#'   \itemize{
#'     \item \code{"inverse"}: weights = 1 / (distance^power)
#'     \item \code{"gaussian"}: weights = exp(-(distance^2)/(2*sigma^2))
#'   }
#'   Default: \code{"inverse"}.
#' @param weight_power Numeric. Power exponent used in inverse weighting (default = 3).
#' @param sigma Numeric. Standard deviation parameter for Gaussian weighting
#'   (if NULL, estimated as the mean distance among selected donors).
#' @param verbose Logical. If TRUE, prints progress messages for each imputed curve.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{coords}}{Named list of imputed curves (same format as input).}
#'   \item{\code{provenance}}{List of data.frames detailing donor curves and their weights for each imputation.}
#' }
#'
#' @details
#' The algorithm proceeds as follows for each incomplete curve:
#' \enumerate{
#'   \item Parse specimen and curve identifiers.
#'   \item Identify donor curves from the complete dataset following the hierarchy.
#'   \item Superimpose incomplete and donor curves using Generalized Procrustes Analysis (GPA).
#'   \item Compute Procrustes distances to all donors.
#'   \item Select the best matching donors (up to \code{max_donors}).
#'   \item Interpolate missing landmarks using TPS transformations from each donor.
#'   \item Compute a weighted mean of all donor TPS predictions.
#' }
#'
#' Curves with too few known landmarks (<3) are skipped with a warning.
#'
#' @examples
#' \dontrun{
#' reconstructed_ldm <- impute_weighted_TPS(complete = complete_curves,
#'                                          incomplete = incomplete_curves,
#'                                          metadata = spec_tab,
#'                                          hierarchy = list("Species", "Genus"),
#'                                          weighting_function = "gaussian",
#'                                          max_donors = 5)
#' }
#'
#' @importFrom Morpho tps3d
#' @importFrom geomorph gpagen
#' @importFrom dplyr left_join filter bind_rows
#' @importFrom tibble tibble
#' @export

impute_weighted_TPS <- function(complete, incomplete, metadata, hierarchy = list(),
                                max_donors = 3,
                                weighting_function = c("inverse", "gaussian"),
                                weight_power = 3, sigma = NULL,
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
      fallback <- complete_meta %>% filter(!(name %in% used_names))
      priority_set <- bind_rows(priority_set, fallback)
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
      message("Imputed: ", nm, " with ", length(donor_names), " donors.")
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
#' Reconstruct missing aperture landmarks by hierarchical weighted TPS
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
#'
#' @return A list with two components:
#'   \describe{
#'     \item{\code{curve_info}}{Updated curve information (identical to input for now).}
#'     \item{\code{landmarks}}{Updated nested landmark list including reconstructed curves.}
#'   }
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item **Separate** complete and incomplete curves within each specimen
#'         (ventral and umbilical curves are ignored).
#'   \item **Impute** missing landmarks on incomplete curves using
#'         \link{impute_weighted_TPS} with hierarchical donor selection.
#'   \item **Replace** the reconstructed curves in the global dataset,
#'         removing any outdated left/right versions to avoid duplication.
#'   \item **Return** a fully reconstructed dataset ready for geometric analysis.
#' }
#'
#' This function is intended to run after the symmetrization step, ensuring
#' that missing aperture sides are geometrically reconstructed while preserving
#' the biological structure of the dataset.
#'
#'
#' @examples
#' \dontrun{
#' reconstructed_data <- reconstruct_missingldm(
#'   merged_data = data_sym,
#'   metadata = spec_tab,
#'   hierarchy = list("Specimen", "Species", "Genus"),
#'   weighting_function = "gaussian",
#'   max_donors = 3,
#'   verbose = TRUE
#' )
#' }
#'
#' @seealso \link{impute_weighted_TPS}
#' @export
#'
#' @importFrom dplyr left_join filter bind_rows
#' @importFrom tibble tibble
#' @importFrom Morpho computeTransform applyTransform
#' @importFrom geomorph gpagen
#'
reconstruct_missingldm <- function(merged_data,
                                   metadata,
                                   hierarchy = list(),
                                   weight_power = 1,
                                   max_donors = 5,
                                   weighting_function = "gaussian",
                                   verbose = TRUE,
                                   sigma = NULL) {

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
    sigma = NULL,
    verbose = verbose
  )

  # Make a copy of the data object to modify safely
  data_sym_original <- merged_data

  # --------------------------------------------------
  # Step 3 - Replace reconstructed curves in global data
  # --------------------------------------------------
  for (nm in names(data_withoutmissing$coords)) {
    specimen_curve <- strsplit(nm, ":")[[1]]
    specimen_id <- specimen_curve[1]
    curve_name <- specimen_curve[2]

    coords <- data_withoutmissing$coords[[nm]]
    full_df <- as.data.frame(coords)
    colnames(full_df) <- c("X", "Y", "Z")

    # Remove any existing L/R curves for this specimen before reinsertion
    existing_names <- names(data_sym_original$landmarks[[specimen_id]])
    if (any(grepl("^[LR]", existing_names))) {
      data_sym_original$landmarks[[specimen_id]] <-
        data_sym_original$landmarks[[specimen_id]][!grepl("^[LR]", existing_names)]
    }

    # Create empty list if specimen not yet present
    if (is.null(merged_data$landmarks[[specimen_id]])) {
      data_sym_original$landmarks[[specimen_id]] <- list()
    }

    # Store reconstructed curve
    data_sym_original$landmarks[[specimen_id]][[curve_name]] <- full_df
  }

  # --------------------------------------------------
  # Step 4 - Return reconstructed dataset
  # --------------------------------------------------
  reconstructed_data <- list(
    curve_info = data_sym_original$curve_info,
    landmarks  = data_sym_original$landmarks
  )

  return(reconstructed_data)
}



