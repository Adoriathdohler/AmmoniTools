
# --------------------------------------------------
# MISSING DATA RECONSTRUCTION TEST
#
# filter_curves : Filter and extract valid paired L/R curves from landmark data
# test_symdist : Compute mean distances for symmetry reconstruction (used as an error metric)
# test_reconstruct_missingldm :
# --------------------------------------------------

utils::globalVariables(c("curve_num ", "spec_tab"))

# ---------------------------------------------------------------------
# Function: filter_curves
# ---------------------------------------------------------------------

#' @title Filter and extract valid paired L/R curves from landmark data
#'
#' @description
#' This function filters curves from a reordered landmark dataset (`reordered_data`)
#' to retain only valid *left/right* (L/R) aperture pairs for each specimen, and optionally
#' only those without missing landmarks. Umbilical (U) and ventral (V) curves are also
#' retained for contextual information but are not used for pairing.
#'
#' @param reordered_data A list containing:
#'   \itemize{
#'     \item \code{landmarks}: a nested list of 3D landmark coordinates for each specimen and curve.
#'     \item \code{curve_info}: a data frame with metadata for each curve
#'           (columns should include at least \code{Specimen}, \code{Curve}, \code{Type}, and \code{missing_ldm}).
#'   }
#' @param complete_only Logical, default = \code{FALSE}.
#'   If \code{TRUE}, only curves without missing landmarks (i.e., with an empty \code{missing_ldm}) are retained.
#'
#' @return
#' A list with the same structure as the input, containing:
#'   \describe{
#'     \item{\code{landmarks}}{Filtered landmarks for specimens that have both L and R curves.}
#'     \item{\code{curve_info}}{Subset of the corresponding metadata for retained curves.}
#'   }
#'
#' @details
#' For each specimen:
#'   \enumerate{
#'     \item All L and R curves are checked for completeness (based on \code{missing_ldm}).
#'     \item Curves of type U and V are always retained (they are excluded from pairing).
#'     \item Only L/R curves that share the same suffix (e.g., L1 and R1) are considered valid pairs.
#'     \item The function builds a new dataset with only paired L/R curves (and optionally U/V).
#'   }
#'
#' @examples
#' \dontrun{
#' filtered <- filter_curves(reordered_data, complete_only = TRUE)
#' }
#'
#' @export


filter_curves <- function(reordered_data, complete_only = FALSE) {

  filtered_landmarks <- list()  # initialize output container

  # -------------------------------------------------------------------
  # Loop over each specimen in the dataset
  # -------------------------------------------------------------------
  for (specimen_id in names(reordered_data$landmarks)) {

    specimen_curves <- reordered_data$landmarks[[specimen_id]]  # extract curves
    valid_L <- list()   # valid left-side curves
    valid_R <- list()   # valid right-side curves
    UV_curves <- list() # ventral and umbilical curves

    # -----------------------------------------------------------------
    # Identify valid L/R curves and store U/V curves separately
    # -----------------------------------------------------------------
    for (curve_id in names(specimen_curves)) {

      # find the corresponding row in curve_info
      match_row <- which(
        reordered_data$curve_info$Specimen == specimen_id &
          reordered_data$curve_info$Curve == curve_id
      )

      if (length(match_row) == 1) {
        info <- reordered_data$curve_info[match_row, ]

        # Keep U and V curves separately (for ventral guidance)
        if (info$Type %in% c("V", "U")) {
          UV_curves[[curve_id]] <- specimen_curves[[curve_id]]
        }

        # Check completeness condition (if required)
        if (!complete_only) {
          # keep all L/R curves regardless of missing landmarks
          if (info$Type == "L") valid_L[[curve_id]] <- specimen_curves[[curve_id]]
          if (info$Type == "R") valid_R[[curve_id]] <- specimen_curves[[curve_id]]
        } else {
          # keep only curves with no missing landmarks
          if (info$Type == "L" && length(info$missing_ldm[[1]]) == 0)
            valid_L[[curve_id]] <- specimen_curves[[curve_id]]
          if (info$Type == "R" && length(info$missing_ldm[[1]]) == 0)
            valid_R[[curve_id]] <- specimen_curves[[curve_id]]
        }
      }
    }

    # -----------------------------------------------------------------
    # Identify and retain only *paired* L/R curves (sharing same suffix)
    # -----------------------------------------------------------------
    common_suffixes <- intersect(
      gsub("^L", "", names(valid_L)),
      gsub("^R", "", names(valid_R))
    )

    paired_curves <- list()
    for (suffix in common_suffixes) {
      Lname <- paste0("L", suffix)
      Rname <- paste0("R", suffix)
      paired_curves[[Lname]] <- valid_L[[Lname]]
      paired_curves[[Rname]] <- valid_R[[Rname]]
    }

    # -----------------------------------------------------------------
    # Store filtered curves and metadata if pairs exist
    # -----------------------------------------------------------------
    if (length(paired_curves) > 0) {

      # combine paired L/R curves and U/V reference curves
      filtered_landmarks$landmarks[[specimen_id]] <- c(paired_curves, UV_curves)

      # identify the metadata rows to retain
      kept_curves <- c(names(paired_curves), names(UV_curves))
      rows_to_keep <- reordered_data$curve_info$Specimen == specimen_id &
        reordered_data$curve_info$Curve %in% kept_curves

      # add metadata to filtered dataset
      if (is.null(filtered_landmarks$curve_info)) {
        filtered_landmarks$curve_info <- reordered_data$curve_info[rows_to_keep, ]
      } else {
        filtered_landmarks$curve_info <- rbind(
          filtered_landmarks$curve_info,
          reordered_data$curve_info[rows_to_keep, ]
        )
      }
    }
  }

  return(filtered_landmarks)
}


# ------------------------------------------------------------------------------
# Function: test_symdist
# ------------------------------------------------------------------------------

#' Compute mean distances for symmetry reconstruction (used as an error metric)
#'
#' For each curve (excluding U/V curves) in each specimen, this function temporarily
#' removes the curve, reconstructs it using the mirrored side and ventral curve
#' with different weighting methods and lambda values, and computes mean distances
#' between the original and reconstructed curves (raw, standardized by centroid size
#' and whorl height).
#'
#' @param filtered_landmarks Nested list of landmarks as returned by extract_landmarks().
#' @param weighting_methods Character vector of weighting methods to test (e.g., c("none", "linear", "exp")).
#' @param lambda_values Numeric vector of lambda values to test (used only for "exp" weighting).
#'
#' @return A data.frame with columns:
#' \itemize{
#'   \item Specimen
#'   \item Curve
#'   \item weighting
#'   \item lambda
#'   \item Mean_Distance
#'   \item Mean_Distance_CS
#'   \item Mean_Distance_H
#' }
#' Each row corresponds to one curve reconstructed under a given weighting/lambda combination.
#'
#' @export

test_symdist <- function(filtered_landmarks,
                         weighting_methods = c("none", "linear", "exp"),
                         lambda_values = 5) {

  dist_lists <- list()

  for (w_method in weighting_methods) {
    for (lambda in lambda_values) {
      message("Processing weighting = ", w_method, " | lambda = ", lambda)

      for (specimen_id in names(filtered_landmarks$landmarks)) {
        for (curve_id in names(filtered_landmarks$landmarks[[specimen_id]])) {

          # Skip U/V curves
          if (grepl("U|V|R", curve_id)) next

          # --- Store original curve and temporarily remove it safely ---
          specimen_data <- filtered_landmarks$landmarks[[specimen_id]]
          original_pts <- specimen_data[[curve_id]]

          # --- Compute centroid size of the aperture ---
          centroid <- colMeans(original_pts)
          CS <- sqrt(sum(rowSums((original_pts - centroid)^2)))

          # --- Compute whorl height (distance between endpoints in XY plane) ---
          H <- sqrt(sum((original_pts[1, 1:2] - original_pts[nrow(original_pts), 1:2])^2))

          # Skip if curve is missing
          if (is.null(original_pts)) next

          # Remove the current curve (without deleting the specimen list)
          specimen_data[[curve_id]] <- NULL

          # Find opposite curve
          curve_num <- sub("^[A-Z]+", "", curve_id)
          curve_id_sym <- if (grepl("^R", curve_id)) paste0("L", curve_num) else paste0("R", curve_num)

          # Prepare data for symmetrization
          data_to_sym <- list(landmarks = list())
          data_to_sym$landmarks[[specimen_id]] <- list()
          data_to_sym$landmarks[[specimen_id]][[curve_id_sym]] <- specimen_data[[curve_id_sym]]
          data_to_sym$landmarks[[specimen_id]][["V"]] <- specimen_data[["V"]]
          data_to_sym$landmarks[[specimen_id]][["U"]] <- specimen_data[["U"]]

          curve_info <- filtered_landmarks$curve_info
          data_to_sym$curve_info <- filtered_landmarks$curve_info

          # --- Reconstruction ---
          reconstructed <- symmetrize_sides(data_to_sym, weighting = w_method, lambda = lambda)
          reconstructed_pts <- as.matrix(reconstructed$landmarks[[specimen_id]][[curve_id]])

          # --- Delete missing landmarks ---
          # Identify missing landmark regions from metadata
          n_ldm <- nrow(data_to_sym$landmarks[[specimen_id]][[curve_id_sym]])
          missing_ldm <- get_missing_indices(specimen_id, curve_id_sym, n_ldm, curve_info)

          # --- Compute distances ---
          # But excluded missing landmarks
          non_missing_ldm <- setdiff(1:n_ldm, missing_ldm)
          distances <- sqrt(rowSums((original_pts[non_missing_ldm,] - reconstructed_pts[non_missing_ldm,])^2))
          mean_dist <- mean(distances, na.rm = TRUE)
          mean_dist_CS <- mean(distances / CS, na.rm = TRUE)
          mean_dist_H <- mean(distances / H, na.rm = TRUE)

          # --- Store aggregated results ---
          dist_lists[[length(dist_lists) + 1]] <- data.frame(
            Specimen = specimen_id,
            Curve = curve_id,
            weighting = w_method,
            lambda = lambda,
            Mean_Distance = mean_dist,
            Mean_Distance_CS = mean_dist_CS,
            Mean_Distance_H = mean_dist_H
          )

          # Restore specimen data safely
          specimen_data[[curve_id]] <- original_pts
          filtered_landmarks$landmarks[[specimen_id]] <- specimen_data
        }
      }
    }
  }
  # Combine into one data frame
  dist_df <- do.call(rbind, dist_lists)
  return(dist_df)
}


# --------------------------------------------------
#  Function: test_reconstruct_missingldm
# --------------------------------------------------
#' Test the performance of missing-landmark reconstruction
#'
#' @description
#' This function evaluates the accuracy of the `reconstruct_missingldm()` pipeline
#' by artificially removing landmarks in selected apertures and reconstructing them
#' under different methodological conditions.
#' It allows testing the effect of:
#' * the proportion of missing landmarks,
#' * the maximum number of donor apertures,
#' * the weighting function used for donor weighting,
#' * different hierarchical structures in the donor selection.
#'
#' For each test configuration, missing landmarks are injected as `NA`, reconstructed,
#' and compared to the original coordinates. Several error metrics are computed,
#' including raw Euclidean distance, distance relative to centroid size, and distance
#' relative to whorl height.
#'
#' @param data A list containing landmark data structured as in the output of
#'   the preprocessing pipeline (`extract_landmark`, `reorder`, `symmetrise_sides`,
#'   `reconstruct_missingldm`). Must include `data$landmarks[[specimen]][[curve]]`
#'   as matrices of coordinates.
#' @param n_test Integer. Number of repeated tests to perform (used only when
#'   `param_to_test = "n_missing_ldm"`). Defaults to 10.
#' @param param_to_test Character. Specifies which parameter to explore.
#'   One of:
#'   * `"n_missing_ldm"` — tests sensitivity to varying amounts of missing data.
#'   * `"n_donors"` — tests reconstruction accuracy as a function of `max_donors`.
#'   * `"weighting_function"` — compares `"inverse"` and `"gaussian"` weighting.
#'   * `"Hierarchy"` — evaluates different hierarchical structures in donor filtering.
#'   * `"none"` — performs a baseline test with fixed parameters.
#' @param verbose Logical. If TRUE, prints progress messages. Default is TRUE.
#'
#' @details
#' For each specimen and each aperture curve not starting with `"U"` or `"V"`,
#' the function:
#' \enumerate{
#'   \item Extracts the original landmark coordinates.
#'   \item Computes two scale references:
#'     \itemize{
#'       \item Centroid size (CS)
#'       \item Whorl height (distance between first and last point in the XY-plane)
#'     }
#'   \item Injects missing landmarks (`NA`) depending on the chosen test mode.
#'   \item Reconstructs missing points using `reconstruct_missingldm()`.
#'   \item Computes Euclidean distances between reconstructed and true landmarks.
#'   \item Stores summary statistics for comparison across methods.
#' }
#'
#' The function returns a data frame summarizing reconstruction error for each test
#' configuration, enabling calibration and benchmarking of the reconstruction method.
#'
#' @return
#' A list with two elements:
#' \describe{
#'   \item{\code{results}}{A data frame containing, for each test:
#'     \itemize{
#'       \item the specimen and curve ID,
#'       \item number and percentage of reconstructed landmarks,
#'       \item mean Euclidean reconstruction error,
#'       \item error scaled by centroid size,
#'       \item error scaled by whorl height,
#'       \item values of the parameter being tested (e.g. donors, weighting type, hierarchy index).
#'     }
#'   }
#'   \item{\code{dist_by_ldm}}{(Currently empty) Reserved for future storage of
#'   landmark-wise error distributions.}
#' }
#'
#' @examples
#' \dontrun{
#' # Baseline test
#' res <- test_reconstruct_missingldm(data = data_list, param_to_test = "none")
#'
#' # Explore effect of number of donors
#' res_donors <- test_reconstruct_missingldm(data_list, param_to_test = "n_donors")
#'
#' # Compare weighting functions
#' res_weights <- test_reconstruct_missingldm(data_list, param_to_test = "weighting_function")
#' }
#'
#' @seealso
#' \code{\link{reconstruct_missingldm}} for the reconstruction method being tested.


test_reconstruct_missingldm <- function(data,
                                        n_test = 10,
                                        param_to_test = c("n_missing_ldm", "n_donors", "weighting_function", "Hierarchy", "none"),
                                        verbose = TRUE) {

  results_list <- list()
  dist_by_ldm_res <- list()
  test_index <- 1

  for (specimen_id in names(data$landmarks)) {
    curve_ids <- names(data$landmarks[[specimen_id]])

    for (curve_id in curve_ids) {
      if (grepl("^[UV]", curve_id)) next

      # Original reference
      original_pts <- data$landmarks[[specimen_id]][[curve_id]]
      n_ldm <- nrow(original_pts)

      # Compute centroid size of the aperture
      centroid <- colMeans(original_pts)
      CS <- sqrt(sum(rowSums((original_pts - centroid)^2)))

      # Compute whorl height (distance between endpoints in XY plane)
      H <- sqrt(sum((original_pts[1, 1:2] - original_pts[nrow(original_pts), 1:2])^2))


      ## ============================
      ##   1 — Vary missing ldm
      ## ============================
      if (param_to_test == "n_missing_ldm") {

        for (t in 1:n_test) {

          # two possible missing index regions
          idx_missing_1 <- sample(1:(round((n_ldm+1)/4)-3), 1)
          idx_missing_2 <- sample(round((n_ldm+2)/4):(n_ldm/2+2), 1)

          # Missing dorsal region
          idx_missing_dorsal <- c(1:idx_missing_1,
                                  (n_ldm-idx_missing_1+1):n_ldm)

          # Missing ventral region
          idx_missing_ventral <- idx_missing_2:(round(n_ldm/2)+idx_missing_2)

          # Choose 1 or 2 missing regions
          if (sample(1:2, 1) == 1) {
            idx_missing <- sample(list(idx_missing_dorsal, idx_missing_ventral), 1)[[1]]
          } else {
            idx_missing <- c(idx_missing_dorsal, idx_missing_ventral)
          }

          ## Inject NA landmarks
          temp_data <- data
          curve_data <- temp_data$landmarks[[specimen_id]][[curve_id]]
          to_na <- idx_missing[idx_missing >= 1 & idx_missing <= n_ldm]

          curve_data[to_na, ] <- NA
          temp_data$landmarks[[specimen_id]][[curve_id]] <- curve_data

          ## Reconstruction
          reconstructed <- reconstruct_missingldm(
            temp_data, metadata = spec_tab,
            max_donors = 3, hierarchy = list(),
            weighting_function = "gaussian", verbose = verbose
          )

          reconstructed_pts <- reconstructed$landmarks[[specimen_id]][[curve_id]]
          dists <- sqrt(rowSums((reconstructed_pts[to_na, ] - original_pts[to_na, ])^2))

          summary_df <- data.frame(
            Specimen = specimen_id,
            Curve = curve_id,
            Nb_reconstructed = length(to_na),
            Pourcentage = length(to_na) / n_ldm * 100,
            Mean_dist = mean(dists),
            Mean_dist_CS = mean(dists / CS, na.rm = TRUE),
            Mean_dist_H = mean(dists / H, na.rm = TRUE)
          )

          results_list[[test_index]] <- summary_df
          test_index <- test_index + 1
        }
      }

      ## ============================
      ##   2 — Vary donors
      ## ============================
      if (param_to_test == "n_donors") {

        pourcent_deleted <- 20
        nb_missing_ldm <- round(n_ldm * pourcent_deleted / 100)
        half <- floor(nb_missing_ldm / 2)
        idx_missing_dorsal <- c(1:half, (n_ldm - half + 1):n_ldm)
        idx_missing_ventral <- ((n_ldm + 1) / 2 - half):((n_ldm + 1) / 2 + half)
        idx_missing <- sample(list(idx_missing_dorsal, idx_missing_ventral), 1)[[1]]

        ## Inject NA
        temp_data <- data
        curve_data <- temp_data$landmarks[[specimen_id]][[curve_id]]
        to_na <- idx_missing[idx_missing >= 1 & idx_missing <= n_ldm]

        curve_data[to_na, ] <- NA
        temp_data$landmarks[[specimen_id]][[curve_id]] <- curve_data

        nb_curve <- sum(!grepl("^[VU]", unlist(lapply(data$landmarks, names))))
        # steps <- c(1:5, seq(10, nb_curve, by = 5))
        steps <- c(seq(1, 10, by=1))

        for (i in steps) {

          reconstructed <- reconstruct_missingldm(
            temp_data, metadata = spec_tab,
            max_donors = i, hierarchy = list(),
            weighting_function = "gaussian", verbose = FALSE
          )

          reconstructed_pts <- reconstructed$landmarks[[specimen_id]][[curve_id]]
          dists <- sqrt(rowSums((reconstructed_pts[to_na, ] - original_pts[to_na, ])^2))

          summary_df <- data.frame(
            Specimen = specimen_id,
            Curve = curve_id,
            nb_donors = i,
            Nb_reconstructed = length(to_na),
            Pourcentage = length(to_na) / n_ldm * 100,
            Mean_dist = mean(dists),
            Mean_dist_CS = mean(dists / CS, na.rm = TRUE),
            Mean_dist_H = mean(dists / H, na.rm = TRUE)
          )

          results_list[[test_index]] <- summary_df
          test_index <- test_index + 1
        }
      }

      ## ============================
      ##   3 — Vary weighting function
      ## ============================

      if (param_to_test == "weighting_function") {
        pourcent_deleted <- 20
        nb_missing_ldm <- round(n_ldm * pourcent_deleted / 100)
        half <- floor(nb_missing_ldm / 2)
        idx_missing_dorsal <- c(1:half, (n_ldm - half + 1):n_ldm)
        idx_missing_ventral <- ((n_ldm + 1) / 2 - half):((n_ldm + 1) / 2 + half)
        idx_missing <- sample(list(idx_missing_dorsal, idx_missing_ventral), 1)[[1]]

        ## Inject NA
        temp_data <- data
        curve_data <- temp_data$landmarks[[specimen_id]][[curve_id]]
        to_na <- idx_missing[idx_missing >= 1 & idx_missing <= n_ldm]

        curve_data[to_na, ] <- NA
        temp_data$landmarks[[specimen_id]][[curve_id]] <- curve_data

        weighting_function_to_test = c("inverse", "gaussian")

        for (i in weighting_function_to_test) {
          reconstructed <- reconstruct_missingldm(temp_data, metadata = spec_tab,
                                                  max_donors = 3, hierarchy = list(),
                                                  weighting_function = i, verbose = FALSE)
          reconstructed_pts <- reconstructed$landmarks[[specimen_id]][[curve_id]]
          dists <- sqrt(rowSums((reconstructed_pts[to_na, ] - original_pts[to_na, ])^2))

          summary_df <- data.frame(
            Specimen = specimen_id,
            Curve = curve_id,
            weighting_function = i,
            Nb_reconstructed = length(to_na),
            Pourcentage = length(to_na) / n_ldm * 100,
            Mean_dist = mean(dists),
            Mean_dist_CS = mean(dists / CS, na.rm = TRUE),
            Mean_dist_H = mean(dists / H, na.rm = TRUE)
          )

          results_list[[test_index]] <- summary_df
          test_index <- test_index + 1
        }
      }

      ## ============================
      ##   4 — Vary hierarchy structure
      ## ============================

      if (param_to_test == "Hierarchy") {

        pourcent_deleted <- 20
        nb_missing_ldm <- round(n_ldm * pourcent_deleted / 100)
        half <- floor(nb_missing_ldm / 2)
        idx_missing_dorsal <- c(1:half, (n_ldm - half + 1):n_ldm)
        idx_missing_ventral <- ((n_ldm + 1) / 2 - half):((n_ldm + 1) / 2 + half)
        idx_missing <- sample(list(idx_missing_dorsal, idx_missing_ventral), 1)[[1]]

        ## Inject NA
        temp_data <- data
        curve_data <- temp_data$landmarks[[specimen_id]][[curve_id]]
        to_na <- idx_missing[idx_missing >= 1 & idx_missing <= n_ldm]

        curve_data[to_na, ] <- NA
        temp_data$landmarks[[specimen_id]][[curve_id]] <- curve_data

        all_hierarchy_lists <- list(
          list(),
          list("Specimen"),
          list("Specimen", "Species"),
          list("Specimen", "Species", "Genus")
        )

        for (i in seq_along(all_hierarchy_lists)) {
          reconstructed <- reconstruct_missingldm(temp_data, metadata = spec_tab,
                                                  max_donors = 3, hierarchy = all_hierarchy_lists[[i]],
                                                  weighting_function = "gaussian", verbose = FALSE)

          reconstructed_pts <- reconstructed$landmarks[[specimen_id]][[curve_id]]
          dists <- sqrt(rowSums((reconstructed_pts[to_na, ] - original_pts[to_na, ])^2))

          summary_df <- data.frame(
            Specimen = specimen_id,
            Curve = curve_id,
            hierarchy = i,
            Nb_reconstructed = length(to_na),
            Pourcentage = length(to_na) / n_ldm * 100,
            Mean_dist = mean(dists),
            Mean_dist_CS = mean(dists / CS, na.rm = TRUE),
            Mean_dist_H = mean(dists / H, na.rm = TRUE)
          )

          results_list[[test_index]] <- summary_df
          test_index <- test_index + 1
        }
      }

      ## ============================
      ##   5 — Vary nothing
      ## ============================

      if (param_to_test == "none") {
        pourcent_deleted <- 20
        nb_missing_ldm <- round(n_ldm * pourcent_deleted / 100)
        half <- floor(nb_missing_ldm / 2)
        idx_missing_dorsal <- c(1:half, (n_ldm - half + 1):n_ldm)
        idx_missing_ventral <- ((n_ldm + 1) / 2 - half):((n_ldm + 1) / 2 + half)
        idx_missing <- sample(list(idx_missing_dorsal, idx_missing_ventral), 1)[[1]]

        ## Inject NA
        temp_data <- data
        curve_data <- temp_data$landmarks[[specimen_id]][[curve_id]]
        to_na <- idx_missing[idx_missing >= 1 & idx_missing <= n_ldm]

        curve_data[to_na, ] <- NA
        temp_data$landmarks[[specimen_id]][[curve_id]] <- curve_data

        reconstructed <- reconstruct_missingldm(temp_data, metadata = spec_tab,
                                                max_donors = 5, hierarchy = list(),
                                                weighting_function = "gaussian", verbose = FALSE)

        reconstructed_pts <- reconstructed$landmarks[[specimen_id]][[curve_id]]
        dists <- sqrt(rowSums((reconstructed_pts[to_na, ] - original_pts[to_na, ])^2))

        summary_df <- data.frame(
          Specimen = specimen_id,
          Curve = curve_id,
          Nb_reconstructed = length(to_na),
          Pourcentage = length(to_na) / n_ldm * 100,
          Mean_dist = mean(dists),
          Mean_dist_CS = mean(dists / CS, na.rm = TRUE),
          Mean_dist_H = mean(dists / H, na.rm = TRUE)
        )

        results_list[[test_index]] <- summary_df
        test_index <- test_index + 1

      }
    }
  }
  return(list(
    results = do.call(rbind, results_list),
    dist_by_ldm = do.call(rbind, dist_by_ldm_res)
  ))
}



