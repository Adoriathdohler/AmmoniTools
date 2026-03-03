
# --------------------------------------------------
# MISSING DATA RECONSTRUCTION TEST
#
# filter_curves : Filters aperture curves based on landmark completeness to define reference and target datasets for reconstruction performance assessment
# test_symmetrize_sides : Quantifies reconstruction accuracy obtained by symmetry-based completion of aperture landmarks
# test_reconstruct_missingldm : Quantifies reconstruction accuracy obtained by interperistome TPS-based completion of aperture landmarks
# --------------------------------------------------

utils::globalVariables(c("curve_num ", "spec_tab"))

# ---------------------------------------------------------------------
# Function: filter_curves
# ---------------------------------------------------------------------

#' Filters aperture curves based on landmark completeness to define reference and target datasets for reconstruction performance assessment
#'
#' This function filters a reordered landmark dataset to retain only valid
#' left/right (L/R) aperture curve pairs for each specimen.
#' Umbilical (U) and ventral (V) curves are always retained for contextual
#' and geometric reference, but are not used for left/right pairing.
#'
#' Optionally, the function can restrict the dataset to curves without
#' missing landmarks.
#'
#' @param reordered_data A list, typically from \code{reordered_landmarks}, containing:
#' \itemize{
#'   \item \code{landmarks}: a nested list of 3D landmark coordinates
#'   for each specimen and each curve.
#'   \item \code{curve_info}: a data frame containing metadata for each curve.
#'   Required columns include \code{Specimen}, \code{Curve}, \code{Type}
#'   (L, R, U, V), and \code{missing_ldm}.
#' }
#'
#' @param complete_only Logical, default = \code{FALSE}.
#' If \code{TRUE}, only L/R curves without missing landmarks
#' (i.e., empty \code{missing_ldm}) are retained.
#'
#' @return
#' A list with the same structure as \code{reordered_data}, containing:
#' \describe{
#'   \item{\code{landmarks}}{Filtered landmark coordinates for specimens
#'   that possess at least one valid L/R curve pair.}
#'   \item{\code{curve_info}}{Corresponding subset of curve metadata
#'   for the retained curves.}
#' }
#'
#' @details
#' For each specimen, the function proceeds as follows:
#' \enumerate{
#'   \item Ventral (V) and umbilical (U) curves are always retained.
#'   \item Left (L) and right (R) curves are optionally filtered based on
#'   completeness.
#'   \item Only L/R curves sharing the same numerical suffix
#'   (e.g., L1 and R1) are considered valid pairs.
#'   \item Specimens without any valid L/R pair are excluded from the output.
#' }
#'
#' This filtering step is typically used prior to symmetry-based
#' reconstruction or validation procedures, where paired lateral curves
#' are required.
#'
#' @examples
#' \dontrun{
#' # Load example dataset
#' data("example_data", package = "AmmoniTools")
#'
#' # Resample and reorder landmarks
#' resampled_data <- resample_curves(
#'   example_data,
#'   n_landmarks = list(R = 25, L = 25, U = 25, V = 25)
#' )
#' reordered_data <- reorder_landmarks(resampled_data)
#'
#' # Keep all paired L/R curves (including incomplete ones)
#' filtered_all <- filter_curves(reordered_data, complete_only = FALSE)
#'
#' # Keep only fully complete L/R curve pairs
#' filtered_complete <- filter_curves(reordered_data, complete_only = TRUE)
#'
#' # Inspect retained curves for one specimen
#' names(filtered_complete$landmarks[[1]])
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
# Function: test_symmetrize_sides
# ------------------------------------------------------------------------------

#' Quantifies reconstruction accuracy obtained by symmetry-based completion of aperture landmarks
#'
#' This function assesses the accuracy of landmark reconstruction based on
#' bilateral symmetry (\link{symmetrize_sides})by performing a leave-one-curve-out procedure.
#' For each specimen and for each lateral curve (excluding U, V and R curves),
#' the curve is temporarily removed, reconstructed using its mirrored counterpart
#' together with the ventral (V) and umbilical (U) curves, and then compared to
#' the original curve.
#'
#' Reconstruction is performed using different weighting strategies, and the
#' discrepancy between original and reconstructed landmarks is quantified using
#' Euclidean distances, expressed as:
#' \itemize{
#'   \item raw distances,
#'   \item distances standardized by centroid size (CS),
#'   \item distances standardized by whorl height (H).
#' }
#'
#' Missing landmark regions defined in the metadata are excluded from the
#' distance calculations.
#'
#' @param filtered_landmarks A list, typically from \link{filter_curves} outputs, containing:
#' \itemize{
#'   \item \code{landmarks}: nested list of landmark coordinates per specimen and curve,
#'   \item \code{curve_info}: data frame describing curves and missing landmarks.
#' }
#' Typically obtained after preprocessing with
#' \code{extract_landmarks()} and related functions.
#'
#' @param weighting_methods Character vector specifying the weighting schemes
#' used for symmetry-based reconstruction.
#' Possible values include \code{"none"}, \code{"linear"}, and \code{"exp"}.
#'
#' @param weight_power_values Numeric vector of weight power values.
#' Only used when \code{weighting_methods = "exp"}.
#'
#' @return A data.frame where each row corresponds to one reconstructed curve
#' under a given weighting configuration, with the following columns:
#' \itemize{
#'   \item \code{Specimen}: specimen identifier
#'   \item \code{Curve}: reconstructed curve name
#'   \item \code{weighting}: weighting method used
#'   \item \code{weight_power}: weight power value
#'   \item \code{Mean_Distance}: mean Euclidean distance between original and reconstructed landmarks
#'   \item \code{Mean_Distance_CS}: mean distance standardized by centroid size
#'   \item \code{Mean_Distance_H}: mean distance standardized by whorl height
#' }
#'
#' @details
#' This function is intended as a validation tool for symmetry-based
#' reconstruction pipelines. By reconstructing curves that are originally
#' complete, it provides an empirical estimate of reconstruction error
#' under different weighting assumptions.
#'
#' @seealso \link{symmetrize_sides}
#'
#' @examples
#' \dontrun{
#' # Load example dataset
#' data("example_data", package = "AmmoniTools")
#' data("example_spec_tab", package = "AmmoniTools")
#'
#' # Preprocess landmarks
#' reordered_data <- resample_curves(example_data,
#'                              n_landmarks = list(R = 25, L = 25, U = 25, V = 25))
#' reordered_data <- reorder_landmarks(reordered_data)
#'
#' filtered_data <- filter_curves(reordered_data, complete_only=FALSE)
#'
#' # Test symmetry-based reconstruction accuracy
#' dist_df <- test_symmetrize_sides(
#'   filtered_landmarks = filtered_data,
#'   weighting_methods = c("none", "linear", "exp"),
#'   weight_power_values = c(2, 5, 10)
#' )
#'
#' # Visualize reconstruction error
#' boxplot(Mean_Distance_CS ~ weighting, data = dist_df,
#'         ylab = "Mean distance (CS-standardized)",
#'         xlab = "Weighting method")
#' }
#'
#' @export

test_symmetrize_sides <- function(filtered_landmarks,
                         weighting_methods = c("none", "linear", "exp"),
                         weight_power_values = 5) {

  dist_lists <- list()

  for (w_method in weighting_methods) {
    for (weight_power in weight_power_values) {
      message("Processing weighting = ", w_method, " | weight_power = ", weight_power)

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
          reconstructed <- symmetrize_sides(data_to_sym, weighting = w_method, weight_power = weight_power)
          reconstructed_pts <- as.matrix(reconstructed$landmarks[[specimen_id]][[curve_id]])

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

          # --- Delete missing landmarks ---
          # Identify missing landmark regions from metadata
          n_ldm <- nrow(data_to_sym$landmarks[[specimen_id]][[curve_id_sym]])
          missing_ldm <- get_missing_indices(specimen_id, curve_id_sym, n_ldm)

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
            weight_power = weight_power,
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
#' Quantifies reconstruction accuracy obtained by interperistome TPS-based completion of aperture landmarks
#'
#' @description
#' This function evaluates the accuracy of the \code{\link{reconstruct_missingldm}} pipeline
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
#'   \item Reconstructs missing points using \code{\link{reconstruct_missingldm}}.
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
#' @seealso \code{\link{reconstruct_missingldm}}
#' @examples
#' \dontrun{
#'
#' # Load example dataset
#' data("example_data", package = "AmmoniTools")
#' data("example_spec_tab", package = "AmmoniTools")
#'
#' # Resample curves
#' resampled_data <- resample_curves(example_data, n_landmarks = list(R = 25, L = 25, U = 25, V = 25))
#'
#' # Reorder ldm
#' reordered_data <- reorder_landmarks(resampled_data)
#'
#' # Symmetrize
#' data_sym <- symmetrize_sides(reordered_data, weighting = "exp")
#'
#' #  Filter peristome without missing landmarks
#' filtered_landmarks_nonmerged <- filter_curves(data_sym, complete_only=TRUE)
#'
#' # Merge R and L sides
#' data_complete <- merge_sides(filtered_landmarks_nonmerged)
#'
#' # Testing the effet of proportion of missing landmarks
#' results_n_missing_ldm <- test_reconstruct_missingldm(data_complete, n_test = 30,
#' param_to_test = "n_missing_ldm")
#'
#' # Testing the effet of maximum number of donor apertures
#' results_n_donors <- test_reconstruct_missingldm(data_complete,
#'                                                 param_to_test = "n_donors")
#'
#' # Testing the effect of weighting function used for donor weighting
#' results_weighting_function <- test_reconstruct_missingldm(data_complete,
#'                                   param_to_test = "weighting_function")
#'
#' # Testing the effect of different hierarchical structures in the donor selection
#' results_Hierarchy <- test_reconstruct_missingldm(data_complete, param_to_test = "Hierarchy")
#' }
#' @seealso \code{\link{reconstruct_missingldm}} for the reconstruction method being tested.
#' @export

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



