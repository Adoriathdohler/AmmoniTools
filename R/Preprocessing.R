
# --------------------------------------------------
#  PREPROCESSING
#
# extract_landmark : Extracts et formats 3D landmarks from 3D Slicer Markups (.json) files
# resample_curves : Resamples 3D aperture curves to a fixed number of points
# view_specimens : Visualizes 3D landmarks of ammonite specimens
# reorder_landmarks : Reorders 3D aperture landmark based on ammonite anatomy
# subset_data : Subsets 3D aperture data and metadata
# --------------------------------------------------

utils::globalVariables("spec_tab")

# --------------------------------------------------
#  Function: extract_landmarks
# --------------------------------------------------

#' Extracts et formats 3D landmarks from 3D Slicer Markups (.json) files
#'
#' This function reads 3D landmark coordinates from Slicer Markups JSON files
#' organized in subfolders, each representing a specimen. It identifies missing
#' landmarks and checks for curve completeness.
#' Each subfolder must contain one or more `.json` files describing landmark curves.
#' The expected number of landmarks for each curve type (R, L, U, V, etc.) can be
#' specified via the `expected_ldm` argument.
#'
#' @param folder Path to the parent folder containing one subfolder per specimen.
#'   Each subfolder should contain the `.json` markup files (3D Slicer Markup ouputs) for that specimen.
#' @param spec_tab A data.frame containing specimen metadata (must include a `Specimen` column).
#'   Only specimens listed in `spec_tab` and passing optional filters are processed.
#' @param expected_ldm A named list giving the expected number of landmarks for each
#'   curve type (e.g., `list(R = 200, L = 200, U = 100, V = 100)`).
#' @param verbose Logical; if `TRUE`, progress and specimen IDs are printed during processing.
#' @param ... Optional filtering arguments matching columns in `spec_tab`.
#'   For example, `Species = "opalinum", Genus = "Harpoceras"`.
#'
#' @return A list containing:
#' \describe{
#'   \item{curve_info}{A `data.frame` summarizing all curves and detected missing landmarks.}
#'   \item{landmarks}{A nested list of 3D landmark coordinates by specimen and curve.}
#'   \item{messages}{A list of text warnings or notes for each specimen.}
#' }
#'
#' @details
#' For each specimen:
#' \enumerate{
#'   \item All `.json` markup files in its folder are read.
#'   \item Landmark coordinates are extracted as 3D matrices (X, Y, Z).
#'   \item Curves missing landmarks labeled `x` or `X` are flagged.
#'   \item A warning is issued if required curves (e.g., "U", "V") are missing.
#' }
#'
#' @examples
#' \dontrun{
#' # Extract all specimen from my_spec_tab
#' raw_data <- extract_landmarks(
#'     "path/to/data",
#'     spec_tab = my_spec_tab,
#'     Species = "opalinum"   # filter applied on 'Species' colonn of spec_tab
#'     )
#'
#' # Extract one particular species
#' raw_data <- extract_landmarks(
#'     "path/to/data",
#'     spec_tab = my_spec_tab,
#'     Species = "opalinum"   # filter applied on 'Species' colonn of spec_tab
#'     )
#'
#' # Extract two particular genus and one particular species
#' raw_data <- extract_landmarks(
#'     "path/to/data",
#'     spec_tab = my_spec_tab,
#'     Species = "opalinum", Genus = c("Harpoceras", "Polypectus")  # combined filters
#'     )
#'}
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr bind_rows
#' @importFrom purrr map map_chr
#' @export

extract_landmarks <- function(folder, spec_tab,
                              expected_ldm = list(R = 200, L = 200, U = 100, V = 100),
                              verbose = TRUE, ...) {
  # --------------------------------------------------
  # Step 1: Apply optional filters to select specimens
  # --------------------------------------------------
  # Filters can be passed like Species = "opalinum", Genus = "Harpoceras"
  filters <- list(...)
  filtered_spec_tab <- spec_tab

  # Apply all provided filters to the specimen metadata
  for (filter_name in names(filters)) {
    filtered_spec_tab <- filtered_spec_tab[
      filtered_spec_tab[[filter_name]] %in% filters[[filter_name]],
    ]
  }

  # Extract the specimen IDs that will actually be processed
  kept_ids <- as.character(filtered_spec_tab$Specimen)

  # --------------------------------------------------
  # Step 2: Identify all specimen folders in the parent directory
  # --------------------------------------------------
  # Each subfolder should correspond to one specimen ID
  sub_folders <- list.dirs(folder, full.names = TRUE, recursive = FALSE)
  sub_folders <- sub_folders[basename(sub_folders) %in% kept_ids]

  # --------------------------------------------------
  # Step 3: Define the function that processes a single specimen
  # --------------------------------------------------
  process_specimen <- function(sub_folder) {
    specimen_id <- basename(sub_folder)
    messages_local <- character()  # Collect local warnings and messages

    # Utility function to append formatted messages
    add_msg <- function(...) {
      msg <- paste0(...)
      messages_local <<- c(messages_local, msg)
    }

    # List all JSON files in the specimen folder
    json_files <- list.files(sub_folder, pattern = "\\.json$", full.names = TRUE)
    if (length(json_files) == 0) {
      add_msg("\u26A0\uFE0F No JSON files found for specimen ", specimen_id)
      return(list(
        specimen_id = specimen_id,
        landmarks = NULL,
        curve_info = NULL,
        messages = messages_local
      ))
    }

    curve_landmarks <- list()         # Store coordinates for each curve
    curve_info_list <- vector("list", length(json_files))  # Metadata per curve

    # Iterate through all .json files for this specimen
    for (j in seq_along(json_files)) {
      json_file <- json_files[j]

      # Safe read: catch invalid JSONs
      json_data <- tryCatch(
        jsonlite::fromJSON(json_file),
        error = function(e) {
          add_msg("\u26A0\uFE0F Could not read JSON file: ", json_file)
          return(NULL)
        }
      )
      if (is.null(json_data)) next

      # Check if the JSON actually contains landmark coordinates
      if (is.null(json_data$markups$controlPoints[[1]]$position)) {
        add_msg("\u26A0\uFE0F Skipped curve ", json_file, " (invalid JSON structure - no landmarks)")
        next
      }

      # Extract coordinates and label information
      coordinates <- as.data.frame(do.call(rbind, json_data$markups$controlPoints[[1]]$position))
      colnames(coordinates) <- c("X", "Y", "Z")

      # Extract curve name and type (R, L, U, V, etc.)
      curve_name <- sub("\\.mrk$", "", tools::file_path_sans_ext(basename(json_file)))
      curve_type <- substr(curve_name, 1, 1)

      # Compare number of landmarks to expected value
      expected <- expected_ldm[[curve_type]]
      n_ldm <- nrow(coordinates)
      if (!is.null(expected) && n_ldm != expected) {
        add_msg("\u26A0\uFE0F Specimen ", specimen_id, " - ", curve_name,
                ": got ", n_ldm, " landmarks (expected ", expected, ")")
      }

      # Identify missing landmarks labeled 'x' or 'X', excluding endpoints
      missing_ldm <- which(grepl("x", json_data$markups$controlPoints[[1]]$label, ignore.case = TRUE))
      missing_ldm <- missing_ldm[missing_ldm != expected & missing_ldm != 1]

      # Store the coordinates and curve metadata
      curve_landmarks[[curve_name]] <- coordinates
      curve_info_list[[j]] <- data.frame(
        Specimen = specimen_id,
        Curve = curve_name,
        Type = curve_type,
        missing_ldm = I(list(missing_ldm)),
        stringsAsFactors = FALSE
      )
    }

    # --------------------------------------------------
    #  Check for essential curves
    # --------------------------------------------------
    curve_names <- names(curve_landmarks)
    has_U <- any(curve_names %in% c("U", "Ubis"))
    has_V <- "V" %in% curve_names
    if (!has_U || !has_V) {
      add_msg("\u26A0\uFE0F Missing required curves for specimen ", specimen_id, ": ",
              if (!has_U) "U or Ubis " else "",
              if (!has_U & !has_V) "and " else "",
              if (!has_V) "V" else "")
    }

    # Return structured output for this specimen
    list(
      specimen_id = specimen_id,
      landmarks = curve_landmarks,
      curve_info = dplyr::bind_rows(curve_info_list),
      messages = messages_local
    )
  }

  # --------------------------------------------------
  # Step 4: Run sequentially (with progress display)
  # --------------------------------------------------
  total <- length(sub_folders)
  if (verbose) message("\U1F9E9 Running in sequential mode (", total, " specimens)")

  # Process each specimen sequentially and print progress
  results <- purrr::map2(sub_folders, seq_along(sub_folders), function(x, i) {
    specimen_id <- basename(x)
    if (verbose) cat("\u27A1\uFE0F [", i, "/", total, "] Processing specimen", specimen_id, "\n")
    process_specimen(x)
  })

  # --------------------------------------------------
  # Step 5: Aggregate all results
  # --------------------------------------------------
  landmarks <- purrr::map(results, "landmarks")
  names(landmarks) <- purrr::map_chr(results, "specimen_id")

  curve_info <- dplyr::bind_rows(purrr::map(results, "curve_info"))
  messages_all <- purrr::map(results, "messages") |> stats::setNames(names(landmarks))

  # --------------------------------------------------
  # Step 6: Display all messages at the end
  # --------------------------------------------------
  cat("\n \U1F4CB --- Summary of warnings/messages ---\n")
  for (id in names(messages_all)) {
    msgs <- messages_all[[id]]
    if (length(msgs) > 0) {
      cat("\n\U1F539", id, ":\n")
      cat(paste0("   ", msgs, collapse = "\n"), "\n")
    }
  }
  cat("\n\u2705 Done. Processed", length(results), "specimens.\n")

  # --------------------------------------------------
  # Step 7: Return structured output
  # --------------------------------------------------
  list(
    curve_info = curve_info,
    landmarks = landmarks,
    messages = messages_all
  )
}


# --------------------------------------------------
#  Function : view_specimens
# --------------------------------------------------

#' Visualizes 3D landmarks of ammonite specimens
#'
#' This function allows you to visualize one or several ammonite specimens
#' in 3D using **rgl**. Curves are colored according to their type
#' (R = red, L = blue, U/V = green). You can either display a single specimen
#' (by ID), all specimens one after the other, or interactively choose from a menu.
#'
#' @param landmarks_list A named list of specimens. Each specimen must be a list
#'   of curves, where each curve is a matrix of 3D coordinates (X, Y, Z).
#' @param specimen_id A character string specifying the ID of a specimen to display.
#'   Must match one of the names of \code{landmarks_list}. If \code{NULL} (default),
#'   an interactive menu will be displayed unless \code{all = TRUE}.
#' @param all Logical, if \code{TRUE} all specimens in \code{landmarks_list} are displayed
#'   sequentially. The user must press \code{Enter} to continue to the next specimen and \code{Escape}
#'   to escape the plotting.
#' @param show_points Logical. If \code{TRUE}, landmark points are displayed
#'   using \code{sphere3d()}.
#' @param point_size Numeric. Size of landmark points (default = 0.2).
#'
#' @details
#' Colors are assigned automatically based on the first letter of the curve name:
#' \itemize{
#'   \item R = red
#'   \item L = blue
#'   \item U and V = green
#'   \item Any other curve (e.g. peristome curve after merging left and right sides
#'   ) = black
#' }
#'
#' @return This function is called for its side effect: opening an interactive
#'   3D window via \pkg{rgl}. It does not return a value.
#'
#' @examples
#' \dontrun{
#' # Example using AmmoniTools example data
#' data("example_data", package = "AmmoniTools")
#'
#' # Display one specimen
#' view_specimens(example_data$landmarks, specimen_id = "182")
#'
#' # Browse all specimens
#' view_specimens(example_data$landmarks, all = TRUE)
#'
#' # Choose interactively
#' view_specimens(example_data$landmarks)
#' }
#'
#' @import rgl
#' @export
view_specimens <- function(landmarks_list,
                           specimen_id = NULL,
                           all = FALSE,
                           show_points = FALSE,
                           point_size = 0.2) {

  curve_colors <- c("R" = "red",
                    "L" = "blue2",
                    "U" = "chartreuse2",
                    "V" = "chartreuse2")

  plot_one <- function(specimen, id) {

    par3d(windowRect = c(50, 50, 800, 800))
    clear3d()
    title3d(main = id)

    for (curve_name in names(specimen)) {

      curve <- specimen[[curve_name]]
      if (is.null(curve)) next

      curve_type <- substr(curve_name, 1, 1)

      color <- ifelse(curve_type %in% names(curve_colors),
                      curve_colors[[curve_type]],
                      "black")

      # Draw curve
      lines3d(curve, col = color, lwd = 2)

      # Optionally draw landmarks
      if (show_points) {
        spheres3d(curve,
                 col = color,
                 radius = point_size)
      }
    }
  }

  if (all) {
    for (id in names(landmarks_list)) {
      plot_one(landmarks_list[[id]], id)
      readline(prompt = "Press [Enter] to continue...")
      rgl::close3d()
    }
  } else if (!is.null(specimen_id)) {

    if (!(specimen_id %in% names(landmarks_list))) {
      stop("Specimen ID not found in landmarks_list.")
    }

    plot_one(landmarks_list[[specimen_id]], specimen_id)

  } else {

    choice <- utils::menu(names(landmarks_list),
                          title = "Select a specimen to display:")

    if (choice > 0) {
      id <- names(landmarks_list)[choice]
      plot_one(landmarks_list[[id]], id)
    }
  }

  invisible(NULL)
}
# --------------------------------------------------
#  Function: resample_curves
# --------------------------------------------------

#' Resamples 3D aperture curves to a fixed number of points
#'
#' This function uniformly resamples the landmarks along each curve of each specimen
#' in a dataset such as \code{raw_data}. It allows to reduce or increase the number
#' of landmarks per curve while preserving their spatial distribution.
#'
#' @param raw_data A list containing at least:
#'   \itemize{
#'     \item \code{landmarks}: a named list of specimens, where each specimen is itself
#'     a list of curves, and each curve is a numeric matrix of 3D coordinates (X, Y, Z).
#'     \item \code{curve_info}: a data frame describing the curves, with at least the columns
#'     \code{Specimen}, \code{Curve}, \code{Type}, and \code{missing_ldm}.
#'   }
#' @param n_landmarks A named list defining the desired number of landmarks per curve type
#'   (e.g., \code{list(R = 50, L = 50, U = 50, V = 50)}). Curve types are identified by
#'   the first letter of each curve name. If a type is not provided, its number of landmarks
#'   remains unchanged.
#'
#' @details
#' The function resamples each curve using linear interpolation between existing landmarks.
#' The resampling preserves the geometric shape of the curve while adjusting the number of
#' landmarks to match the specified target values. The interpolation is performed independently
#' for each curve in each specimen.
#'
#' If the input \code{raw_data} contains a \code{curve_info} table, this function will:
#' \itemize{
#'   \item Update the column \code{n_ldm} to reflect the new number of landmarks.
#'   \item Rescale the indices of missing landmarks (\code{missing_ldm}) proportionally to the
#'   new sampling density.
#'   \item Preserve the structure of \code{curve_info}, ensuring \code{missing_ldm} remains
#'   a list of integer vectors.
#' }
#' It can be used, for instance, to standardize the sampling resolution between
#' left/right apertures or to reduce file size for rapid analyses.
#'
#' @return
#' A list of the same structure as \code{raw_data}, with:
#' \itemize{
#'   \item \code{landmarks}: resampled 3D coordinates for each curve.
#'   \item \code{curve_info}: an updated data frame with adjusted \code{n_ldm} and
#'   rescaled \code{missing_ldm}.
#' }
#'
#' @examples
#' \dontrun{
#' # Example using AmmoniTools example data
#' data("example_data", package = "AmmoniTools")
#'
#' # Resample all curves to 50 landmarks each
#' resampled_data <- resample_curves(example_data, n_landmarks = list(R = 50, L = 50, U = 50, V = 50))
#'
#' # Visualize one specimen before and after
#' view_specimens(example_data$landmarks, specimen_id = "182")
#' view_specimens(resampled$landmarks, specimen_id = "182")
#' }
#'
#' @importFrom stats approx
#' @export

resample_curves <- function(raw_data,
                            n_landmarks = list(R = 50, L = 50, U = 100, V = 100)) {

  # --- Helper function: resample one curve ---
  resample_curve <- function(curve, n_points) {
    curve <- as.data.frame(curve)
    curve <- curve[, sapply(curve, is.numeric), drop = FALSE]
    if (ncol(curve) < 3 || nrow(curve) < 2) return(curve)
    curve <- stats::na.omit(curve)

    dists <- sqrt(rowSums(diff(as.matrix(curve))^2))
    cum_dists <- c(0, cumsum(dists))
    total_length <- max(cum_dists)
    new_dists <- seq(0, total_length, length.out = n_points)

    new_curve <- cbind(
      approx(cum_dists, curve[, 1], xout = new_dists)$y,
      approx(cum_dists, curve[, 2], xout = new_dists)$y,
      approx(cum_dists, curve[, 3], xout = new_dists)$y
    )
    colnames(new_curve) <- colnames(curve)
    return(new_curve)
  }

  # --- Initialize output ---
  new_data <- raw_data
  new_data$landmarks <- list()

  # --- Loop over all specimens ---
  for (specimen_name in names(raw_data$landmarks)) {
    specimen <- raw_data$landmarks[[specimen_name]]
    new_specimen <- list()

    for (curve_name in names(specimen)) {
      curve <- specimen[[curve_name]]
      if (is.null(curve) || nrow(as.data.frame(curve)) < 2) {
        new_specimen[[curve_name]] <- curve
        next
      }

      curve_type <- substr(curve_name, 1, 1)
      n_target <- ifelse(curve_type %in% names(n_landmarks),
                         n_landmarks[[curve_type]],
                         nrow(curve))

      new_specimen[[curve_name]] <- resample_curve(curve, n_target)
    }

    new_data$landmarks[[specimen_name]] <- new_specimen
  }

  # --- Update curve_info with new landmark counts and rescaled missing indices ---
  if (!is.null(raw_data$curve_info)) {
    new_info <- raw_data$curve_info

    # Ensure we have a column for the number of landmarks
    if (!"n_ldm" %in% names(new_info)) {
      new_info$n_ldm <- NA_integer_
    }

    for (i in seq_len(nrow(new_info))) {
      specimen_id <- as.character(new_info$Specimen[i])
      curve_name  <- as.character(new_info$Curve[i])

      # Determine curve type (first letter)
      curve_type <- substr(curve_name, 1, 1)

      # Old number of landmarks
      old_n <- nrow(raw_data$landmarks[[specimen_id]][[curve_name]])

      # Determine new number of landmarks for this curve type
      if (curve_type %in% names(n_landmarks)) {
        new_n <- n_landmarks[[curve_type]]
      } else {
        new_n <- old_n  # fallback if not specified
      }

      # Update n_ldm
      new_info$n_ldm[i] <- new_n

      # --- Handle missing landmarks (list of integer vectors) ---
      old_missing <- new_info$missing_ldm[[i]]

      if (is.numeric(old_missing) && length(old_missing) > 0) {
        # Rescale missing indices proportionally to new sampling
        scaled_missing <- round(old_missing / old_n * new_n)
        scaled_missing <- unique(pmin(pmax(scaled_missing, 1), new_n)) # stay within bounds
        new_info$missing_ldm[[i]] <- scaled_missing
      } else {
        # If no missing landmarks, ensure it's an empty integer vector
        new_info$missing_ldm[[i]] <- integer(0)
      }
    }

    # Assign updated info to output
    new_data$curve_info <- new_info
  }

  return(new_data)
}

# --------------------------------------------------
#  Function: reorder_landmarks
# --------------------------------------------------

#' Reorders 3D aperture landmark based on ammonite anatomy
#'
#' This function reorders (i.e., reverses) the curves of each specimen when necessary,
#' based on their spatial relationship with the ventral curve (`V`). This ensures
#' all curves follow a consistent anatomical direction.
#'
#' The function also adjusts the indices of missing landmarks in the `curve_info`
#' table accordingly when a curve is reversed.
#'
#' @param data A list containing:
#' \describe{
#'   \item{landmarks}{A named list of specimens. Each specimen is a list of curves
#'   (matrices of 2D or 3D coordinates), including at least the `V` curve and optionally
#'   `L1`, `R1`, `U`, `Ubis`, etc.}
#'   \item{curve_info}{A data frame describing curve metadata. Must include the following columns:
#'     \itemize{
#'       \item \code{Specimen}: specimen ID
#'       \item \code{Curve}: curve name
#'       \item \code{missing_ldm}: list column with indices of missing landmarks for each curve
#'     }
#'   }
#' }
#'
#'
#' @return A list with:
#' \describe{
#'   \item{`landmarks`}{The reordered list of landmark curves, with corrected orientation.}
#'   \item{`curve_info`}{The updated `curve_info` data frame, with corrected indices
#'   for reversed curves.}
#' }
#'
#' @details
#' The logic follows anatomical expectations:
#' \itemize{
#'   \item For `R` curves: the last point should be closer to the `V` curve than the first.
#'   \item For `L` curves: the first point should be closer to the `V` curve than the last.
#'   \item The ventral curve (`V`) is oriented based on the attachment point of `L1` or `R1`.
#'   \item Umbilical curves (`U`, `Ubis`, etc.) are also oriented using the same reference point.
#' }
#'
#' A progress bar is displayed during processing. At the end, a message lists all reordered curves.
#'
#' @examples
#' \dontrun{
#' # Example using AmmoniTools example data
#' data("example_data", package = "AmmoniTools")
#'
#' reordered_data <- reorder_landmarks(example_data)
#' }
#'
#' @export

reorder_landmarks <- function(data) {

  # Extract curve info and landmarks from input data
  curve_info_ord <- data$curve_info
  landmarks_ord <- data$landmarks
  reordered_log <- list()  # inversed curve storage

  # progress bar setup
  specimen_ids <- names(landmarks_ord)
  pb <- progress::progress_bar$new(
    format = "  Reordering [:bar] :percent (:current/:total) | Specimen: :specimen",
    total = length(specimen_ids),
    clear = FALSE,
    width = 80
  )

  # --------------------------------------------------
  # Helper function: find closest ventral point (from V curve)
  # to a given 3D point
  # --------------------------------------------------
  find_closest_ventral_point_to_point <- function(point, V_curve) {
    distances <- apply(V_curve, 1, function(v_point) {
      sqrt(sum((point - v_point)^2))
    })
    return(which.min(distances))
  }

  # --------------------------------------------------
  # Helper function: reverse (flip) a curve
  # - Updates curve order (reverse row order)
  # - Logs the reordering
  # - Adjusts indices of missing landmarks accordingly
  # --------------------------------------------------
  reverse_curve <- function(curve, specimen_id, curve_name, V_curve, curve_info_ord) {
    curve <- curve[nrow(curve):1, ] # reverse curve
    reordered_log[[length(reordered_log) + 1]] <<- list(specimen = specimen_id, curve = curve_name)

    # Update missing landmark indices if present
    curve_info_idx <- which(curve_info_ord$Specimen == specimen_id & curve_info_ord$Curve == curve_name)
    if (length(curve_info_idx) > 0 && length(curve_info_ord$missing_ldm[[curve_info_idx]]) > 0) {
      max_ldm <- max(curve_info_ord$missing_ldm[[curve_info_idx]])
      curve_info_ord$missing_ldm[[curve_info_idx]] <- max_ldm + 1 - curve_info_ord$missing_ldm[[curve_info_idx]]
    }

    return(list(curve = curve, curve_info_ord = curve_info_ord))
  }

  updated_curve_info <- curve_info_ord


  # --------------------------------------------------
  # Main loop over specimens
  # --------------------------------------------------
  landmarks_ord <- lapply(names(landmarks_ord), function(specimen_id) {

    pb$tick(tokens = list(specimen = specimen_id)) # update progression bar

    specimen_data <- landmarks_ord[[specimen_id]]

    # Extract curves: L/R, V, U
    LR_curves <- specimen_data[grep("^[LR][0-9]+$", names(specimen_data))]
    V_curve <- specimen_data[["V"]]
    U_curves <- specimen_data[grep("^U", names(specimen_data))]

    # Sort L/R curve names (L1, L2... R1, R2...)
    LR_names <- grep("^[LR][0-9]+$", names(specimen_data), value = TRUE)
    LR_names <- sort(LR_names)

    # --------------------------------------------------
    # Reorder each L/R curve if necessary
    # --------------------------------------------------
    LR_curves <- lapply(LR_names, function(LR_name) {
      curve <- specimen_data[[LR_name]]
      n <- nrow(curve)

      if (grepl("^R", LR_name)) {
        # For R curves, check orientation with respect to V
        closest_V_to_ldmN_idx <- find_closest_ventral_point_to_point(curve[n, ], V_curve)
        closest_V_to_ldmN <- V_curve[closest_V_to_ldmN_idx, ]
        dist_ldmN <- sqrt(sum((curve[n, ] - closest_V_to_ldmN)^2))
        dist_ldm1 <- sqrt(sum((curve[1, ] - closest_V_to_ldmN)^2))
        dist_ldm1_ldmN <- sqrt(sum((curve[1, ] - curve[n, ])^2))

        # Reverse if end landmark is farther from V than start
        if ((dist_ldmN > dist_ldm1) || (dist_ldmN > dist_ldm1_ldmN / 2)) {
          result <- reverse_curve(curve, specimen_id, LR_name, V_curve, updated_curve_info)
          curve <- result$curve
          updated_curve_info <- result$curve_info_ord
        }

      } else if (grepl("^L", LR_name)) {
        # For L curves, check orientation with respect to V
        closest_V_to_ldm1_idx <- find_closest_ventral_point_to_point(curve[1, ], V_curve)
        closest_V_to_ldm1 <- V_curve[closest_V_to_ldm1_idx, ]
        dist_ldm1 <- sqrt(sum((curve[1, ] - closest_V_to_ldm1)^2))
        dist_ldmN <- sqrt(sum((curve[n, ] - closest_V_to_ldm1)^2))
        dist_ldm1_ldmN <- sqrt(sum((curve[1, ] - curve[n, ])^2))

        # Reverse if start landmark is farther from V than end
        if ((dist_ldm1 > dist_ldmN) || (dist_ldm1 > dist_ldm1_ldmN / 2)) {
          result <- reverse_curve(curve, specimen_id, LR_name, V_curve, updated_curve_info)
          curve <- result$curve
          updated_curve_info <- result$curve_info_ord
        }
      }
      return(curve)
    })
    names(LR_curves) <- LR_names

    # --------------------------------------------------
    # Reorder V and U curves relative to L1/R1
    # --------------------------------------------------
    ap1 <- if ("L1" %in% names(LR_curves)) {
      LR_curves$L1[1, ]
    } else if ("R1" %in% names(LR_curves)) {
      LR_curves$R1[nrow(LR_curves$R1), ]
    } else {
      warning("No L1 or R1 found for specimen ", specimen_id, ". Skipping V/U reordering.")
      return(specimen_data)
    }

    # Reverse V if necessary (compare ap1 to start vs end of V)
    if (sum((ap1 - V_curve[nrow(V_curve), ])^2) < sum((ap1 - V_curve[1, ])^2)) {
      result <- reverse_curve(V_curve, specimen_id, "V", V_curve, updated_curve_info)
      V_curve <- result$curve
      updated_curve_info <- result$curve_info_ord
    }

    # Reverse U curves if necessary
    U_curves <- lapply(names(U_curves), function(U_name) {
      U_curve <- U_curves[[U_name]]

      if (sum((ap1 - U_curve[nrow(U_curve), ])^2) < sum((ap1 - U_curve[1, ])^2)) {
        result <- reverse_curve(U_curve, specimen_id, U_name, V_curve, updated_curve_info)
        U_curve <- result$curve
        updated_curve_info <- result$curve_info_ord
      }
      return(U_curve)
    })
    names(U_curves) <- names(specimen_data[grep("^U", names(specimen_data))])

    # Replace curves in specimen data with reordered versions
    specimen_data[names(LR_curves)] <- LR_curves
    specimen_data[["V"]] <- V_curve
    specimen_data[names(U_curves)] <- U_curves

    return(specimen_data)
  })

  names(landmarks_ord) <- names(data$landmarks)

  # --------------------------------------------------
  # Report which curves were reordered
  # --------------------------------------------------
  if (length(reordered_log) > 0) {
    reordered_df <- do.call(rbind, lapply(reordered_log, as.data.frame))
    reordered_df <- reordered_df[order(reordered_df$specimen, reordered_df$curve), ]

    # Regroup curves by specimen
    grouped <- split(reordered_df$curve, reordered_df$specimen)

    msg <- paste0(
      "The following curves were reordered:\n",
      paste0("- specimen ", names(grouped), " curve ", sapply(grouped, function(x) paste(x, collapse = ", ")), collapse = "\n")
    )
    message(msg)
  } else {
    message("No curves were reordered.")
  }

  # Return updated landmarks and curve_info
  return(list(landmarks = landmarks_ord, curve_info = updated_curve_info))
}


# --------------------------------------------------
#  Function: subset_data
# --------------------------------------------------

#' Subsets 3D aperture data and metadata
#'
#' @description
#' This function filters both the landmark data and associated metadata (`spec_tab`)
#' according to one or several user-defined criteria (e.g., species, genus, morphology).
#' It returns a reduced version of the input dataset containing only the selected specimens.
#'
#' @param data A list object containing at least:
#' \itemize{
#'   \item `landmarks`: a named list of specimen landmark data (each element named by specimen ID)
#'   \item `curve_info`: a data frame containing additional information about curves (e.g. curve names, indices)
#' }
#'
#' @param spec_tab A data frame containing specimen metadata, with at least one column
#' named `Specimen` (the unique specimen identifier used in `data$landmarks`).
#'
#' @param ... Named arguments defining filtering criteria.
#' Each argument name must correspond to a column name in `spec_tab`,
#' and its value must be a vector of accepted values.
#'
#' @details
#' The function filters the metadata table (`spec_tab`) according to the provided criteria,
#' then keeps only the specimens whose IDs are present in both the filtered metadata and
#' the landmark dataset (`data$landmarks`).
#'
#' If some filtered specimen IDs are not found in `data$landmarks`, they are silently removed.
#'
#' @return
#' A list containing two elements:
#' \itemize{
#'   \item `landmarks` - a subset of the original landmark list, containing only the selected specimens.
#'   \item `curve_info` - a filtered version of `data$curve_info` containing only the selected specimens.
#' }
#'
#' @examples
#' \dontrun{
#'
#' data("example_spec_tab", package = "AmmoniTools")
#' data("example_data", package = "AmmoniTools")
#'
#' # Example: keep only specimens of genus "Grammoceras" and species "P. tiziani"
#' subset_result <- subset_data(
#'   data = example_data,
#'   spec_tab = example_spec_tab,
#'   Genus = "Grammoceras",
#'   Species = "tiziani"
#' )
#'
#' # Example: keep only specimens of species "aalensis" and "opalinum"
#' subset_result <- subset_data(
#'   data = example_data,
#'   spec_tab = example_spec_tab,
#'   Species = c("aalensis", "opalinum")
#' )
#'
#' }
#'
#' @importFrom dplyr filter
#' @export

subset_data <- function(data, spec_tab, ...) {
  # ---- STEP 1: Collect filtering criteria ----
  # The user provides filters as named arguments, e.g. Genus = "Perisphinctes", Species = "tiziani"
  filters <- list(...)

  # ---- STEP 2: Apply filters sequentially to the specimen metadata ----
  filtered_spec_tab <- spec_tab
  for (filter_name in names(filters)) {
    # Keep only rows where the column value matches one of the filter values
    filtered_spec_tab <- filtered_spec_tab[filtered_spec_tab[[filter_name]] %in% filters[[filter_name]], ]
  }

  # ---- STEP 3: Identify which specimen IDs to keep ----
  kept_ids <- as.character(filtered_spec_tab$Specimen)

  # Handle potential mismatches between metadata and landmark data
  missing_ids <- setdiff(kept_ids, names(data$landmarks))
  kept_ids <- intersect(kept_ids, names(data$landmarks))

  # ---- STEP 4: Subset the landmark and curve info data ----
  data_subset <- list(
    landmarks = data$landmarks[kept_ids],
    curve_info = dplyr::filter(data$curve_info, Specimen %in% kept_ids)
  )

  # ---- STEP 5: Return the subsetted data ----
  return(data_subset)
}







