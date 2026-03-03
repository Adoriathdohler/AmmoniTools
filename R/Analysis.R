# --------------------------------------------------
# ANALYSIS
#
# Proc_peristome: Performs a Generalized Procrustes Analysis on 3D peristome curves
# compile_pca: Compiles PCA on aligned 3D peristomes curves
# plot_pca_shapes3d: Plots 3D peristome shape deformation along PCA axes
# compute_lda_accuracy: Compute Linear Discriminant Analysis classification accuracy across increasing numbers of PCA components, with optional cross-validation
# run_lda: Performs Linear Discriminant Analysis on 3D peristome PCA scores
# reconstruct_CAC_shape: Reconstructs 3D shapes along a CAC axis
# plot_cac_shapes3d: Plots 3D shapes along a CAC axis
# --------------------------------------------------


# --------------------------------------------------
#  Function: Proc_peristome
# --------------------------------------------------
#' Performs a Generalized Procrustes Analysis on 3D peristome curves
#'
#' @description
#' This function performs a Generalized Procrustes Analysis (GPA) on a set of
#' apertural curves (e.g., ammonite peristomes) extracted from multiple specimens
#'  using \code{\link[geomorph]{gpagen}}.It aligns all curves to a common coordinate system,
#'  optionally allowing sliding landmarks.
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
#' @param ProcD 	A logical value indicating whether or not Procrustes distance
#' should be used as the criterion for optimizing the positions of semilandmarks
#' (if not, bending energy is used). Parameter from \code{\link[geomorph]{gpagen}}.
#'
#' @details
#' The function first extracts all apertural curves from each specimen, stores them in
#' a unified 3D array `(landmarks x dimensions x apertures)`, and then performs
#' a Procrustes superimposition using \code{\link[geomorph]{gpagen}}.
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
#'   \item and other standard outputs from \code{\link[geomorph]{gpagen}}.
#' }
#'
#' @examples
#' \dontrun{
#' gpa <- Proc_peristome(reconstructed_data, slide = TRUE, fixed_ldm = c(1, 25, 49))
#' summary(gpa)
#' plot(gpa)
#' }
#' @seealso \link[geomorph]{gpagen}
#' @importFrom geomorph gpagen
#' @importFrom stats complete.cases
#' @export

Proc_peristome <- function(data, slide = FALSE, fixed_ldm = c(), ProcD = FALSE) {

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
    gpa_result <- gpagen(apertures_array, PrinAxes = TRUE, curves = curve_matrix, ProcD = ProcD)
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
#  Function: compile_pca
# --------------------------------------------------

#' Compile PCA on aligned 3D peristome curves
#'
#' @description
#' Performs a principal component analysis (PCA) from GPA-aligned
#' landmark coordinates using \code{geomorph::gm.prcomp}.
#' The function optionally removes components explaining less than 1 percent
#' of total variance, adds centroid size, and merges specimen metadata.
#'
#' @param gpa A GPA object produced by \code{geomorph::gpagen}
#'   containing aligned coordinates in \code{gpa$coords}
#'   and centroid sizes in \code{gpa$Csize}.
#' @param data A list containing reconstructed landmark configurations.
#' @param spec_tab A data frame containing specimen metadata.
#'   Must include a column named \code{Specimen}.
#' @param filter_pc01 Logical. If TRUE, principal components explaining
#'   less than 1 percent of total variance are removed.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{pca_scores}{Data frame of PCA scores with centroid size and metadata.}
#'     \item{variance}{Summary of explained variance.}
#'     \item{pca}{The full \code{gm.prcomp} object.}
#'   }
#'
#' @examples
#' \dontrun{
#' res <- compile_pca(gpa, data, spec_tab, filter_pc01 = TRUE)
#' }
#'
#' @importFrom geomorph gm.prcomp
#' @importFrom dplyr left_join
#' @export

compile_pca <- function(gpa, data, spec_tab, filter_pc01 = FALSE) {

  # ---------------------------------------------------------------
  # 1. Run PCA
  # ---------------------------------------------------------------
  pca <- geomorph::gm.prcomp(gpa$coords)
  pca_scores <- as.data.frame(pca$x)

  # Capture full PCA variance
  summary_pca <- summary(pca)
  invisible(summary_pca)

  # ---------------------------------------------------------------
  # 2. Expand Specimen ID for each aperture
  # ---------------------------------------------------------------
  specimen_vec <- character()

  for (specimen_id in names(data$landmarks)) {
    aperture_names <- grep("^\\d+$",
                           names(data$landmarks[[specimen_id]]),
                           value = TRUE)
    if (length(aperture_names) > 0) {
      specimen_vec <- c(specimen_vec,
                        rep(specimen_id, length(aperture_names)))
    }
  }

  pca_scores$Specimen <- as.numeric(specimen_vec)
  pca_scores$Csize <- gpa$Csize

  # ---------------------------------------------------------------
  # 3. Optionally filter PCs
  # ---------------------------------------------------------------
  # with variance < 1%
  if (filter_pc01) {

    # variance explained (%)
    var_expl <- summary_pca$PC.summary[2, ]

    # PCs above threshold
    keep_pc <- names(var_expl)[var_expl > 0.01]

    comp_cols <- intersect(colnames(pca_scores), keep_pc)
    other_cols <- c("Specimen", "Csize")

    pca_scores <- pca_scores[, c(comp_cols, other_cols), drop = FALSE]

    # Also filter the variance object
    variance <- summary_pca
    variance$PC.summary <- variance$PC.summary[, keep_pc, drop = FALSE]

  } else {
    variance <- summary_pca
  }

  # ---------------------------------------------------------------
  # 4. Merge with specimen metadata
  # ---------------------------------------------------------------
  pca_scores <- dplyr::left_join(pca_scores,
                                 spec_tab,
                                 by = "Specimen")

  # ---------------------------------------------------------------
  # 5. Return
  # ---------------------------------------------------------------
  list(
    pca_scores = pca_scores,
    variance = variance,
    pca = pca
  )
}



# --------------------------------------------------
#  Function: plot_pca_shapes3d
# --------------------------------------------------

#' Plots 3D peristome shape deformation along PCA axes
#'
#' @description
#' This function visualizes 3D shape deformations associated with selected
#' principal components from a geometric morphometric analysis.
#' For each requested PCA axis, three warped shapes are displayed in 3D space:
#' the minimal deformation, the consensus shape, and the maximal deformation
#' relative to the consensus configuration.
#'
#' Shape visualization relies internally on
#' \link[geomorph]{plotRefToTarget}, which displays landmark-wise displacements
#' from a reference shape (the GPA consensus) to a target shape using
#' point-based deformation plots in an \code{rgl} scene.
#'
#' The resulting visualization is arranged as a grid with one row per principal
#' component and three columns corresponding to the \emph{min}, \emph{mean}, and
#' \emph{max} warped shapes. The output is returned as an interactive
#' \code{rglwidget}, suitable for use in RMarkdown or Jupyter environments.
#'
#' @param gpa A GPA object obtained from \link[geomorph]{gpagen}, containing
#' at least the consensus shape (\code{gpa$consensus}).
#' @param pca A PCA object obtained from \link[geomorph]{gm.prcomp}, containing
#' warped shapes in \code{pca$shapes}.
#' @param pcs Numeric vector of principal component indices to visualize
#' (default: \code{1:3}).
#' @param min_col Color used for the minimal deformation shape.
#' @param mid_col Color used for the consensus deformation shape
#' @param max_col Color used for the maximal deformation shape.
#' @param size_min Point size for landmarks in the minimal deformation plot.
#' @param size_mid Point size for landmarks in the consensus deformation plot.
#' @param size_max Point size for landmarks in the maximal deformation plot.
#' @param width Width of the returned \code{rglwidget}.
#' @param height Height of the returned \code{rglwidget}.
#'
#' @details
#' For each selected principal component:
#' \enumerate{
#' \item The minimum and maximum warped shapes are extracted from
#' \code{pca$shapes[[pc]]}.
#' \item An consensus shape is computed as the arithmetic mean of the
#' minimum and maximum shapes.
#' \item Each shape is visualized relative to the GPA consensus using
#' \link[geomorph]{plotRefToTarget} with \code{method = "points"}.
#' }
#'
#' No interpolation along the PC axis is performed beyond the min–mean–max
#' configurations provided by \link[geomorph]{gm.prcomp}.
#'
#' @return
#' An interactive \code{rglwidget} displaying a grid of 3D deformation plots.
#'
#' @seealso \link[geomorph]{gm.prcomp}
#' @seealso \link[geomorph]{gpagen}
#' @seealso \link[geomorph]{plotRefToTarget}
#'
#' @importFrom geomorph gridPar
#' @importFrom geomorph plotRefToTarget
#' @export
plot_pca_shapes3d <- function(gpa, pca, pcs = 1:3,
                              min_col = "blue", max_col = "red", mid_col = "black",
                              size_min = 1.5, size_mid = 1.5, size_max = 1.5,
                              width = 600, height = 600) {

  # --------------------------
  # Grid parameters (colors)
  # --------------------------
  gp_min <- gridPar(pt.bg = "gray", pt.size = 0,
                    tar.pt.bg = min_col, tar.pt.size = size_min)

  gp_mid <- gridPar(pt.bg = "gray", pt.size = 0,
                    tar.pt.bg = mid_col, tar.pt.size = size_mid)

  gp_max <- gridPar(pt.bg = "gray", pt.size = 0,
                    tar.pt.bg = max_col, tar.pt.size = size_max)

  # --------------------------
  # Open 3D viewer
  # --------------------------
  rgl::open3d(silent = TRUE)
  par3d(mfrow3d(length(pcs), 3))

  # --------------------------
  # Loop over selected PCs
  # --------------------------
  for (pc in pcs) {

    # Shapes
    shp_min <- pca$shapes[[pc]]$min
    shp_max <- pca$shapes[[pc]]$max
    shp_mid <- (shp_min + shp_max) / 2   # consensus shape

    # Plot min
    plotRefToTarget(gpa$consensus, shp_min, method = "points", gridPars = gp_min)

    # Plot mid
    plotRefToTarget(gpa$consensus, shp_mid, method = "points", gridPars = gp_mid)

    # Plot max
    plotRefToTarget(gpa$consensus, shp_max, method = "points", gridPars = gp_max)
  }

  # --------------------------
  # Return widget for RMarkdown / Jupyter
  # --------------------------
  rgl::rglwidget(width = width, height = height)
}


# --------------------------------------------------
# Function: reconstruct_CAC_shape
# --------------------------------------------------

#'
#' Reconstructs theoretical 3D shapes along a Centroid Axis Component (CAC)
#'
#' @description
#' This function reconstructs theoretical 3D landmark configurations
#' corresponding to a given Centroid Axis Component (CAC), as defined by
#' \link[Morpho]{CAC}. For the selected CAC axis (by default the first),
#' three shapes are generated: the minimal, mean, and maximal shapes observed
#' along that axis.
#'
#' Shape reconstruction is performed in the original 3D landmark space by
#' combining the GPA consensus shape with a linear deformation derived from
#' PCA loadings and CAC coefficients.
#'
#' @param gpa A \code{geomorph} GPA object obtained from
#'   \link[geomorph]{gpagen}, containing the aligned coordinates and the
#'   consensus shape (\code{gpa$consensus}).
#' @param pca A PCA object (typically from \link[geomorph]{gm.prcomp})
#'   providing the rotation/loadings matrix used as input for the CAC analysis.
#' @param cac_res Output of \link[Morpho]{CAC}, containing CAC loadings
#'   (\code{cac_res$CAC}) and specimen scores (\code{cac_res$CACscores}).
#'
#' @return
#' A named list of 3D landmark configurations (matrices with dimensions
#' \code{n_landmarks x 3}) containing:
#' \describe{
#'   \item{min}{Shape corresponding to the minimal observed CAC score.}
#'   \item{mean}{Shape corresponding to the mean CAC score
#'               (typically the consensus shape, CAC score = 0).}
#'   \item{max}{Shape corresponding to the maximal observed CAC score.}
#' }
#'
#' @details
#' The reconstruction follows a linear model:
#' \deqn{
#'   X_{CAC} = X_{consensus} + s \cdot (P \times c)
#' }
#' where:
#' \itemize{
#'   \item \eqn{X_{consensus}} is the GPA consensus shape,
#'   \item \eqn{P} is the matrix of PCA loadings used for the CAC analysis,
#'   \item \eqn{c} is the vector of CAC loadings for the selected axis,
#'   \item \eqn{s} is the chosen CAC score (minimum, mean = 0, or maximum).
#' }
#'
#' The minimal and maximal shapes correspond to the extreme CAC scores observed
#' in the dataset, while the mean shape corresponds to a CAC score of zero.
#'
#' This function is primarily intended for visualization and interpretation
#' of CAC axes, for example in combination with 3D plotting utilities.
#'
#' @seealso
#' \link[Morpho]{CAC},
#' \link[geomorph]{gm.prcomp},
#' \link[geomorph]{gpagen}
#'
#' @export
reconstruct_CAC_shape <- function(gpa, pca, cac_res) {

  cac_axis <- 1                       # default: first CAC axis
  consensus <- gpa$consensus
  n_landmarks <- nrow(consensus)

  # --------------------------
  # Extract PCs used for CAC
  # --------------------------
  pc_indices <- 1:length(cac_res$CAC[,1])   # row indices of CAC loadings
  pcs <- pca$rotation[, pc_indices]         # n_landmarks*3 x num_PCs

  # --------------------------
  # Extract CAC loadings for this axis
  # --------------------------
  cac_loadings <- cac_res$CAC[, cac_axis]   # vector of loadings

  # --------------------------
  # Theoretical CAC scores
  # --------------------------
  score_min  <- min(cac_res$CACscores[, cac_axis])
  score_max  <- max(cac_res$CACscores[, cac_axis])
  score_mean <- 0                            # CAC scores are centered

  # --------------------------
  # Internal function: reconstruct shape given a CAC score
  # --------------------------
  reconstruct <- function(score) {
    # Linear combination in PC space
    delta_vec <- pcs %*% cac_loadings        # vector of length n_landmarks*3
    delta <- matrix(delta_vec * score, ncol = 3, byrow = TRUE)

    consensus + delta
  }

  # Return min, mean, max shapes
  list(
    min  = reconstruct(score_min),
    mean = reconstruct(score_mean),
    max  = reconstruct(score_max)
  )
}


# --------------------------------------------------
# Function: plot_cac_shapes3d
# --------------------------------------------------

#' Plots 3D shapes along a CAC axis
#'
#' This function visualizes theoretical 3D shapes along a specified CAC axis,
#' showing minimal, maximal, and mean shapes in a 3D rgl viewer.
#'
#' @param gpa Geomorph GPA object containing the consensus shape.
#' @param pca PCA object used to compute CAC.
#' @param cac_res Output of \code{CAC()} containing loadings and scores.
#' @param min_col Color for the minimal shape (default = "blue").
#' @param max_col Color for the maximal shape (default = "red").
#' @param mid_col Color for the mean shape (default = "black").
#' @param size_min, size_mid, size_max Point sizes for min/mean/max shapes.
#' @param size_mid Point size for landmarks in the consensus deformation plot.
#' @param size_max Point size for landmarks in the maximal deformation plot.
#' @param width,height Width and height (pixels) of the rgl 3D widget.
#'
#' @return An interactive rgl 3D widget displaying the shapes.
#'
#' @importFrom dplyr left_join
#' @importFrom geomorph gridPar
#' @importFrom geomorph plotRefToTarget
#' @export
plot_cac_shapes3d <- function(gpa, pca, cac_res,
                              min_col = "blue", max_col = "red", mid_col = "black",
                              size_min = 1.5, size_mid = 1.5, size_max = 1.5,
                              width = 600, height = 600) {

  # --------------------------
  # Reconstruct shapes along CAC axis
  # --------------------------
  shapes <- reconstruct_CAC_shape(gpa, pca, cac_res)

  # --------------------------
  # Set grid parameters for coloring
  # --------------------------
  gp_min <- gridPar(pt.bg = "gray", pt.size = 0,
                    tar.pt.bg = min_col, tar.pt.size = size_min)
  gp_mid <- gridPar(pt.bg = "gray", pt.size = 0,
                    tar.pt.bg = mid_col, tar.pt.size = size_mid)
  gp_max <- gridPar(pt.bg = "gray", pt.size = 0,
                    tar.pt.bg = max_col, tar.pt.size = size_max)

  # --------------------------
  # Open 3D viewer
  # --------------------------
  rgl::open3d(silent = TRUE)
  par3d(mfrow3d(1, 3))   # layout: 1 row, 3 columns

  # --------------------------
  # Plot minimal, mean, maximal shapes
  # --------------------------
  plotRefToTarget(gpa$consensus, shapes$min,  method = "points", gridPars = gp_min)
  plotRefToTarget(gpa$consensus, shapes$mean, method = "points", gridPars = gp_mid)
  plotRefToTarget(gpa$consensus, shapes$max,  method = "points", gridPars = gp_max)

  # --------------------------
  # Return 3D widget
  # --------------------------
  rgl::rglwidget(width = width, height = height)
}


# --------------------------------------------------
# Function: compile_cac
# --------------------------------------------------
' Compile CAC Scores from Procrustes Coordinates
#'
#' Computes the Corrected Aperture Contribution (CAC) for a set of specimens based on
#' Procrustes-aligned landmark coordinates, and returns a data frame of CAC scores
#' combined with specimen metadata.
#'
#' @param gpa An object returned by \code{geomorph::gpagen()}, containing Procrustes-aligned coordinates (\code{coords}) and centroid sizes (\code{Csize}).
#' @param data A list containing specimen landmarks, with structure \code{data$landmarks}.
#' @param spec_tab A data frame containing specimen metadata. Must include a column \code{Specimen} matching the specimen IDs in \code{data$landmarks}.
#' @param log Logical. If \code{TRUE}, logging messages from the CAC computation are printed. Default is \code{FALSE}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{cac_scores}}{A data frame containing CAC scores (\code{CACscore}), optionally RSC scores (\code{RSC1}, \code{RSC2}, …), specimen IDs, centroid sizes, and metadata from \code{spec_tab}.}
#'   \item{\code{cac}}{The full CAC object returned by \code{CAC()}, containing all intermediate results.}
#' }
#'
#' @details
#' The function first vectorizes the Procrustes coordinates using \code{geomorph::two.d.array}.
#' CAC scores are computed via the \code{CAC()} function (external, must be available in the environment).
#' Metadata from \code{spec_tab} is then joined to the CAC scores based on specimen IDs.
#'
#' @examples
#' \dontrun{
#' library(geomorph)
#' data("example_data", package = "AmmoniTools")
#' data("example_spec_tab", package = "AmmoniTools")
#'
#' # GPA alignment
#' gpa <- gpagen(example_data$coords)
#'
#' # Compute CAC scores
#' results <- compile_cac(gpa, data = example_data, spec_tab = example_spec_tab)
#'
#' # View scores
#' head(results$cac_scores)
#' }
#'
#' @importFrom Morpho CAC
#' @importFrom geomorph two.d.array
#' @importFrom dplyr left_join
#' @export
compile_cac <- function(gpa, data, spec_tab, log = FALSE) {

  # --------------------------------------------------
  # 1. Vectorize Procrustes coordinates
  # --------------------------------------------------
  coords_mat <- geomorph::two.d.array(gpa$coords)

  # --------------------------------------------------
  # 2. Compute CAC on Procrustes coordinates
  # --------------------------------------------------
  cac <- CAC(coords_mat, gpa$Csize, log = log)

  # --------------------------------------------------
  # 3. Build CAC score dataframe
  # --------------------------------------------------
  cac_scores <- as.data.frame(cac$CACscore)
  colnames(cac_scores) <- "CACscore"

  # Add RSC scores if present
  if (!is.null(cac$RSCscores)) {
    rsc_df <- as.data.frame(cac$RSCscores)
    colnames(rsc_df) <- paste0("RSC", seq_len(ncol(rsc_df)))
    cac_scores <- cbind(cac_scores, rsc_df)
  }

  # --------------------------------------------------
  # 4. Rebuild Specimen vector (same logic as PCA)
  # --------------------------------------------------
  specimen_vec <- character()

  for (specimen_id in names(data$landmarks)) {
    aperture_names <- grep("^\\d+$",
                           names(data$landmarks[[specimen_id]]),
                           value = TRUE)
    if (length(aperture_names) > 0) {
      specimen_vec <- c(specimen_vec,
                        rep(specimen_id, length(aperture_names)))
    }
  }

  cac_scores$Specimen <- as.numeric(specimen_vec)
  cac_scores$Csize <- gpa$Csize

  # --------------------------------------------------
  # 5. Join metadata
  # --------------------------------------------------
  cac_scores <- dplyr::left_join(cac_scores,
                                 spec_tab,
                                 by = "Specimen")

  # --------------------------------------------------
  # 6. Output
  # --------------------------------------------------
  list(
    cac_scores = cac_scores,
    cac = cac
  )
}

# --------------------------------------------------
# Function: compute_lda_accuracy
# --------------------------------------------------
#' Compute Linear Discriminant Analysis classification accuracy across increasing numbers of PCA components, with optional cross-validation
#'
#' This function evaluates the classification performance of a Linear
#' Discriminant Analysis (LDA) applied to PCA scores by progressively
#' increasing the number of principal components used as predictors.
#' Classification accuracy is computed either using leave-one-out
#' cross-validation (default) or apparent (training) accuracy.
#'
#' The function is particularly useful for assessing how many principal
#' components are required to optimally discriminate groups (e.g. genera,
#' species) in morphometric datasets.
#'
#' @param pca_scores A data frame containing PCA scores. Principal components
#'   must be named using the prefix `"Comp"` (e.g. `"Comp1"`, `"Comp2"`, ...),
#'   and a grouping variable must be present.
#' @param max_pc Integer specifying the maximum number of principal components
#'   to include in the LDA. Defaults to 30. If larger than the number of
#'   available components, all components are used.
#' @param group_var Character string giving the name of the grouping variable
#'   used for LDA classification (e.g. `"Genus"` or `"Species"`).
#' @param cv Logical indicating whether to compute classification accuracy
#'   using leave-one-out cross-validation (`TRUE`, default) or apparent
#'   (training) accuracy (`FALSE`).
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{n_pc}{Number of principal components used in the LDA.}
#'   \item{accuracy}{Classification accuracy for the corresponding number of PCs.}
#' }
#'
#' @details
#' Rows containing missing values are removed prior to analysis.
#' If fewer than two groups are present for a given number of PCs,
#' the accuracy is returned as `NA`.
#'
#' @examples
#' \dontrun{
#' # Assuming pca_scores is a data frame containing Comp1, Comp2, ..., and Genus
#' acc <- compute_lda_accuracy(
#'   pca_scores = pca_scores,
#'   max_pc = 10,
#'   group_var = "Genus",
#'   cv = TRUE
#' )
#'
#' plot(acc$n_pc, acc$accuracy, type = "b",
#'      xlab = "Number of PCs", ylab = "LDA accuracy")
#' }
#'
#' @seealso
#' MASS::lda
#' @importFrom MASS lda
#' @export
compute_lda_accuracy <- function(pca_scores,
                                 max_pc = 30,
                                 group_var = "Genus",
                                 cv = TRUE) {

  # -----------------------------
  # 0. Checks
  # -----------------------------
  if (!group_var %in% names(pca_scores)) {
    stop(paste0("Grouping variable '", group_var, "' not found in pca_scores"))
  }

  pca_vars_all <- grep("^Comp", names(pca_scores), value = TRUE)

  if (length(pca_vars_all) < max_pc) {
    warning("max_pc larger than available PCs; using all available PCs")
    max_pc <- length(pca_vars_all)
  }

  acc <- data.frame(
    n_pc = integer(),
    accuracy = numeric()
  )

  # -----------------------------
  # 1. Loop over number of PCs
  # -----------------------------
  for (n_pc in seq_len(max_pc)) {

    pca_vars <- pca_vars_all[1:n_pc]

    lda_data <- pca_scores[, c(pca_vars, group_var)]
    lda_data <- stats::na.omit(lda_data)

    lda_data[[group_var]] <- as.factor(lda_data[[group_var]])

    # Skip if only one group (LDA impossible)
    if (nlevels(lda_data[[group_var]]) < 2) {
      acc <- rbind(
        acc,
        data.frame(n_pc = n_pc, accuracy = NA_real_)
      )
      next
    }

    # Build formula dynamically
    lda_formula <- stats::as.formula(
      paste(group_var, "~ .")
    )

    if (cv) {
      # Leave-One-Out Cross-Validation
      lda_fit <- MASS::lda(lda_formula, data = lda_data, CV = TRUE)
      acc_value <- mean(lda_fit$class == lda_data[[group_var]])
    } else {
      # Apparent (training) accuracy
      lda_fit <- lda(lda_formula, data = lda_data)
      pred <- predict(lda_fit)$class
      acc_value <- mean(pred == lda_data[[group_var]])
    }

    acc <- rbind(
      acc,
      data.frame(
        n_pc = n_pc,
        accuracy = acc_value
      )
    )
  }

  return(acc)
}

# --------------------------------------------------
# Function: run_lda
# --------------------------------------------------
#' Performs Linear Discriminant Analysis on 3D peristome PCA scores
#'
#' This function runs a Linear Discriminant Analysis (LDA) using a specified
#' number of principal components derived from a PCA. It returns the scores
#' along the first two linear discriminant axes together with the original
#' grouping variable.
#'
#' The function is intended for visualizing and analysing group separation
#' (e.g. among genera or species) in morphometric datasets after dimensionality
#' reduction by PCA.
#'
#' @param pca_scores A data frame containing PCA scores. Principal components
#'   must be named using the prefix `"Comp"` (e.g. `"Comp1"`, `"Comp2"`, ...),
#'   and must include a grouping variable.
#' @param n_pc Integer specifying the number of principal components to include
#'   in the LDA. Defaults to 30. If larger than the number of available PCs,
#'   all available components are used.
#' @param group_var Character string giving the name of the grouping variable
#'   used for LDA classification (e.g. `"Genus"` or `"Species"`).
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{LD1}{Scores along the first linear discriminant axis.}
#'   \item{LD2}{Scores along the second linear discriminant axis.}
#'   \item{group_var}{Grouping variable used for the LDA (factor).}
#' }
#'
#' @details
#' Rows containing missing values are removed prior to analysis.
#' The function requires at least two levels in the grouping variable.
#' Only the first two linear discriminant axes are returned.
#'
#' @examples
#' \dontrun{
#' # Assuming pca_scores is a data frame containing Comp1, Comp2, ..., and Genus
#' lda_scores <- run_lda(
#'   pca_scores = pca_scores,
#'   n_pc = 15,
#'   group_var = "Genus"
#' )
#'
#' plot(lda_scores$LD1, lda_scores$LD2,
#'      col = lda_scores$Genus,
#'      xlab = "LD1", ylab = "LD2")
#' }
#'
#' @seealso
#' MASS::lda
#'
#' @importFrom MASS lda
#' @export
run_lda <- function(pca_scores,
                    n_pc = 30,
                    group_var = "Genus") {

  # -----------------------------
  # 0. Checks
  # -----------------------------
  if (!group_var %in% names(pca_scores)) {
    stop(paste0("Grouping variable '", group_var, "' not found in pca_scores"))
  }

  pca_vars_all <- grep("^Comp", names(pca_scores), value = TRUE)

  if (length(pca_vars_all) < n_pc) {
    warning("n_pc larger than available PCs; using all available PCs")
    n_pc <- length(pca_vars_all)
  }

  pca_vars <- pca_vars_all[1:n_pc]

  lda_data <- pca_scores[, c(pca_vars, group_var)]
  lda_data <- stats::na.omit(lda_data)

  lda_data[[group_var]] <- as.factor(lda_data[[group_var]])

  # LDA requires at least 2 groups
  if (nlevels(lda_data[[group_var]]) < 2) {
    stop("LDA requires at least two groups in grouping variable")
  }

  # -----------------------------
  # 1. Run LDA
  # -----------------------------
  lda_formula <- stats::as.formula(
    paste(group_var, "~ .")
  )

  lda_model <- MASS::lda(lda_formula, data = lda_data)
  lda_pred  <- predict(lda_model)

  # -----------------------------
  # 2. Output
  # -----------------------------
  out <- data.frame(
    LD1 = lda_pred$x[, 1],
    LD2 = lda_pred$x[, 2],
    Group = lda_data[[group_var]]
  )

  names(out)[3] <- group_var  # keep original name

  return(out)
}
