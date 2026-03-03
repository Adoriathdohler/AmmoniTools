
# --------------------------------------------------
#  Data: example_spec_tab
# --------------------------------------------------

#' Example specimen metadata associated with example_data
#'
#' @description
#' `example_spec_tab` contains specimen-level metadata corresponding to
#' the landmark dataset \code{example_data}.
#'
#' The dataset includes taxonomic, stratigraphic, geographic and acquisition
#' information for 162 ammonite specimens belonging to Hildoceratoidea.
#'
#' This table is primarily intended to:
#'   - provide taxonomic grouping variables for morphometric analyses,
#'   - allow filtering by genus, species, family or stratigraphic age,
#'   - support hierarchical missing-landmark estimation procedures,
#'   - link biological interpretation to geometric morphometric results.
#'
#' @details
#' Each row corresponds to a single specimen identified by a unique
#' specimen ID matching the identifiers used in \code{example_data$landmarks}
#' and \code{example_data$curve_info}.
#'
#' The dataset contains 20 variables including:
#'
#' \itemize{
#'   \item \code{Specimen}: unique specimen identifier
#'   \item \code{Collection}: holding institution or source
#'   \item \code{Collection_code}: institutional catalog number
#'   \item \code{Family}, \code{Genus}, \code{Subgenus}, \code{Species}: taxonomic information
#'   \item \code{Age}: stratigraphic stage (e.g., Toarcian, Aalenian)
#'   \item \code{Locality}: geographic origin
#'   \item \code{Collector}: specimen collector
#'   \item \code{Diameter_mm}: maximum shell diameter (in mm)
#'   \item \code{Imaging_technic}: 3D acquisition method
#'   \item \code{3Dreconstruction_software}: software used for model reconstruction
#'   \item \code{Landmarkings_credit}: person responsible for landmarking
#' }
#'
#' @format
#' A tibble with 162 rows and 20 columns.
#'
#' @examples
#' data("example_spec_tab", package = "AmmoniTools")
#' str(example_spec_tab)

"example_spec_tab"

# --------------------------------------------------
#  Data: example_data
# --------------------------------------------------

#' Example dataset of peristome landmark curves
#'
#' @description
#' `example_data` is a structured example dataset containing raw 3D landmark
#' coordinates extracted from ammonite peristomes.
#'
#' The dataset includes peristomes from 162 specimens of Hildoceratoidea and
#' is intended for testing and demonstration of functions implemented in
#' the package (e.g. landmark extraction, interpolation, missing landmark
#' estimation, and morphometric analyses).
#'
#' The object is a list with the following main components:
#'
#' • `landmarks`: a nested list of specimens.
#'   Each specimen contains multiple landmark curves:
#'     - L: left lateral curves
#'     - R: right lateral curves
#'     - U/Ubis: umbilical curves
#'     - V: ventral curve
#'
#'   Each curve is stored as a data.frame of 3D coordinates (X, Y, Z).
#'   L and R curves contain 200 equidistant points.
#'   U, Ubis and V curves contain 100 equidistant points.
#'
#' • `curve_info`: a data.frame describing each curve with:
#'     - `Specimen`: specimen identifier
#'     - `Curve`: curve name (e.g. L1, L2, R1, R2, U, V)
#'     - `Type`: curve category (L, R, U, V)
#'     - `missing_ldm`: list column indicating indices of missing landmarks
#'
#' • `example_data_tab`: specimen-level metadata (taxonomy and identifiers).
#'
#' More information on the specimens is available in `example_spec_tab`.
#'
#' @details
#' Coordinates are expressed in 3D Cartesian space and correspond to
#' semi-landmarks sampled along each peristome curve after applying
#' \code{extract_landmarks()}.
#'
#' @format
#' A list containing nested landmark curves and associated metadata.
#'
#' @examples
#' data("example_data", package = "AmmoniTools")
#' str(example_data)

"example_data"
