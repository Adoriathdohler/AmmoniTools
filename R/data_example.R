
# --------------------------------------------------
#  Data: example_spec_tab
# --------------------------------------------------
#' Example Specimen Data
#'
#' @description
#' A small subset of specimen metadata used for testing and demonstration.
#'
#' @format A data frame with specimen metadata, including columns such as:
#'   \itemize{
#'     \item \code{N}: Specimen identifier
#'     \item \code{Genus}: Genus name
#'     \item \code{Species}: Species name
#'   }
#' @examples
#' data("example_spec_tab", package = "AmmoniTools")
#' str(example_spec_tab)

"example_spec_tab"

# --------------------------------------------------
#  Data: example_data
# --------------------------------------------------
#' Example Landmark Data
#'
#' @description
#' A small subset of raw landmark data for testing and demonstration
#' (data after using \code{extract_landmarks} function)
#'
#' @format A list of specimen landmark curves.
#' @examples
#' data("example_data", package = "AmmoniTools")
#' str(example_data)

"example_data"
