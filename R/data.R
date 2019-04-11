#' The Ensembl dictionary.
#'
#' A dataframe containing the Ensembl 96 release dictionary. Updated at 2019-04-11.
#'
#' @format A data frame with 184 rows and 8 variables:
#' \describe{
#'   \item{Prefix}{Ensembl ID prefix for each available species.}
#'   \item{Species.Full.name}{Species scientific and common names, when available.}
#'   \item{Species.Short.name}{Species scientific binomial names, separated by spaces.}
#'   \item{Species.fasta}{Species scientific binomial names, separated by underlines, appropriate for saving as fasta identifier.}
#'   \item{Common}{Common species names.}
#'   \item{dataset}{Ensembl dataset names, used as an input for biomaRt.}
#'   \item{description}{Brief description of the dataset.}
#'   \item{version}{Dataset/genome version.}
#' }
#' @source \url{https://www.ensembl.org/info/genome/stable_ids/prefixes.html}
"dict"
