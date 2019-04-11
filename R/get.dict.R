#' Get the Ensembl dictionary
#'
#' This function retrieves our creates an Ensembl dictionary, which includes species names in different formats, as well as IDs prefixes and dataset names for all available species names.
#'
#' @param url A character containing the URL of the official Ensembl site containing the prefix table. Defaults to the European mirror (www.ensembl.org/info/genome/stable_ids/prefixes.html)
#' @param withdatasets A logical indicating whether it should query Ensembl for dataset names. This is necessary for \code{get.orthoIDs} and \code{get.seq} functions. Defaults to TRUE.
#' @param host Ensembl host. Defaults to the European mirror.
#'
#'
#' @return A data frame containing the Ensembl dictionary.
#' @export
#'
#' @examples
#' get.dict()
#'
get.dict <- function(url = "https://www.ensembl.org/info/genome/stable_ids/prefixes.html", withdatasets = TRUE, host = "www.ensembl.org"){
#  library(rvest)
#  library(biomaRt)
#  library(dplyr)
  # Scrap data from the website
  webpage <- xml2::read_html(url)
  webtable <- webpage %>% rvest::html_nodes("tbody")
  contents <- webtable %>% rvest::html_nodes("td") %>% rvest::html_text(trim = TRUE)
  Prefix <- contents[1:length(contents) %% 2 == 1]
  Species <- contents[1:length(contents) %% 2 == 0]
  prefix_table <- data.frame(Prefix = Prefix, Species.Full.name = Species, stringsAsFactors = F)
  # This list does not include D. melanogaster, we manually include it
  Dmel <- c("FBgn", "Drosophila melanogaster (Vinegar fly)")
  prefix_table <- rbind(prefix_table, Dmel)
  # Check for duplicated identifiers
  # prefix_table$Prefix[duplicated(prefix_table$Prefix)]
  # Correcting duplicate identifier situation and additional mismatches
  prefix_table[grep("Homo sapiens", prefix_table$Species.Full.name),]$Prefix <- "ENS.0"
  prefix_table[grep("Dingo", prefix_table$Species.Full.name),]$Prefix <- "ENSCAF.0002"
  prefix_table[grep("Dog", prefix_table$Species.Full.name),]$Prefix <- "ENSCAF.0000"
  prefix_table[grep("CriGri", prefix_table$Species.Full.name),]$Prefix <- "ENSCGR.000000"
  prefix_table[grep("CHOK1GS", prefix_table$Species.Full.name),]$Prefix <- "ENSCGR.00001"
  prefix_table[grep("PICR", prefix_table$Species.Full.name),]$Prefix <-"ENSCGR.00015"
  prefix_table[grep("HdrR", prefix_table$Species.Full.name),]$Prefix <-"ENSORL.00000"
  prefix_table[grep("HNI", prefix_table$Species.Full.name),]$Prefix <-"ENSORL.0002"
  prefix_table[grep("HSOK", prefix_table$Species.Full.name),]$Prefix <-"ENSORL.00015"
  # We add an additional column with underline-separated species names
  prefix_table$Species.Short.name <- gsub("^(.+) \\(.*\\)", "\\1", prefix_table$Species.Full.name)
  prefix_table$Species.fasta <- gsub(" ", "_", prefix_table$Species.Short.name)

  if (withdatasets){
  # We retrieve the list of Ensembl datasets to have a full prefix-dataset-species dataset
  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = host, ensemblRedirect = F)
  DatasetList <- biomaRt::listDatasets(mart)
  # We make an additional column with common names, which we'll use for joining both datasets
  DatasetList$Common <- gsub("(.*) genes.*", "\\1", DatasetList$description)
  prefix_table$Common <- gsub(".* \\((.*)\\)", "\\1", prefix_table$Species.Full.name)
  # We correct some mismatches in common names in the Datasets dataset
  DatasetList$Common <- gsub("C.savignyi", "Ciona savignyi", DatasetList$Common)
  DatasetList$Common <- gsub("C.intestinalis", "Ciona intestinalis", DatasetList$Common)
  DatasetList$Common <- gsub("Drosophila melanogaster", "Vinegar fly", DatasetList$Common)
  DatasetList$Common <- as.character(DatasetList$Common)
  prefix_table <- dplyr::right_join(prefix_table, DatasetList, by = "Common")
  }
  return(prefix_table)
}
