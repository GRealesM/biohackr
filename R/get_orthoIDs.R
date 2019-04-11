#' Get ortholog gene IDs from Ensembl
#'
#' Get Ensembl IDs of ortholog genes in desired species, that can subsequently used to retrieve their sequences using \code{get.seq} function.
#' @param genes A character vector of input genes. By default, Human gene names.
#' @param query_species A character indicating the query species name from which we will search orthologs. Defaults to Homo sapiens.
#' @param input A character indicating the type of gene identifiers. By default, gene names (external_gene_names). See \code{biomaRt::listFilters()} for more details.
#' @param set A character vector indicating names of species to search for orthologs. By default, all available species. However, I recommend to limit this to a number of species of interest.
#' @param host A character specifying the Ensembl mirror host. Defaults to the European mirror.
#' @param check A logical. If TRUE, it returns the Ensembl dictionary and stops. Useful for correcting mispelled species names
#' @param online A logical indicating if the dictionary should be created from online resources. If TRUE, it will run \code{get.dict}. Useful to get the most up-to-date dictionary, albeit slower. If FALSE, it will import dict dataset. Updated on April 2019. Defults to FALSE.
#'
#' @return A data frame containing Ensembl IDs and species names, among other information. If \code{check = TRUE}, it returns the Ensembl dictionary
#' @export
#'
#' @examples
#'
#' IDs <- get.orthoIDs(genes = c("APOE", "HTR3A", "PAX9"), set = c("Mus musculus", "Felis catus"))
#'
#'
get.orthoIDs <- function(genes = NULL, query_species = "Homo sapiens", input = "external_gene_name", set = "all", host = "www.ensembl.org",check = FALSE, online = FALSE){

  # Create or import the dictionary
  if(!online){
    dict <- dict
  } else {
      # Create dictionary from European server and use other mirrors if site is down
      cat("Creating dictionary...\n")
      dict <- try(get.dict(),silent = T)
      if("try-error" %in% class(dict)) {
        host = "useast.ensembl.org"
        dict <- try(get.dict("https://useast.ensembl.org/info/genome/stable_ids/prefixes.html", host = host), silent = T)

      }
      if("try-error" %in% class(dict)) {
        host = "uswest.ensembl.org"
        dict <- try(get.dict("https://uswest.ensembl.org/info/genome/stable_ids/prefixes.html", host = host), silent = T)
      }
      if("try-error" %in% class(dict)) {
        host = "asia.ensembl.org"
        dict <- try(get.dict("https://asia.ensembl.org/info/genome/stable_ids/prefixes.html", host = host), silent = T)
      }
      if("try-error" %in% class(dict)) stop("Dictionary creation failed. Server unreachable.")
      cat("Host is:", host, "\n")
  }
  # Sanity checks
  if(check) return(dict)
  if(is.null(genes)) stop("Please provide a list of gene names/symbols (e.g. APOE, HTR3A)\n")

  #Load mart
  if(is.null(query_species) || !query_species %in% dict$Species.Short.name) stop("Please, provide a valid query species name. For a list of valid species names, use flag check = TRUE, and check Species.Short.name column.\n")
  query_dataset <- dict[grep(query_species, dict$Species.Short.name),]$dataset
  CurlHandle <- RCurl::getCurlHandle()
  cat("Loading query mart: ",query_species,"...", sep = "")
  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = query_dataset, host = host, verbose = F, ensemblRedirect = F)
  cat("OK!\n")
  #Create Query list
  cat("Working on query dataset: ", query_species, "\n", sep="")
  bmquery <- biomaRt::getBM(c("ensembl_gene_id","external_gene_name"), filters = input, values = genes, mart = mart, curl = CurlHandle)
  bmquery$Confidence <- NA
  bmquery$Species.fasta <- gsub(" ", "_",query_species)
  bmquery <- bmquery[,c(4,2,1,3)]
  bmquery <- dplyr::left_join(bmquery, dict, by = "Species.fasta") [, c(1:4, 6, 8, 10:11)]
  names(bmquery)[1:3] <- c("Species.name","Gene.name", "Ensembl.Gene.ID")

  # We subset our species depending on our needs
  if(length(set) == 1 && set == "all"){
    cat("You chose all available species.\n")
  } else if(class(set) == "character"){
    cat("You chose a custom set of species. Checking species list...")
    set <- c(query_species, set)
    if (all(set %in% dict$Species.Short.name)){
      dict <- dict[match(set, dict$Species.Short.name),]
      cat("OK!\n")
    } else {
      notfound <- which(!set %in% dict$Species)
      stop(paste("Invalid species name(s): ", paste(set[notfound], collapse = ", "), "\nPlease, check spelling.\nFor a list of valid species names, use flag check = TRUE\n"))
    }
  } else {
    stop("Please, use a valid input.\n")
  }
  # Remove query species from list, so it will only take different species in its search for orthologs
  dict_noquery <- dict[-match(query_species, dict$Species.Short.name),]
  #Retrieve data from Ensembl
  cat("Retrieving IDs for orthologs in the list of species provided.\nDepending of the number of genes queried, this could take a while. Time for some coffee...\n")
  for(i in seq_along(dict_noquery$Species.Short.name)){
    cat("Working on", dict_noquery$Species.Short.name[i], "dataset\n", sep = " ")
    abb <- strsplit(dict_noquery$dataset[i], split = "_")[[1]][1]
    Bm <- biomaRt::getBM(attributes = c("external_gene_name", paste(abb, "_homolog_ensembl_gene", sep = ""), paste(abb, "_homolog_orthology_confidence", sep = "")), filters = input, values = genes, mart = mart, curl = CurlHandle)
    assign(paste(dict_noquery$Species.fasta[i], "-idset", sep = ""), Bm)
  }
  id_list <- mget(ls(pattern = "-idset"))
  id_list <- Map(cbind, id_list, Species.fasta = vapply(strsplit(names(id_list),"-"), `[`, 1, FUN.VALUE=character(1)))
  id_list <- data.table::rbindlist(id_list, use.names = F)
  id_list <- suppressWarnings(dplyr::left_join(id_list, dict, by = "Species.fasta")) [, c(4,1:3, 6, 8, 10:11)]
  names(id_list)[1:4] <- c("Species.name", "Gene.name", "Ensembl.Gene.ID", "Confidence")
  id_list[id_list == ""] <- NA
  id_list <- rbind(bmquery, id_list)
  cat("Done retrieving IDs!\n")
  cat("Done!\n")
  return(id_list)
}
