#' Get sequences from Ensembl IDs
#'
#' This function retrieves the coding or protein sequence from a list of identifiers.
#' @param IDs A character vector containing the Ensembl identifiers, consistent with input_type
#' @param seqtype Type of output sequence, coding or protein ("coding", "peptide", respectively). See \code{biomaRt::listAttributes()} for more details.
#' @param input_type A character indicating the input type. By default Ensembl gene IDs (ensembl_gene_id). See \code{biomaRt::listFilters} for more details.
#' @param longest A logical indicating whether, in case a gene ID has several available transcripts, only the longest one should be kept. If FALSE, it will keep all available transcript sequences.
#' @param host A character indicating the Ensembl mirror host. Defaults to the European mirror.
#' @param online A logical indicating if the dictionary should be created from online resources. If TRUE, it will run \code{get.dict}. Useful to get the most up-to-date dictionary, albeit slower. If FALSE, it will import dict dataset. Updated on April 2019. Defults to FALSE.
#'
#' @return A data.frame with the desired sequences, corresponding IDs, and species names, among others
#' @export
#'
#' @examples
#' Seqfile <- get.seq(IDs = c("ENSG00000166736", "ENSMUSG00000032269"), seqtype = "peptide", host = "uswest.ensembl.org")
#' @importFrom magrittr %>%
#'
#'
get.seq <- function(IDs, seqtype = NULL, input_type = "ensembl_gene_id", longest = TRUE, host = "www.ensembl.org", online = FALSE){

  # Sanity checks
  if(is.null(seqtype)|| !seqtype %in% c("coding", "peptide")) stop("Please, choose a valid seqtype: coding, protein")
  if(is.null(IDs)) stop("Please, provide a vector of gene identifiers")

  # Get dictionary
  if(!online){
    dict <- dict
  } else{
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
  # Check if IDs are alright
  geneslist <- data.frame(x = IDs, stringsAsFactors = F)
  genestable <- geneslist %>% fuzzyjoin::regex_inner_join(dict, by = c(x = "Prefix"))
  if (length(setdiff(geneslist$x, genestable$x)) != 0)
    warning("The following IDs were not recognized: \n", paste(setdiff(geneslist$x, genestable$x), collapse = ", "))

  # Load mart and get curl handle, which improves function speed
  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = host)
  CurlHandle <- RCurl::getCurlHandle()
  list_of_datasets <- unique(genestable$dataset)
  cat("You provided", length(geneslist$x), "Gene IDs and",length(list_of_datasets), "species.\n" , sep = " ")
  cat("You selected", seqtype,"sequence type!\nGathering requested information, this could take a while...\n")
  for(i in 1:length(list_of_datasets)){
    DS <- biomaRt::useDataset(list_of_datasets[i], mart = mart)
    genes_for_BM <- genestable %>% dplyr::filter(dataset == list_of_datasets[i]) %>% .$x
    species_for_BM <- genestable %>% dplyr::filter(dataset == list_of_datasets[i]) %>% .$Species.Short.name %>% .[1]
    cat("Working on", species_for_BM, "dataset\n", sep = " ")
    BM <- biomaRt::getBM(c("ensembl_gene_id", "ensembl_transcript_id", seqtype), filters = input_type, values = genes_for_BM, mart = DS, curl = CurlHandle, quote = "")
    names(BM) <- gsub("coding|peptide", "seq",names(BM))
    BM$length <- nchar(BM$seq)
    BM$Species_name <- species_for_BM
    assign(paste(species_for_BM, "_seqs"), BM)
    rm(DS)
  }
  listBM <- mget(ls(pattern = "_seqs"))
  seq_df <- as.data.frame(data.table::rbindlist(listBM))[,c(5,3,1,4,2)]
  names(seq_df) <- c("Species.name", "Ensembl.Gene.ID", "Ensembl.Transcript.ID", "Length", "Sequence")

  # Sequence actions
  if(longest == T){
    cat("You chose to keep the longest sequence per gene. Processing...\n")
    seq_df <- seq_df %>% dplyr::group_by(Ensembl.Gene.ID) %>% dplyr::mutate(rank = rank(dplyr::desc(Length), ties.method = "random")) %>% dplyr::filter(rank == 1) %>% dplyr::select(-rank)
  }
  return(seq_df)
}
