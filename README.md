
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biohackr

<!-- badges: start -->

<!-- badges: end -->

**Biohackr** is a package devoted to fix some
[biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
limitations, and it is [available at
GitHub](https://github.com/GRealesM/biohackr/). Ensembl makes genomic
data freely accessible through biomaRt R package, which allows to query
to Ensembl servers from R. However, Ensembl datasets are organized by
species, and therefore it is not straightforward to extract information
and sequences (i.e. sequences from ortholog genes) across different
species in just one call. Biohackr intends to simplify things by
providing a simple hack for retrieving homolog IDs and sequences from
many different species at once.

## Installation

**Important note**: To use biohackr, you need to have some packages
installed first. Most of them are installed automatically from
[CRAN](https://cran.r-project.org/) when installing biohackr. However,
biomaRt (arguably the most important one), being a
[Bioconductor](https://bioconductor.org) package, it must be installed a
little bit differently:

``` r
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install() # For installing core Bioconductor packages
BiocManager::install("biomaRt") # For installing biomaRt
```

Once you have biomaRt installed, you may install the released version of
biohackr from [GitHub](https://github.com/GRealesM/) with:

``` r
library(devtools)
devtools::install_github("GRealesM/biohackr")
```

## Example

We may have a list of genes of interest and may want to quickly download
the coding or protein sequences of several ortholog genes in other
species of interest for comparative purporses and further
(e.g. molecular evolution) analyses. First, we define a list of query
gene names. By default get.orthoIDs take human gene names as an input,
but that behavior can be changed by modifying query\_species and input
parameters. By default, we can use human gene names.

``` r
genes <- c("HTR3A", "HTR3B", "PAX9", "APOE")
```

We also define the list of our species of interest by their scientific
name. In case any of them is mispelled, biohackr will throw an error
indicating which names were misspelled. By default it retrieves
***all*** available species/strains in Ensembl (184 in release 96, April
2019), which may take long and eventually crash if connection with the
server is lost, so I recommend to use a custom set of species of
interest. In this case I chose a random sample of
species.

``` r
species_set <- c("Vulpes vulpes","Takifugu rubripes","Canis lupus familiaris","Rhinopithecus bieti","Anser brachyrhynchus","Coturnix japonica","Peromyscus maniculatus bairdii","Urocitellus parryii","Canis lupus dingo","Ochotona princeps","Chinchilla lanigera")
```

First, we’ll need to retrieve all Ensembl IDs (e.g. `ENSG00000166736`)
for querying the server for sequences afterwards. We’ll use
`get.orthoIDs()` for that first
step

``` r
orthoIDs <- get.orthoIDs(genes, host = "uswest.ensembl.org", set = species_set)
```

The default host is www.ensembl.org (UK mirror). However, sometimes
different mirrors may be down. The function will try to find a
functioning mirror if default is down. However, we can override this
behavior by specifying a different mirror ([US
West](https://uswest.ensembl.org) in this case). Other mirrors are
[useast.ensembl.org](https://useast.ensembl.org) (US East) and
[asia.ensembl.org](https://asia.ensembl.org) (Asia).

Once we have our list of identifiers, we remove duplicated identifiers,
and instances in which no ortholog was found, then we use `get.seq` for
retrieving the seqs.

``` r
IDs <- unique(na.omit(orthoIDs$Ensembl.Gene.ID))
orthocds <- get.seq(IDs, seqtype = "coding", host = "uswest.ensembl.org")
```

By default, `get.seq` will keep only the longest transcript for each
gene ID, but this behavior can be changed by setting `longest = FALSE`.
This will keep all transcripts. Now we can merge both datasets and
export the result as fasta using the `seqinr` package. The following
code will result in a fasta file named `biohackR_example.fas`, with
headers in the format `>Species_name|TranscriptID|Original_gene_name`.

``` r
library(seqinr)
orthomerged <- dplyr::left_join(orthoIDs, orthocds, "Ensembl.Gene.ID")
orthomerged$Fasta.name <- paste(orthomerged$Species.name.x, orthomerged$Ensembl.Transcript.ID, orthomerged$Gene.name, sep = "|")
orthomerged <- na.omit(orthomerged)
seqinr::write.fasta(sequences = as.list(orthomerged$Sequence), names = orthomerged$Fasta.name, file.out = "biohackR_example.fas")
```
