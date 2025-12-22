# msdialSE

<!-- badges: start -->
<!-- badges: end -->

`msdialSE` is an R package that helps import and handle MS-DIAL outputs as
`SummarizedExperiment` objects, making it easier to integrate downstream
analysis and visualization workflows in R.

## Installation

You can install the development version from GitHub:

``` r
install.packages("remotes")
remotes::install_github("takakioka/msdialSE")

#library(msdialSE)
# se <- load_lipidomics_se("path/to/msdial_alignment.csv")
# se
