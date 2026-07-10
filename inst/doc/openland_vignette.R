## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----include=FALSE------------------------------------------------------------
url <- "https://zenodo.org/record/3685230/files/SaoLourencoBasin.rda?download=1"
dataset_available <- OpenLand:::.openland_try_download_and_load_rda(url,
  object = "SaoLourencoBasin",
  timeout = 10,
  cache = FALSE,
  quiet = TRUE
)

## ----include=FALSE, eval=isFALSE(dataset_available)---------------------------
knitr::knit_exit(paste0(
  "## The complete rendering of this document depends on the availability of the dataset in zenodo. ",
  "Please try at another time."
))

