devtools::document()
knitr::knit("vignettes/ncbi_bold.Rmd.orig", output = "vignettes/ncbi_bold.Rmd")
pkgdown::build_site()
