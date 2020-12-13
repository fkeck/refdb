

refdb_merge <-function(...) {

  x <- list(...)

  stopifnot("All reference db objects must be dataframes" =
              any(!sapply(x, is.data.frame))
            )
  stopifnot("xxxx" =
              sapply(x, attributes)
  )


}


