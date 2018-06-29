#' @useDynLib tnorm
#' @importFrom Rcpp sourceCpp
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("tnorm", libpath)
}
