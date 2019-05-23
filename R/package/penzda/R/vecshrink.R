#' Soft threshholding operator.
#'
#' Applies the shrinkage operator to vector v with threshold a.
#' @param v vector
#' @param a threshold
#' @return s soft threshholded vector.
#' @export

vecshrink <- function(v,a){
  s <- sign(v)*pmax(abs(v)- a, 0 )
  return(s)
}