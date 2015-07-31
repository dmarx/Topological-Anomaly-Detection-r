#' Convert a condensed index to a row/col index (for use with dist).
#'
#' @param n The length/width of the (square) matrix (i.e. the number of 
#'   observations passed into dist).
#' @param ix The condensed index to be converted
#' @return a row and column index into the uncondensed square matrix.
#' @seealso \code{\link{dist}}, \code{\link{condensed_index_from_row_col}}
#' @export
#' @examples
#' data(iris)
#' n = nrow(iris)
#' dx = dist(iris[,-5])
#' ix = which.max(dx) # 1964
#' 
#' rc_ix = row_col_from_condensed_index(n, ix) # 14, 119
#' dist(iris[rc_ix,-5]) == max(dx) # TRUE

# from http://stackoverflow.com/a/12643509/819544
row_col_from_condensed_index = function(n, ix){
  nr=ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
  nc=n-(2*n-nr+1)*nr/2+ix+nr
  cbind(nr,nc)
}
