#' Convert a row/col index to a condensed index (for use with dist).
#'
#' @param i The row index of the uncondensed square matrix s.t. i!=j.
#' @param j The col index of the uncondensed square matrix s.t. i!=j.
#' @param n The length/width of the (square) matrix (i.e. the number of 
#'   observations passed into dist).
#' @return A integer value giving the corresponding index into the condensed
#'   matrix (vector).
#' @seealso \code{\link{dist}}, \code{\link{row_col_from_condensed_index}}
#' @export
#' @examples
#' data(iris)
#' n = nrow(iris)
#' dx = dist(iris[,-5])
#' i = 14; j=119
#' 
#' ix = condensed_index_from_row_col(i,j,n) # 1964
#' dist(iris[c(i,j),-5]) == dx[ix] # TRUE

# from http://stackoverflow.com/a/12643509/819544
condensed_index_from_row_col = function(i,j,n){
  n*(i-1) - i*(i-1)/2 + j-i
}
