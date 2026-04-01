#' Negation of `%in%`
#'
#' Returns `TRUE` when `x` is not contained in `table`.
#'
#' @param x A vector of values to match.
#' @param table A vector or list to match against.
#'
#' @return A logical vector.
#' @rdname not-in
#' @name not-in
#' @export
#'
#' @examples
#' 1:3 %!in% c(2, 4)
`%!in%` <- function(x, table) {
  !(x %in% table)
}

#' Reset a temporary directory
#'
#' Deletes and recreates a temporary directory, then updates `TMPDIR` and the
#' `bigstatsr` temporary-directory option.
#'
#' @param tempdir_path Path to the temporary directory.
#'
#' @return Invisibly returns `NULL`.
#' @export

reset_temp_dir <- function(tempdir_path) {
  # If the directory exists, delete it
  if (dir.exists(tempdir_path)) {
    unlink(tempdir_path, recursive = TRUE)  # Delete directory and its contents
  }

  # Recreate the temporary directory
  dir.create(tempdir_path, recursive = TRUE)

  # Set the directory for temporary files in R and bigstatsr
  Sys.setenv(TMPDIR = tempdir_path)
  options(bigstatsr.temporary_directory = tempdir_path)
}
