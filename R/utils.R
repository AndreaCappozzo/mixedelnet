#' @export
control_EM <- function(itermax = 1000,
                       tol = 1e-08,
                       err = .Machine$double.xmax / 2,CD_threshold=1e-8) {
  list(
    itermax = itermax,
    tol = tol,
    err = err,
    CD_threshold = CD_threshold # Convergence threshold for coordinate descent
  )
}
