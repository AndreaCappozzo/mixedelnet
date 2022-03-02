#' @export
control_EM <- function(itermax = 1000,
                       tol = 1e-08,
                       err = .Machine$double.xmax / 2,
                         BETA_zero_tol=1e-7, CD_threshold=1e-8) {
  list(
    itermax = itermax,
    tol = tol,
    err = err,
    BETA_zero_tol = BETA_zero_tol,
    CD_threshold = CD_threshold # Convergence threshold for coordinate descent
  )
}
