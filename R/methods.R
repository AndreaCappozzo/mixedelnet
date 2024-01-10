#' @export
predict.lmm <- function(object, new_data, ...){
  if (!inherits(object, "lmm"))
    stop("object not of class \"lmm\"")

  args <- list(...)

  if(length(args)==0){ # if the grouping variable is not known for newdata, I simply use the average estimates
    return(c(object$beta %*% t(new_data)))
  }

  grouping_variable <- unlist(new_data[, args$grouping_variable_name])
  model_matrix_new_data <-
    subset(new_data, select = -get(args$grouping_variable_name))
  random_effect_component <- object$mu_raneff[as.character(grouping_variable)]
  random_effect_component[is.na(random_effect_component)] <- 0 # if there are obs whose groups are not in the mu_raneff I use the expected value of the ranef that is 0
  c(object$beta %*% t(model_matrix_new_data)) +
    random_effect_component
}
