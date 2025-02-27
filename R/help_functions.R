#' str_eval
#'
#' This function is an inner function to evaluate the models.
#'
#' This function is used within the free order procedure
#'
#' @param x Model to be evaluates
#'
#' @author Julian Urban

#'
#' @return model
#' @export
#'
#'
str_eval <- function(x) {return(eval(parse(text = x), envir = parent.frame()))}
















