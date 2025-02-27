#' para_norm_extraction
#'
#' This function extracts the norm scores of a gamlss norming model.
#'
#' This function provide norm scores.
#'
#' @param model GAMLSS model for which norm scores should be extracted
#' @param scale Defines the scale of the norm scores
#'
#' @author Julian Urban
#'
#' @import gamlss.dist
#'
#' @return The norm scores
#'
#'
para_norm_extraction <- function(model,
                                 scale = "percentiles") {
  family <- model[["family"]][1]
  if(family == "BCPE") {
  norm_scores <- gamlss.dist::pBCPE(q = model[["y"]],
                                    mu = model[["mu.fv"]],
                                    sigma = model[["sigma.fv"]],
                                    nu = model[["nu.fv"]],
                                    tau = model[["tau.fv"]])
  } else if(family == "NO") {
    norm_scores <- gamlss.dist::pNO(q = model[["y"]],
                                    mu = model[["mu.fv"]],
                                    sigma = model[["sigma.fv"]])
  } else if(family == "SHASH") {
    norm_scores <- gamlss.dist::pSHASH(q = model[["y"]],
                                       mu = model[["mu.fv"]],
                                       sigma = model[["sigma.fv"]],
                                       nu = model[["nu.fv"]],
                                       tau = model[["tau.fv"]])
  }
  if(scale == "z") {
    norm_scores <- qnorm(norm_scores)
  } else if(scale == "T") {
    norm_scores <- qnorm(norm_scores) * 10 + 50
  } else if(scale == "IQ") {
    norm_scores <- qnorm(norm_scores) * 15 + 100
  } else if(scale == "Z") {
    norm_scores <- qnorm(norm_scores) * 10 + 100
  } else if(scale == "PISA") {
    norm_scores <- qnorm(norm_scores) * 100 + 500
  }

  return(norm_scores)

}
