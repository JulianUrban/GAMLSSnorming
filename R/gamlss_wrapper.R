#' gamlss_wrapper
#'
#' This function conducts the free order procedure for any PDF using either
#' power polynomials or penalized splines
#'
#' This function is used to select the model.
#'
#' @param model A formula that defines the model for continuous norming.
#' @param data Defines the data set.
#' @param family A character that defines the specific PDF for modelling.
#' @param smooth A character defining the smooth function.
#' @param criterion Either "BIC", "AIC" or "GAIC(3)".
#' @param max.mu An integer. Defines the maximum power for polynomials.
#' @param max.sigma An integer. Defines the maximum power for polynomials.
#' @param max.nu An integer. Defines the maximum power for polynomials.
#' @param max.tau An integer. Defines the maximum power for polynomials.
#'
#' @author Julian Urban
#'
#' @import gamlss
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @importFrom gamlss GAIC
#' @importFrom stats poly
#'
#' @return The gamlss model selected by the free order procedure.
#' @export
#'
#' @examples
#' # Apply the free order procedure for a BCPE and the BIC.
#' gamlss_wrapper(model = g ~ e_age_month,
#'                data = Data_norm,
#'                criterion = "BIC",
#'                family = "NO",
#'                smooth = "pb")
#'
gamlss_wrapper<- function(model, ..., data, family, smooth, criterion, max.mu = 20, max.sigma = 20, max.nu = 20, max.tau = 20) {
  if(smooth == "poly") {
    output <- gamlss_poly(model = model,
                          ... = ...,
                          data = data,
                          family = family,
                          criterion = criterion,
                          max.mu = max.mu,
                          max.sigma = max.sigma,
                          max.nu = max.nu,
                          max.tau = max.tau)
  } else if(smooth == "pb") {
    output <- gamlss_pb(model = model,
                        ... = ...,
                        data = data,
                        family = family,
                        criterion = criterion)
  } else {
    break("Specified smooth function not available!")
  }
}

















