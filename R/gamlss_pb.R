#' gamlss_pb
#'
#' This function is in inner function and conducts the free order procedure
#' with penalized splines.
#'
#' This function is used for penalized splines
#'
#' @param model A formula that defines the model for continuous norming.
#' @param data Defines the data set.
#' @param family A character that defines the specific PDF for modelling.
#' @param criterion Either "BIC", "AIC" or "GAIC(3)".
#'
#' @author Julian Urban
#'
#' @import gamlss
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @importFrom stats GAIC
#' @importFrom stats poly
#'
#' @return The gamlss model selected by the free order procedure.
#'
#' @examples
#' # Apply the free order procedure for a BCPE and the BIC.
#' gamlss_pb(model = g ~ Normgruppe,
#'           data = Data_norm,
#'           criterion = "BIC",
#'           family = "NO")
#'
gamlss_pb <- function(model, ..., data, family, criterion) {
  if(family %in% c("BCPE", "BCT", "EGB2",
                   "GB2", "GT", "JSU",
                   "NET", "SEP1", "SEP2",
                   "SEP3", "SEP4", "SHASH",
                   "SHASHo")) {
    output <- four_parameters_pb_splines(model = model,
                                         ... = ...,
                                         data = data,
                                         family = family,
                                         criterion = criterion)
  } else if(family %in% c("BNB", "BCCG", "DBURR12",
                          "exGAUS", "GG", "LNO",
                          "NOF", "PE", "PE2",
                          "ST1", "ST2", "ST3",
                          "ST4", "ST5", "TF")) {
    output <- three_parameters_pb_splines(model = model,
                                          ... = ...,
                                          data = data,
                                          family = family,
                                          criterion = criterion)
  } else if(family %in% c("DPO", "GA", "GU",
                          "LO", "LOGNO", "NO",
                          "LQNO", "PARETO2", "PARETO2o",
                          "PIG", "RG", "WEI",
                          "WEI3", "ZALG", "ZAZIPF")) {
    output <- two_parameters_pb_splines(model = model,
                                        ... = ...,
                                        data = data,
                                        family = family,
                                        criterion = criterion)
  } else {
    break("Specified family not available!")
  }
}

#' four_parameters_pb_splines
#'
#' This function is in inner function and conducts the free order procedure
#' with penalized splines for PDFs with four parameters.
#'
#' This function is used for penalized splines
#'
#' @param model A formula that defines the model for continuous norming.
#' @param data Defines the data set.
#' @param family A character that defines the specific PDF for modelling.
#' @param criterion Either "BIC", "AIC" or "GAIC(3)".
#'
#' @author Julian Urban
#'
#' @import gamlss
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @importFrom stats GAIC
#' @importFrom stats poly
#'
#' @return The gamlss model selected by the free order procedure.
#' @export
#'
#' @examples
#' # Apply the free order procedure for a BCPE and the BIC.
#' four_parameters_pb_splines(model = g ~ Normgruppe,
#'                            data = Data_norm,
#'                            criterion = "BIC",
#'                            family = "BCPE")
#'
four_parameters_pb_splines <- function(model, ..., data, criterion, family) {

  #===============================================================================
  # STEP 1: Preperation
  #-------------------------------------------------------------------------------
  model_string <- deparse(model)
  variables <- gsub(" ", "", strsplit(model_string, "~")[[1]])
  data_use <- data[rowSums(is.na(data[, variables])) == 0, variables]

  n <- nrow(data_use)
  if (criterion == "AIC"){
    k <- 2
  } else if (criterion == "BIC"){
    k <- log(n)
  } else if (criterion == "GAIC(3)") {
    k <- 3
  }

  if(min(data_use[, variables[1]]) < 0) {
    data_use[, variables[1]] <- data_use[, variables[1]] - min(data_use[, variables[1]]) + 0.000001
  }

  best_model <- "initial"

  mu <- 1
  sigma <- 1
  nu <- 1
  tau <- 1
  i.mu <- mu
  i.sigma <- sigma
  i.nu <- nu
  i.tau <- tau

  param <- c("1", paste("gamlss::pb(", variables[2],",method = 'GAIC', k =", k, ")"))

  #===============================================================================
  # STEP 2: Fit initial model & prepare free order procedure
  #-------------------------------------------------------------------------------
  init.model <- {paste("","gamlss::gamlss(", model_string, ",sigma.formula = ~ 1,
                     nu.formula = ~ 1, tau.formula = ~ 1, family = family",
                       ", data = data_use,",""," method = RS(10000))",
                       sep = "")}
  init.model.output <- str_eval(noquote(init.model))

  if (criterion == "AIC"){
    init.model.output.crit <- AIC(init.model.output)
  } else if (criterion == "BIC"){
    init.model.output.crit <- BIC(init.model.output)
  } else if (criterion == "GAIC(3)"){
    init.model.output.crit <- GAIC(init.model.output, k = 3)
  }

  parameter.output <- matrix(NA, nrow = 100, ncol = 6);
  colnames(parameter.output) <- c("nr", "mu", "sigma", "nu", "tau", "criterion")

  po <- 1
  parameter.output[po,] <- c(po, mu, sigma, nu, tau, init.model.output.crit);

  endloop <- "no"

  #===============================================================================
  # STEP 3: Main part of free order procedure
  #-------------------------------------------------------------------------------
  while (endloop != "yes") {

    output <- rep(NA, times = 9)

    #===============================================================================
    # STEP 3.1: Forward model selection
    #-------------------------------------------------------------------------------
    if(mu == 1) {
      i.mu <- 2
    }
    if(sigma == 1) {
      i.sigma <- 2
    }
    if(nu == 1) {
      i.nu <- 2
    }
    if(tau == 1) {
      i.tau <- 2
    }
    # Define the 4 forward models
    mu.model.forward <- {paste('',
                               'gamlss::gamlss(', variables[1],
                               '~', param[i.mu] ,
                               ', sigma.formula = ~ ', param[sigma],
                               ', nu.formula = ~ ', param[nu],
                               ', tau.formula = ~ ', param[tau],
                               ', family = family ',
                               ', data = data_use,',
                               "",
                               'method = RS(10000))',
                               sep = '')}

    sigma.model.forward <- {paste('',
                                  'gamlss::gamlss(', variables[1],
                                  '~', param[mu] ,
                                  ', sigma.formula = ~ ', param[i.sigma],
                                  ', nu.formula = ~ ', param[nu],
                                  ', tau.formula = ~ ', param[tau],
                                  ', family = family ',
                                  ', data = data_use,',
                                  "",
                                  'method = RS(10000))',
                                  sep = '')}

    nu.model.forward <- {paste('',
                               'gamlss::gamlss(', variables[1],
                               '~', param[mu] ,
                               ', sigma.formula = ~ ', param[sigma],
                               ', nu.formula = ~ ', param[i.nu],
                               ', tau.formula = ~ ', param[tau],
                               ', family = family ',
                               ', data = data_use,',
                               "",
                               'method = RS(10000))',
                               sep = '')}

    tau.model.forward <- {paste('',
                                'gamlss::gamlss(', variables[1],
                                '~', param[mu] ,
                                ', sigma.formula = ~ ', param[sigma],
                                ', nu.formula = ~ ', param[nu],
                                ', tau.formula = ~ ', param[i.tau],
                                ', family = family ',
                                ', data = data_use,',
                                "",
                                'method = RS(10000))',
                                sep = '')}

    mu.model.forward.output <- str_eval(noquote(mu.model.forward))
    sigma.model.forward.output <- str_eval(noquote(sigma.model.forward));
    nu.model.forward.output <- str_eval(noquote(nu.model.forward));
    tau.model.forward.output <- str_eval(noquote(tau.model.forward))

    if(criterion == "AIC"){
      mu.model.forward.output.crit <- AIC(mu.model.forward.output);
      sigma.model.forward.output.crit <- AIC(sigma.model.forward.output);
      nu.model.forward.output.crit <- AIC(nu.model.forward.output);
      tau.model.forward.output.crit <- AIC(tau.model.forward.output)
    } else if(criterion == "BIC"){
      mu.model.forward.output.crit <- BIC(mu.model.forward.output);
      sigma.model.forward.output.crit <- BIC(sigma.model.forward.output);
      nu.model.forward.output.crit <- BIC(nu.model.forward.output);
      tau.model.forward.output.crit <- BIC(tau.model.forward.output)
    } else if(criterion == "GAIC(3)") {
      mu.model.forward.output.crit <- GAIC(mu.model.forward.output, k = 3);
      sigma.model.forward.output.crit <- GAIC(sigma.model.forward.output, k = 3);
      nu.model.forward.output.crit <- GAIC(nu.model.forward.output, k = 3);
      tau.model.forward.output.crit <- GAIC(tau.model.forward.output, k = 3)
    }

    #===============================================================================
    # STEP 3.2: Backward model selection
    #-------------------------------------------------------------------------------
    if(mu == 2) {
      i.mu <- 1
    }
    if(sigma == 2) {
      i.sigma <- 1
    }
    if(nu == 2) {
      i.nu <- 1
    }
    if(tau == 2) {
      i.tau <- 1
    }

    if(mu == 2) {
      i.mu <- 1
      mu.model.backward <- {paste('',
                                  'gamlss::gamlss(', variables[1],
                                  '~', param[i.mu] ,
                                  ', sigma.formula = ~ ', param[sigma],
                                  ', nu.formula = ~ ', param[nu],
                                  ', tau.formula = ~ ', param[tau],
                                  ', family = family ',
                                  ', data = data_use,',
                                  "",
                                  'method = RS(10000))',
                                  sep = '')}

      mu.model.backward.output <- str_eval(noquote(mu.model.backward))

      if (criterion == "AIC"){
        mu.model.backward.output.crit <- AIC(mu.model.backward.output)
      } else if (criterion == "BIC"){
        mu.model.backward.output.crit <- BIC(mu.model.backward.output)
      } else if (criterion == "GAIC(3)"){
        mu.model.backward.output.crit <- GAIC(mu.model.backward.output, k = 3)
      }
    } else {
      mu.model.backward.output.crit <- Inf
    }

    if(sigma == 2) {
      i.sigma <- 1
      sigma.model.backward <- {paste('',
                                     'gamlss::gamlss(', variables[1],
                                     '~', param[mu] ,
                                     ', sigma.formula = ~ ', param[i.sigma],
                                     ', nu.formula = ~ ', param[nu],
                                     ', tau.formula = ~ ', param[tau],
                                     ', family = family ',
                                     ', data = data_use,',
                                     "",
                                     'method = RS(10000))',
                                     sep = '')}

      sigma.model.backward.output <- str_eval(noquote(sigma.model.backward))

      if (criterion == "AIC"){
        sigma.model.backward.output.crit <- AIC(sigma.model.backward.output)
      } else if (criterion == "BIC"){
        sigma.model.backward.output.crit <- BIC(sigma.model.backward.output)
      } else if (criterion == "GAIC(3)"){
        sigma.model.backward.output.crit <- GAIC(sigma.model.backward.output, k = 3)
      }
    } else {
      sigma.model.backward.output.crit <- Inf
    }

    if(nu == 2) {
      i.nu <- 1
      nu.model.backward <- {paste('',
                                  'gamlss::gamlss(', variables[1],
                                  '~', param[mu] ,
                                  ', sigma.formula = ~ ', param[sigma],
                                  ', nu.formula = ~ ', param[i.nu],
                                  ', tau.formula = ~ ', param[tau],
                                  ', family = family ',
                                  ', data = data_use,',
                                  "",
                                  'method = RS(10000))',
                                  sep = '')}

      nu.model.backward.output <- str_eval(noquote(nu.model.backward))

      if (criterion == "AIC"){
        nu.model.backward.output.crit <- AIC(nu.model.backward.output)
      } else if (criterion == "BIC"){
        nu.model.backward.output.crit <- BIC(nu.model.backward.output)
      } else if (criterion == "GAIC(3)"){
        nu.model.backward.output.crit <- GAIC(nu.model.backward.output, k = 3)
      }
    } else {
      nu.model.backward.output.crit <- Inf
    }

    if(tau == 2) {
      i.tau <- 1
      tau.model.backward <- {paste('',
                                   'gamlss::gamlss(', variables[1],
                                   '~', param[mu] ,
                                   ', sigma.formula = ~ ', param[sigma],
                                   ', nu.formula = ~ ', param[nu],
                                   ', tau.formula = ~ ', param[i.tau],
                                   ', family = family ',
                                   ', data = data_use,',
                                   "",
                                   'method = RS(10000))',
                                   sep = '')}

      tau.model.backward.output <- str_eval(noquote(tau.model.backward))

      if (criterion == "AIC"){
        tau.model.backward.output.crit <- AIC(tau.model.backward.output)
      } else if (criterion == "BIC"){
        tau.model.backward.output.crit <- BIC(tau.model.backward.output)
      } else if (criterion == "GAIC(3)"){
        tau.model.backward.output.crit <- GAIC(tau.model.backward.output, k = 3)
      }
    } else {
      tau.model.backward.output.crit <- Inf
    }


    #===============================================================================
    # STEP 4: Output creation & Evaluation
    #-------------------------------------------------------------------------------
    output <- {c(init.model.output.crit, mu.model.forward.output.crit,
                 sigma.model.forward.output.crit, nu.model.forward.output.crit,
                 tau.model.forward.output.crit, mu.model.backward.output.crit,
                 sigma.model.backward.output.crit, nu.model.backward.output.crit,
                 tau.model.backward.output.crit)}
    names(output) <- {c("initial", "mu.forward", "sigma.forward", "nu.forward",
                        "tau.forward", "mu.backward", "sigma.backward", "nu.backward", "tau.backward")}


    if ({output["mu.forward"] < output["initial"] & output["mu.forward"] <
        output["sigma.forward"] & output["mu.forward"] < output["nu.forward"] &
        output["mu.forward"] < output["tau.forward"] & output["mu.forward"] <
        output["mu.backward"] & output["mu.forward"] < output["sigma.backward"] &
        output["mu.forward"] < output["nu.backward"] & output["mu.forward"] <
        output["tau.backward"]})
    {mu <- 2
    init.model.output.crit <- mu.model.forward.output.crit
    best_model <- "mu.forward"
    } else {
      if ({output["sigma.forward"] < output["initial"] &
          output["sigma.forward"] < output["mu.forward"] &
          output["sigma.forward"] < output["nu.forward"] &
          output["sigma.forward"] < output["tau.forward"] &
          output["sigma.forward"] < output["mu.backward"] &
          output["sigma.forward"] < output["sigma.backward"] &
          output["sigma.forward"] < output["nu.backward"] &
          output["sigma.forward"] < output["tau.backward"]}) {
        sigma <- 2
        init.model.output.crit <- sigma.model.forward.output.crit
        best_model <- "sigma.forward"
      } else {
        if ({output["nu.forward"] < output["initial"] & output["nu.forward"] <
            output["mu.forward"] & output["nu.forward"] < output["sigma.forward"] &
            output["nu.forward"] < output["tau.forward"] & output["nu.forward"] <
            output["mu.backward"] & output["nu.forward"] < output["sigma.backward"] &
            output["nu.forward"] < output["nu.backward"] & output["nu.forward"] <
            output["tau.backward"]}) {
          nu <- 2
          init.model.output.crit <- nu.model.forward.output.crit
          best_model <- "nu.forward"
        } else {
          if ({output["tau.forward"] < output["initial"] &
              output["tau.forward"] < output["mu.forward"] &
              output["tau.forward"] < output["sigma.forward"] &
              output["tau.forward"] < output["nu.forward"] &
              output["tau.forward"] < output["mu.backward"] &
              output["tau.forward"] < output["sigma.backward"] &
              output["tau.forward"] < output["nu.backward"] &
              output["tau.forward"] < output["tau.backward"]}) {
            tau <- 2
            init.model.output.crit <- tau.model.forward.output.crit
            best_model <- "tau.forward"
          } else {
            if ({output["mu.backward"] < output["initial"] &
                output["mu.backward"] < output["mu.forward"] &
                output["mu.backward"] < output["sigma.forward"] &
                output["mu.backward"] < output["nu.forward"] &
                output["mu.backward"] < output["tau.forward"] &
                output["mu.backward"] < output["sigma.backward"] &
                output["mu.backward"] < output["nu.backward"] &
                output["mu.backward"] < output["tau.backward"]}) {
              mu <- 1
              init.model.output.crit <- mu.model.backward.output.crit
              best_model <- "mu.backward"
            } else {
              if ({output["sigma.backward"] < output["initial"] &
                  output["sigma.backward"] < output["mu.forward"] &
                  output["sigma.backward"] < output["sigma.forward"] &
                  output["sigma.backward"] < output["nu.forward"] &
                  output["sigma.backward"] < output["tau.forward"] &
                  output["sigma.backward"] < output["mu.backward"] &
                  output["sigma.backward"] < output["nu.backward"] &
                  output["sigma.backward"] < output["tau.backward"]}) {
                sigma <- 1
                init.model.output.crit <- sigma.model.backward.output.crit
                best_model <- "sigma.backward"
              } else {
                if ({output["nu.backward"] < output["initial"] &
                    output["nu.backward"] < output["mu.forward"] &
                    output["nu.backward"] < output["sigma.forward"] &
                    output["nu.backward"] < output["nu.forward"] &
                    output["nu.backward"] < output["tau.forward"] &
                    output["nu.backward"] < output["mu.backward"] &
                    output["nu.backward"] < output["sigma.backward"] &
                    output["nu.backward"] < output["tau.backward"]}) {
                  nu <- 1
                  init.model.output.crit <- nu.model.backward.output.crit
                  best_model <- "nu.backward"
                } else {
                  if ({output["tau.backward"] < output["initial"] &
                      output["tau.backward"] < output["mu.forward"] &
                      output["tau.backward"] < output["sigma.forward"] &
                      output["tau.backward"] < output["nu.forward"] &
                      output["tau.backward"] < output["tau.forward"] &
                      output["tau.backward"] < output["mu.backward"] &
                      output["tau.backward"] < output["sigma.backward"] &
                      output["tau.backward"] < output["nu.backward"]}) {
                    tau <- 1
                    init.model.output.crit <- tau.model.backward.output.crit
                    best_model <- "tau.backward"
                  } else {
                    endloop <- "yes" }}}}}}}}
  } # end while loop

  #===============================================================================
  # STEP 5: Build output
  #-------------------------------------------------------------------------------
  parameter.output <- matrix(NA, nrow = 100, ncol = 6);
  colnames(parameter.output) <- c("nr", "mu", "sigma", "nu", "tau", "criterion")

  po <- 1
  parameter.output[po,] <- c(po, mu, sigma, nu, tau, init.model.output.crit)

  po <- po + 1
  if (init.model.output.crit != "NA"){
    {parameter.output[po,] <- c(po, mu, sigma, nu, tau,
                                round(init.model.output.crit,3))}
  }
  if (init.model.output.crit == "NA"){
    {parameter.output[po,] <- c(po, mu, sigma, nu, tau, init.model.output.crit)}
  }
  found.parameters.full <- parameter.output[, 2:5]
  ind <- !is.na(found.parameters.full)
  found.parameters <- {tapply(found.parameters.full[ind],
                              col(found.parameters.full)[ind], tail, 1)}
  names(found.parameters) <- colnames(found.parameters.full)

  if(best_model == "initial") {
    model_result <- init.model.output
  } else {
    best_model_split <- strsplit(best_model, "\\.")[[1]]
    model_result <- eval(parse(text = paste(best_model_split[1], "model", best_model_split[2], "output", sep = ".", collapse = "")))
  }

  return(model_result)
} # end function

#' three_parameters_pb_splines
#'
#' This function is in inner function and conducts the free order procedure
#' with penalized splines for PDFs with three parameters.
#'
#' This function is used for penalized splines
#'
#' @param model A formula that defines the model for continuous norming.
#' @param data Defines the data set.
#' @param family A character that defines the specific PDF for modelling.
#' @param criterion Either "BIC", "AIC" or "GAIC(3)".
#'
#' @author Julian Urban
#'
#' @import gamlss
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @importFrom stats GAIC
#' @importFrom stats poly
#'
#' @return The gamlss model selected by the free order procedure.
#' @export
#'
#' @examples
#' # Apply the free order procedure for a BCPE and the BIC.
#' three_parameters_pb_splines(model = g ~ Normgruppe,
#'                            data = Data_norm,
#'                            criterion = "BIC",
#'                            family = "BCPE")
three_parameters_pb_splines <- function(model, ..., data, criterion, family) {

  #===============================================================================
  # STEP 1: Preperation
  #-------------------------------------------------------------------------------
  model_string <- deparse(model)
  variables <- gsub(" ", "", strsplit(model_string, "~")[[1]])
  data_use <- data[rowSums(is.na(data[, variables])) == 0, variables]

  n <- nrow(data_use)
  if (criterion == "AIC"){
    k <- 2
  } else if (criterion == "BIC"){
    k <- log(n)
  } else if (criterion == "GAIC(3)") {
    k <- 3
  }

  if(min(data_use[, variables[1]]) < 0) {
    data_use[, variables[1]] <- data_use[, variables[1]] - min(data_use[, variables[1]]) + 0.000001
  }

  best_model <- "initial"

  mu <- 1
  sigma <- 1
  nu <- 1
  i.mu <- mu
  i.sigma <- sigma
  i.nu <- nu

  param <- c("1", paste("gamlss::pb(", variables[2],",method = 'GAIC', k =", k, ")"))

  #===============================================================================
  # STEP 2: Fit initial model & prepare free order procedure
  #-------------------------------------------------------------------------------
  init.model <- {paste("","gamlss::gamlss(", model_string, ",sigma.formula = ~ 1,
                     nu.formula = ~ 1, family = family",
                       ", data = data_use,",""," method = RS(10000))",
                       sep = "")}
  init.model.output <- str_eval(noquote(init.model))

  if (criterion == "AIC"){
    init.model.output.crit <- AIC(init.model.output)
  } else if (criterion == "BIC"){
    init.model.output.crit <- BIC(init.model.output)
  } else if (criterion == "GAIC(3)"){
    init.model.output.crit <- GAIC(init.model.output, k = 3)
  }

  parameter.output <- matrix(NA, nrow = 100, ncol = 5);
  colnames(parameter.output) <- c("nr", "mu", "sigma", "nu", "criterion")

  po <- 1
  parameter.output[po,] <- c(po, mu, sigma, nu, init.model.output.crit);

  endloop <- "no"

  #===============================================================================
  # STEP 3: Main part of free order procedure
  #-------------------------------------------------------------------------------
  while (endloop != "yes") {

    output <- rep(NA, times = 9)

    #===============================================================================
    # STEP 3.1: Forward model selection
    #-------------------------------------------------------------------------------
    if(mu == 1) {
      i.mu <- 2
    }
    if(sigma == 1) {
      i.sigma <- 2
    }
    if(nu == 1) {
      i.nu <- 2
    }

    # Define the 4 forward models
    mu.model.forward <- {paste('',
                               'gamlss::gamlss(', variables[1],
                               '~', param[i.mu] ,
                               ', sigma.formula = ~ ', param[sigma],
                               ', nu.formula = ~ ', param[nu],
                               ', family = family ',
                               ', data = data_use,',
                               "",
                               'method = RS(10000))',
                               sep = '')}

    sigma.model.forward <- {paste('',
                                  'gamlss::gamlss(', variables[1],
                                  '~', param[mu] ,
                                  ', sigma.formula = ~ ', param[i.sigma],
                                  ', nu.formula = ~ ', param[nu],
                                  ', family = family ',
                                  ', data = data_use,',
                                  "",
                                  'method = RS(10000))',
                                  sep = '')}

    nu.model.forward <- {paste('',
                               'gamlss::gamlss(', variables[1],
                               '~', param[mu] ,
                               ', sigma.formula = ~ ', param[sigma],
                               ', nu.formula = ~ ', param[i.nu],
                               ', family = family ',
                               ', data = data_use,',
                               "",
                               'method = RS(10000))',
                               sep = '')}

    mu.model.forward.output <- str_eval(noquote(mu.model.forward))
    sigma.model.forward.output <- str_eval(noquote(sigma.model.forward))
    nu.model.forward.output <- str_eval(noquote(nu.model.forward))

    if(criterion == "AIC"){
      mu.model.forward.output.crit <- AIC(mu.model.forward.output)
      sigma.model.forward.output.crit <- AIC(sigma.model.forward.output)
      nu.model.forward.output.crit <- AIC(nu.model.forward.output)
    } else if(criterion == "BIC"){
      mu.model.forward.output.crit <- BIC(mu.model.forward.output)
      sigma.model.forward.output.crit <- BIC(sigma.model.forward.output)
      nu.model.forward.output.crit <- BIC(nu.model.forward.output)
    } else if(criterion == "GAIC(3)") {
      mu.model.forward.output.crit <- GAIC(mu.model.forward.output, k = 3)
      sigma.model.forward.output.crit <- GAIC(sigma.model.forward.output, k = 3)
      nu.model.forward.output.crit <- GAIC(nu.model.forward.output, k = 3)
    }

    #===============================================================================
    # STEP 3.2: Backward model selection
    #-------------------------------------------------------------------------------
    if(mu == 2) {
      i.mu <- 1
    }
    if(sigma == 2) {
      i.sigma <- 1
    }
    if(nu == 2) {
      i.nu <- 1
    }

    if(mu == 2) {
      i.mu <- 1
      mu.model.backward <- {paste('',
                                  'gamlss::gamlss(', variables[1],
                                  '~', param[i.mu] ,
                                  ', sigma.formula = ~ ', param[sigma],
                                  ', nu.formula = ~ ', param[nu],
                                  ', family = family ',
                                  ', data = data_use,',
                                  "",
                                  'method = RS(10000))',
                                  sep = '')}

      mu.model.backward.output <- str_eval(noquote(mu.model.backward))

      if (criterion == "AIC"){
        mu.model.backward.output.crit <- AIC(mu.model.backward.output)
      } else if (criterion == "BIC"){
        mu.model.backward.output.crit <- BIC(mu.model.backward.output)
      } else if (criterion == "GAIC(3)"){
        mu.model.backward.output.crit <- GAIC(mu.model.backward.output, k = 3)
      }
    } else {
      mu.model.backward.output.crit <- Inf
    }

    if(sigma == 2) {
      i.sigma <- 1
      sigma.model.backward <- {paste('',
                                     'gamlss::gamlss(', variables[1],
                                     '~', param[mu] ,
                                     ', sigma.formula = ~ ', param[i.sigma],
                                     ', nu.formula = ~ ', param[nu],
                                     ', family = family ',
                                     ', data = data_use,',
                                     "",
                                     'method = RS(10000))',
                                     sep = '')}

      sigma.model.backward.output <- str_eval(noquote(sigma.model.backward))

      if (criterion == "AIC"){
        sigma.model.backward.output.crit <- AIC(sigma.model.backward.output)
      } else if (criterion == "BIC"){
        sigma.model.backward.output.crit <- BIC(sigma.model.backward.output)
      } else if (criterion == "GAIC(3)"){
        sigma.model.backward.output.crit <- GAIC(sigma.model.backward.output, k = 3)
      }
    } else {
      sigma.model.backward.output.crit <- Inf
    }

    if(nu == 2) {
      i.nu <- 1
      nu.model.backward <- {paste('',
                                  'gamlss::gamlss(', variables[1],
                                  '~', param[mu] ,
                                  ', sigma.formula = ~ ', param[sigma],
                                  ', nu.formula = ~ ', param[i.nu],
                                  ', family = family ',
                                  ', data = data_use,',
                                  "",
                                  'method = RS(10000))',
                                  sep = '')}

      nu.model.backward.output <- str_eval(noquote(nu.model.backward))

      if (criterion == "AIC"){
        nu.model.backward.output.crit <- AIC(nu.model.backward.output)
      } else if (criterion == "BIC"){
        nu.model.backward.output.crit <- BIC(nu.model.backward.output)
      } else if (criterion == "GAIC(3)"){
        nu.model.backward.output.crit <- GAIC(nu.model.backward.output, k = 3)
      }
    } else {
      nu.model.backward.output.crit <- Inf
    }

    #===============================================================================
    # STEP 4: Output creation & Evaluation
    #-------------------------------------------------------------------------------
    output <- {c(init.model.output.crit, mu.model.forward.output.crit,
                 sigma.model.forward.output.crit, nu.model.forward.output.crit,
                 mu.model.backward.output.crit,
                 sigma.model.backward.output.crit, nu.model.backward.output.crit)}
    names(output) <- {c("initial", "mu.forward", "sigma.forward", "nu.forward",
                        "mu.backward", "sigma.backward", "nu.backward")}


    if ({output["mu.forward"] < output["initial"] &
        output["mu.forward"] < output["sigma.forward"] &
        output["mu.forward"] < output["nu.forward"] &
        output["mu.forward"] <  output["mu.backward"] &
        output["mu.forward"] < output["sigma.backward"] &
        output["mu.forward"] < output["nu.backward"]})
    {mu <- 2
    init.model.output.crit <- mu.model.forward.output.crit
    best_model <- "mu.forward"
    } else {
      if ({output["sigma.forward"] < output["initial"] &
          output["sigma.forward"] < output["mu.forward"] &
          output["sigma.forward"] < output["nu.forward"] &
          output["sigma.forward"] < output["mu.backward"] &
          output["sigma.forward"] < output["sigma.backward"] &
          output["sigma.forward"] < output["nu.backward"]}) {
        sigma <- 2
        init.model.output.crit <- sigma.model.forward.output.crit
        best_model <- "sigma.forward"
      } else {
        if ({output["nu.forward"] < output["initial"] &
            output["nu.forward"] < output["mu.forward"] &
            output["nu.forward"] < output["sigma.forward"] &
            output["nu.forward"] < output["mu.backward"] &
            output["nu.forward"] < output["sigma.backward"] &
            output["nu.forward"] < output["nu.backward"]}) {
          nu <- 2
          init.model.output.crit <- nu.model.forward.output.crit
          best_model <- "nu.forward"
        } else {
          if ({output["mu.backward"] < output["initial"] &
              output["mu.backward"] < output["mu.forward"] &
              output["mu.backward"] < output["sigma.forward"] &
              output["mu.backward"] < output["nu.forward"] &
              output["mu.backward"] < output["sigma.backward"] &
              output["mu.backward"] < output["nu.backward"]}) {
            mu <- 1
            init.model.output.crit <- mu.model.backward.output.crit
            best_model <- "mu.backward"
          } else {
            if ({output["sigma.backward"] < output["initial"] &
                output["sigma.backward"] < output["mu.forward"] &
                output["sigma.backward"] < output["sigma.forward"] &
                output["sigma.backward"] < output["nu.forward"] &
                output["sigma.backward"] < output["mu.backward"] &
                output["sigma.backward"] < output["nu.backward"]}) {
              sigma <- 1
              init.model.output.crit <- sigma.model.backward.output.crit
              best_model <- "sigma.backward"
            } else {
              if ({output["nu.backward"] < output["initial"] &
                  output["nu.backward"] < output["mu.forward"] &
                  output["nu.backward"] < output["sigma.forward"] &
                  output["nu.backward"] < output["nu.forward"] &
                  output["nu.backward"] < output["mu.backward"] &
                  output["nu.backward"] < output["sigma.backward"]}) {
                nu <- 1
                init.model.output.crit <- nu.model.backward.output.crit
                best_model <- "nu.backward"
              }  else {
                endloop <- "yes" }}}}}}
  } # end while loop

  #===============================================================================
  # STEP 5: Build output
  #-------------------------------------------------------------------------------
  parameter.output <- matrix(NA, nrow = 100, ncol = 5);
  colnames(parameter.output) <- c("nr", "mu", "sigma", "nu", "criterion")

  po <- 1
  parameter.output[po,] <- c(po, mu, sigma, nu, init.model.output.crit)

  po <- po + 1
  if (init.model.output.crit != "NA"){
    {parameter.output[po,] <- c(po, mu, sigma, nu,
                                round(init.model.output.crit,3))}
  }
  if (init.model.output.crit == "NA"){
    {parameter.output[po,] <- c(po, mu, sigma, nu, init.model.output.crit)}
  }
  found.parameters.full <- parameter.output[, 2:4]
  ind <- !is.na(found.parameters.full)
  found.parameters <- {tapply(found.parameters.full[ind],
                              col(found.parameters.full)[ind], tail, 1)}
  names(found.parameters) <- colnames(found.parameters.full)

  if(best_model == "initial") {
    model_result <- init.model.output
  } else {
    best_model_split <- strsplit(best_model, "\\.")[[1]]
    model_result <- eval(parse(text = paste(best_model_split[1], "model", best_model_split[2], "output", sep = ".", collapse = "")))
  }

  return(model_result)
} # end function

#' two_parameters_pb_splines
#'
#' This function is in inner function and conducts the free order procedure
#' with penalized splines for PDFs with two parameters.
#'
#' This function is used for penalized splines
#'
#' @param model A formula that defines the model for continuous norming.
#' @param data Defines the data set.
#' @param family A character that defines the specific PDF for modelling.
#' @param criterion Either "BIC", "AIC" or "GAIC(3)".
#'
#' @author Julian Urban
#'
#' @import gamlss
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @importFrom stats GAIC
#' @importFrom stats poly
#'
#' @return The gamlss model selected by the free order procedure.
#' @export
#'
#' @examples
#' # Apply the free order procedure for a BCPE and the BIC.
#' two_parameters_pb_splines(model = g ~ Normgruppe,
#'                            data = Data_norm,
#'                            criterion = "BIC",
#'                            family = "BCPE")
#'
two_parameters_pb_splines <- function(model, ..., data, criterion, family) {

  #===============================================================================
  # STEP 1: Preperation
  #-------------------------------------------------------------------------------
  model_string <- deparse(model)
  variables <- gsub(" ", "", strsplit(model_string, "~")[[1]])
  data_use <- data[rowSums(is.na(data[, variables])) == 0, variables]

  n <- nrow(data_use)
  if (criterion == "AIC"){
    k <- 2
  } else if (criterion == "BIC"){
    k <- log(n)
  } else if (criterion == "GAIC(3)") {
    k <- 3
  }

  if(min(data_use[, variables[1]]) < 0) {
    data_use[, variables[1]] <- data_use[, variables[1]] - min(data_use[, variables[1]]) + 0.000001
  }

  best_model <- "initial"

  mu <- 1
  sigma <- 1
  i.mu <- mu
  i.sigma <- sigma

  param <- c("1", paste("gamlss::pb(", variables[2],",method = 'GAIC', k =", k, ")"))

  #===============================================================================
  # STEP 2: Fit initial model & prepare free order procedure
  #-------------------------------------------------------------------------------
  init.model <- {paste("","gamlss::gamlss(", model_string, ",sigma.formula = ~ 1,
                     nu.formula = ~ 1, tau.formula = ~ 1, family = family",
                       ", data = data_use,",""," method = RS(10000))",
                       sep = "")}
  init.model.output <- str_eval(noquote(init.model))

  if (criterion == "AIC"){
    init.model.output.crit <- AIC(init.model.output)
  } else if (criterion == "BIC"){
    init.model.output.crit <- BIC(init.model.output)
  } else if (criterion == "GAIC(3)"){
    init.model.output.crit <- GAIC(init.model.output, k = 3)
  }

  parameter.output <- matrix(NA, nrow = 100, ncol = 4);
  colnames(parameter.output) <- c("nr", "mu", "sigma", "criterion")

  po <- 1
  parameter.output[po,] <- c(po, mu, sigma, init.model.output.crit);

  endloop <- "no"

  #===============================================================================
  # STEP 3: Main part of free order procedure
  #-------------------------------------------------------------------------------
  while (endloop != "yes") {

    output <- rep(NA, times = 9)

    #===============================================================================
    # STEP 3.1: Forward model selection
    #-------------------------------------------------------------------------------
    if(mu == 1) {
      i.mu <- 2
    }
    if(sigma == 1) {
      i.sigma <- 2
    }

    # Define the 4 forward models
    mu.model.forward <- {paste('',
                               'gamlss::gamlss(', variables[1],
                               '~', param[i.mu] ,
                               ', sigma.formula = ~ ', param[sigma],
                               ', family = family ',
                               ', data = data_use,',
                               "",
                               'method = RS(10000))',
                               sep = '')}

    sigma.model.forward <- {paste('',
                                  'gamlss::gamlss(', variables[1],
                                  '~', param[mu] ,
                                  ', sigma.formula = ~ ', param[i.sigma],
                                  ', family = family ',
                                  ', data = data_use,',
                                  "",
                                  'method = RS(10000))',
                                  sep = '')}

    mu.model.forward.output <- str_eval(noquote(mu.model.forward))
    sigma.model.forward.output <- str_eval(noquote(sigma.model.forward))

    if(criterion == "AIC"){
      mu.model.forward.output.crit <- AIC(mu.model.forward.output)
      sigma.model.forward.output.crit <- AIC(sigma.model.forward.output)
    } else if(criterion == "BIC"){
      mu.model.forward.output.crit <- BIC(mu.model.forward.output)
      sigma.model.forward.output.crit <- BIC(sigma.model.forward.output)
    } else if(criterion == "GAIC(3)") {
      mu.model.forward.output.crit <- GAIC(mu.model.forward.output, k = 3)
      sigma.model.forward.output.crit <- GAIC(sigma.model.forward.output, k = 3)
    }

    #===============================================================================
    # STEP 3.2: Backward model selection
    #-------------------------------------------------------------------------------
    if(mu == 2) {
      i.mu <- 1
    }
    if(sigma == 2) {
      i.sigma <- 1
    }

    if(mu == 2) {
      i.mu <- 1
      mu.model.backward <- {paste('',
                                  'gamlss::gamlss(', variables[1],
                                  '~', param[i.mu] ,
                                  ', sigma.formula = ~ ', param[sigma],
                                  ', family = family ',
                                  ', data = data_use,',
                                  "",
                                  'method = RS(10000))',
                                  sep = '')}

      mu.model.backward.output <- str_eval(noquote(mu.model.backward))

      if (criterion == "AIC"){
        mu.model.backward.output.crit <- AIC(mu.model.backward.output)
      } else if (criterion == "BIC"){
        mu.model.backward.output.crit <- BIC(mu.model.backward.output)
      } else if (criterion == "GAIC(3)"){
        mu.model.backward.output.crit <- GAIC(mu.model.backward.output, k = 3)
      }
    } else {
      mu.model.backward.output.crit <- Inf
    }

    if(sigma == 2) {
      i.sigma <- 1
      sigma.model.backward <- {paste('',
                                     'gamlss::gamlss(', variables[1],
                                     '~', param[mu] ,
                                     ', sigma.formula = ~ ', param[i.sigma],
                                     ', family = family ',
                                     ', data = data_use,',
                                     "",
                                     'method = RS(10000))',
                                     sep = '')}

      sigma.model.backward.output <- str_eval(noquote(sigma.model.backward))

      if (criterion == "AIC"){
        sigma.model.backward.output.crit <- AIC(sigma.model.backward.output)
      } else if (criterion == "BIC"){
        sigma.model.backward.output.crit <- BIC(sigma.model.backward.output)
      } else if (criterion == "GAIC(3)"){
        sigma.model.backward.output.crit <- GAIC(sigma.model.backward.output, k = 3)
      }
    } else {
      sigma.model.backward.output.crit <- Inf
    }


    #===============================================================================
    # STEP 4: Output creation & Evaluation
    #-------------------------------------------------------------------------------
    output <- {c(init.model.output.crit,
                 mu.model.forward.output.crit, sigma.model.forward.output.crit,
                 mu.model.backward.output.crit, sigma.model.backward.output.crit)}
    names(output) <- {c("initial",
                        "mu.forward", "sigma.forward",
                        "mu.backward", "sigma.backward")}


    if ({output["mu.forward"] < output["initial"] &
        output["mu.forward"] < output["sigma.forward"] &
        output["mu.forward"] < output["mu.backward"] &
        output["mu.forward"] < output["sigma.backward"]})
    {mu <- 2
    init.model.output.crit <- mu.model.forward.output.crit
    best_model <- "mu.forward"
    } else {
      if ({output["sigma.forward"] < output["initial"] &
          output["sigma.forward"] < output["mu.forward"] &
          output["sigma.forward"] < output["mu.backward"] &
          output["sigma.forward"] < output["sigma.backward"]}) {
        sigma <- 2
        init.model.output.crit <- sigma.model.forward.output.crit
        best_model <- "sigma.forward"
      } else {
        if ({output["mu.backward"] < output["initial"] &
            output["mu.backward"] < output["mu.forward"] &
            output["mu.backward"] < output["sigma.forward"] &
            output["mu.backward"] < output["sigma.backward"] }) {
          mu <- 1
          init.model.output.crit <- mu.model.backward.output.crit
          best_model <- "mu.backward"
        } else {
          if ({output["sigma.backward"] < output["initial"] &
              output["sigma.backward"] < output["mu.forward"] &
              output["sigma.backward"] < output["sigma.forward"] &
              output["sigma.backward"] < output["mu.backward"]}) {
            sigma <- 1
            init.model.output.crit <- sigma.model.backward.output.crit
            best_model <- "sigma.backward"
          } else {
            endloop <- "yes" }}}}
  } # end while loop

  #===============================================================================
  # STEP 5: Build output
  #-------------------------------------------------------------------------------
  parameter.output <- matrix(NA, nrow = 100, ncol = 4);
  colnames(parameter.output) <- c("nr", "mu", "sigma", "criterion")

  po <- 1
  parameter.output[po,] <- c(po, mu, sigma, init.model.output.crit)

  po <- po + 1
  if (init.model.output.crit != "NA"){
    {parameter.output[po,] <- c(po, mu, sigma,
                                round(init.model.output.crit,3))}
  }
  if (init.model.output.crit == "NA"){
    {parameter.output[po,] <- c(po, mu, sigma, init.model.output.crit)}
  }

  found.parameters.full <- parameter.output[, 2:3]
  ind <- !is.na(found.parameters.full)
  found.parameters <- {tapply(found.parameters.full[ind],
                              col(found.parameters.full)[ind], tail, 1)}
  names(found.parameters) <- colnames(found.parameters.full)

  if(best_model == "initial") {
    model_result <- init.model.output
  } else {
    best_model_split <- strsplit(best_model, "\\.")[[1]]
    model_result <- eval(parse(text = paste(best_model_split[1], "model", best_model_split[2], "output", sep = ".", collapse = "")))
  }

  return(model_result)
} # end function








