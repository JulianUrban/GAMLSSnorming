% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamlss_pb.R
\name{four_parameters_pb_splines}
\alias{four_parameters_pb_splines}
\title{four_parameters_pb_splines}
\usage{
four_parameters_pb_splines(model, ..., data, criterion, family)
}
\arguments{
\item{model}{A formula that defines the model for continuous norming.}

\item{data}{Defines the data set.}

\item{criterion}{Either "BIC", "AIC" or "GAIC(3)".}

\item{family}{A character that defines the specific PDF for modelling.}
}
\value{
The gamlss model selected by the free order procedure.
}
\description{
This function is in inner function and conducts the free order procedure
with penalized splines for PDFs with four parameters.
}
\details{
This function is used for penalized splines
}
\examples{
# Apply the free order procedure for a BCPE and the BIC.
four_parameters_pb_splines(model = g ~ Normgruppe,
                           data = Data_norm,
                           criterion = "BIC",
                           family = "BCPE")

}
\author{
Julian Urban
}
