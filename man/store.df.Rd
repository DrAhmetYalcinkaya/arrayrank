% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/store.df.R
\name{store.df}
\alias{store.df}
\title{Store Data Frame for Export}
\usage{
store.df(data)
}
\arguments{
\item{data}{A data frame with protein expression data where rows represent proteins and columns represent different samples or calculated values. The row names should contain the protein names.}
}
\value{
A data frame with an additional column for protein names, ready for export to Excel or other formats.
}
\description{
This function takes a data frame, adds protein names as a separate column, and prepares it for export, typically into an Excel file.
}
\examples{
\dontrun{
# Example usage of store.df function
data_to_store <- store.df(data = your_data)
}
}
