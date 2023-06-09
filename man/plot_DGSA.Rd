% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DGSA.R
\name{plot_DGSA}
\alias{plot_DGSA}
\title{plot DGSA}
\usage{
plot_DGSA(selected_matrix_disease_component, case_tag)
}
\arguments{
\item{selected_matrix_disease_component}{Disease component matrix of
the selected genes that contains the disease component of all patients.
Output of the function \code{generate_disease_component}.}

\item{case_tag}{Character vector of the same length as the number of
columns of full_data. Patients must be in the same order as in full_data.
It must be indicated for each patient whether he/she is healthy or not.
One value should be used to indicate whether the patient is healthy and
another value should be used to indicate whether the patient's sample is
tumourous. The user will then be asked which one indicates whether
the patient is healthy. Only two values are valid in the vector in total.}
}
\value{
The heatmap of the DGSA result.
}
\description{
It draws the heatmap of the DGSA result by selecting the 100 genes with
the highest variability between samples.
}
