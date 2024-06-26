% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{survival_event}
\alias{survival_event}
\title{Survival event vector}
\format{
Character vector of length 121.
\describe{A value of "0" indicates that the patient did not relapse
during follow-up, a value of "1" indicates that the patient did. Samples
from healthy tissue contain a value of \code{NA}.
}
}
\source{
The data are from the study GSE42568. Information extracted from
the file GSE42568_family.soft.gz available at
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568}.
}
\usage{
data(survival_event, package = "GSSTDA")
}
\description{
Character vector of length 121 containing whether or not the patient
is relapsed.
}
\keyword{datasets}
