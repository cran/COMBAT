\name{ext_simes}
\alias{ext_simes}
\title{
  Extended Simes Procedure
}
\description{
  Combine a vector of test P values by correction for number of independent tests.
}
\usage{
ext_simes(x, cor_r)
}
\arguments{
 \item{x}{a vector of SNP-level P values.}
 \item{cor_r}{correlation among P values.}
}
\value{
P value.
}
\seealso{
\code{\link{COMBAT}}
}
\examples{
# see ?COMBAT
}
