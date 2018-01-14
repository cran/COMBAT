\name{ld.Rsquare}
\alias{ld.Rsquare}
\title{
 Linkage Disequilibrium Estimate as Correlation
}
\description{
  Compute linkage disequilibrium among SNPs using correlation coefficients.
}
\usage{
ld.Rsquare(x)
}
\arguments{
 \item{x}{a matrix of SNP genotypes with samples in the rows.}
}
\value{
A positive definite correlation matrix.
}
\seealso{
\code{\link{gates}}, \code{\link{vegas}}, \code{\link{simpleM}}.
}
