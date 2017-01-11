\name{ext_simes}
\alias{ext_simes}
\alias{gates}
\alias{vegas}
\title{
  Calling Gene-based Association Tests
}
\description{
   Calling Gene-based Association Tests.
}
\usage{
ext_simes(x, cor_r)
gates(x, cor_G)
vegas(x, cor_G, vegas_vec, n_simul)
}
\arguments{
 \item{x}{a vector of SNP-level P values.}
 \item{cor_r}{correlation among P values.}
 \item{cor_G}{SNP-SNP correlation matrix.}
 \item{vegas_vec}{a numeric vector, specifying the fraction of the top SNPs to be used in the VEGAS method.}
 \item{n_simul}{number of simulations.}
}
\details{
 Function ext_simes implements extended Simes method, gates implements the GATES method, and vegas implements VEGAS with different proportion tests.
}
\value{
P value(s).
}
\references{
 Jianfei Huang, Minghui Wang, James B. Potash, Shizhong Han. COMBAT: A Combined Association Test for Genes using Summary Statistics. Manuscript submitted.
}
\seealso{
\code{\link{COMBAT}}
}
\examples{
# read SNP P values
file1 <- paste(path.package("COMBAT"),"extdata","SNP_info.txt.gz",sep="/")
snp.info  <- read.table(file1, header = TRUE, as.is=TRUE)
snp.pvals <- as.matrix(snp.info[,2])

# read reference genotype
file2 <- paste(path.package("COMBAT"),"extdata","SNP_ref.txt.gz",sep="/")
snp.ref   <- read.table(file2, header = TRUE)
snp.ref   <- as.matrix(snp.ref)

#call COMBAT; not run
#COMBAT(snp.pvals, snp.ref, nperm=50)

}
