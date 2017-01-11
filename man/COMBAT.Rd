\name{COMBAT}
\alias{COMBAT}
\title{
  A Combined Gene-based Association Test
}
\description{
   This function implements a combined gene-based association test using SNP-level P valuese and reference genotype data.
}
\usage{
COMBAT(x, snp.ref, vegas.pct = c(0.1,0.2,1), nperm = 75,
   seed=12345, ncores=1)
}
\arguments{
 \item{x}{a vector of SNP-level P values.}
 \item{snp.ref}{a matrix of SNP genotypes (coded as allele counts) from reference samples, with samples in rows and SNPs in columns.}
 \item{vegas.pct}{a numeric vector, specifying the fraction of the top SNPs to be used in the VEGAS method.}
 \item{nperm}{number of permutations for computing the correlation between P values of different tests.}
 \item{seed}{random seed to derive consistent outcome.}
 \item{ncores}{number of CPU cores for parallel computing.}
}
\details{
 Some technical details about the method.
}
\value{
A vector of p-values from differrnt gene-based tests with summary statistics. By default, 7 tests will be computed: GATES, 4 VEGAS tests, SimpleM and one COMBAT.
}
\references{
 Jianfei Huang, Minghui Wang, James B. Potash, Shizhong Han. COMBAT: A Combined Association Test for Genes using Summary Statistics. Manuscript submitted.
}
\seealso{
\code{\link{ext_simes}}, \code{\link{gates}}, \code{\link{vegas}}.
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
