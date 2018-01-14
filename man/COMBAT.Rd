\name{COMBAT}
\alias{COMBAT}
\title{
  A Combined Gene-based Association Test
}
\description{
   This function implements a combined gene-based association test using SNP-level P values and reference genotype data.
}
\usage{
COMBAT(x, snp.ref, vegas.pct = c(0.1,0.2,0.3,0.4,1), pca_cut_perc = 0.995,
   nperm = 100, seed=12345, ncores=1)
}
\arguments{
 \item{x}{a vector of SNP-level P values.}
 \item{snp.ref}{a matrix of SNP genotypes (coded as allele counts) from reference samples, with samples in rows and SNPs in columns.}
 \item{vegas.pct}{a numeric vector, fraction of the top SNPs to be used in the VEGAS method.}
 \item{pca_cut_perc}{numeric, cutoff for percentage of sum of eigen values in the simpleM approach.}
 \item{nperm}{number of permutations for computing the correlation between P values of different tests.}
 \item{seed}{random seed to derive consistent outcome.}
 \item{ncores}{number of CPU cores for parallel computing.}
}
\details{
 COMBAT uses simulation and the extended Simes procedure (\code{\link{ext_simes}}) to combine multiple gene-based association test statistics (currently including \code{\link{gates}}, \code{\link{vegas}}, and \code{\link{simpleM}}) to perform a more powerful association analysis. This method does not require raw genotype or phenotype data, but needs only SNP-level P-values and correlations between SNPs from ancestry-matched samples. The technical details about the method is described in Wang et al (2017) <doi:10.1534/genetics.117.300257>.
}
\value{
A vector of p-values from COMAT and each individual gene-based test.
}
\references{
 Minghui Wang, Jianfei Huang, Yiyuan Liu, Li Ma, James B. Potash, Shizhong Han. COMBAT: A Combined Association Test for Genes using Summary Statistics. Genetics 2017, 207(3): 883-891. https://doi.org/10.1534/genetics.117.300257.
}
\seealso{
\code{\link{ext_simes}}, \code{\link{gates}}, \code{\link{vegas}}, \code{\link{simpleM}}.
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
#call COMBAT
\donttest{
COMBAT(snp.pvals, snp.ref, nperm=100, ncores=2)
}
\dontshow{
COMBAT(snp.pvals, snp.ref, nperm=10, ncores=1)
}
}
