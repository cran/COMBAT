\name{gates}
\alias{gates}
\alias{vegas}
\alias{simpleM}
\title{
  Calling Gene-based Association Tests
}
\description{
   Several gene-based association tests methods are implemented.
}
\usage{
gates(x, cor_G)
vegas(x, cor_G, vegas.pct=c(0.1,0.2,0.3,0.4,1), max.simulation=1e6)
simpleM(x, cor_G, pca_cut_perc=0.995)
}
\arguments{
 \item{x}{a vector of SNP-level P values.}
 \item{cor_G}{SNP-SNP correlation matrix.}
 \item{vegas.pct}{a numeric vector, specifying the fraction of the top SNPs to be used in the VEGAS method.}
 \item{max.simulation}{maximum number of simulations to be performed. Must be at least 1e6.}
 \item{pca_cut_perc}{cutoff for percentage of sum of eigen values.}
}
\details{
 Function \code{gates} implements the GATES method (Li et al 2011, American Journal of Human Genetics 88:283-293), \code{vegas} implements VEGAS with different proportion tests (Liu et al 2010, American Journal of Human Genetics 87:139-145), and \code{simpleM} is the simpleM method (Gao et al 2008, Genetic Epidemiology 32:361-369).
}
\value{
P value(s).
}
\seealso{
\code{\link{COMBAT}}, \code{\link{ld.Rsquare}}.
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

#compute correlation among SNPs
cor_G <- ld.Rsquare(snp.ref)

#call gates
(pval_gates <- gates(x=snp.pvals, cor_G=cor_G))

#call vegas
(pval_vegas <- vegas(x=snp.pvals, cor_G=cor_G))

#call simpleM
(pval_simpleM <- simpleM(x=snp.pvals, cor_G=cor_G))
}
