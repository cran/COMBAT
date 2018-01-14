# COMBAT function
COMBAT = function(x, snp.ref, vegas.pct = c(0.1,0.2,0.3,0.4,1), pca_cut_perc=0.995, nperm = 100, seed=12345, ncores=1){
	pvalues <- as.numeric(x)
	n_snps <- length(pvalues)
	#
	cor_G <- ld.Rsquare(snp.ref)
	set.seed(seed)
	#
	pval_gates <- gates(x=pvalues, cor_G=cor_G)
	pval_vegas <- vegas(x=pvalues, cor_G=cor_G, vegas.pct=vegas.pct)
	pval_simpleM <- simpleM(x=pvalues, cor_G=cor_G, pca_cut_perc=pca_cut_perc)
	gene_pvals <- c(GATES=pval_gates, pval_vegas, simpleM=pval_simpleM)
	#
	# compute p-value correlation matrix
	rd  <- rmvnorm(nperm, mean=rep(0,n_snps),sigma=cor_G)
	rd2 <- rd^2
	simul_pval_mat <- pchisq(rd2,1,lower.tail=FALSE) 
	func1=function(x,cor_G,vegas.pct,pca_cut_perc=0.995){
		p_gates <- gates(x=x, cor_G=cor_G)
		p_vegas <- vegas(x=x, cor_G=cor_G, vegas.pct=vegas.pct)
		p_simpleM <- simpleM(x=x, cor_G=cor_G, pca_cut_perc=pca_cut_perc)
		c(p_gates, p_vegas, p_simpleM)
	}
	if(ncores>1 && requireNamespace("parallel",quietly = TRUE)){
		cl=parallel::makeCluster(ncores)
		parallel::clusterSetRNGStream(cl, .Random.seed)
		gene_pval_mat=parallel::parApply(cl,simul_pval_mat,1,func1,cor_G=cor_G,vegas.pct=vegas.pct,pca_cut_perc=pca_cut_perc)
		parallel::stopCluster(cl)
	}else{
		gene_pval_mat=apply(simul_pval_mat,1,func1,cor_G=cor_G,vegas.pct=vegas.pct,pca_cut_perc=pca_cut_perc)
	}
	gene_pval_mat=t(gene_pval_mat)
	method_cor <- cor(gene_pval_mat)
	
	# combat by simpleM
	#p_combat_simpleM <- simpleM(x=gene_pvals, cor_G=method_cor, pca_cut_perc=pca_cut_perc)
	
	# compute the COMBAT P-value using the extended Simes procedure	
	order_pvals <- order(gene_pvals)
	sort_pvals <- gene_pvals[order_pvals]
	method_cor <- method_cor[order_pvals, order_pvals]
	p_combat_simes <- ext_simes(sort_pvals, method_cor)

	#
	res <- c(COMBAT=p_combat_simes,gene_pvals)
	res
}
