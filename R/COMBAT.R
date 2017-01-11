# ==============================================================================================
# Function of extended Simes
ext_simes = function(x, cor_r){
	eff.snpcount.fun <- function(ldmat) {
		ldmat <- as.matrix(ldmat)
		snpcount.local <- dim(ldmat)[1]
		if (snpcount.local <= 1) return(1)
		ev <- eigen(ldmat, only.values = TRUE)$values
		if (sum(ev < 0) != 0) {
				ev <- ev[ev > 0]
				ev <- ev/sum(ev) * snpcount.local
		}
		ev <- ev[ev > 1]
		snpcount.local - sum(ev - 1)
	}
	eff.snpcount.global <- eff.snpcount.fun(cor_r)
	
	n_values <- length(x)
	candid <- sapply(1:n_values, function(i){
		(eff.snpcount.global * x[i])/eff.snpcount.fun(cor_r[1:i,1:i])
	
	})
	
	p_ext_simes <- min(candid)
	p_ext_simes
}

# Function of GATES
gates = function(x, cor_G){
	pval_sort  <- sort(x)
	pval_order <- order(x)
	n_snps     <- length(x)

	cor_P <- cor_G[pval_order, pval_order]
	cor_P <- 0.2982*cor_P^6 - 0.0127*cor_P^5 + 0.0588*cor_P^4 + 0.0099*cor_P^3 + 0.6281*cor_P^2 - 0.0009*cor_P
        
	p_gates <- ext_simes(pval_sort, cor_P)
	p_gates
}

# Function of VEGAS with different proportion tests and Fisher combination test
vegas = function(x, cor_G, vegas_vec, n_simul){
	chisq_vec <- qchisq(x,1,lower.tail=FALSE)
	chisq_vec[chisq_vec == Inf] <- 60
	n_snps <- length(x)
	n_tests <- length(vegas_vec)

	TS_obs <- rep(NA,n_tests)
	TS_obs[1] <- max(chisq_vec, na.rm=TRUE)
	for (j in 2:n_tests){
		k <- vegas_vec[j]
		TS_obs[j] <- sum(sort(chisq_vec, decreasing = TRUE)[1:k])
	}

	rd  <- rmvnorm(n_simul, mean=rep(0,n_snps),sigma=cor_G)
	rd2 <- rd^2

	pPerm0 <- rep(NA,n_tests)
	T0s <- apply(rd2,1,max)
	pPerm0[1]<- (sum(T0s >= TS_obs[1])+1)/(length(T0s)+1)
	for(j in 2:n_tests){
		k <- vegas_vec[j]
		for (i in 1:n_simul){
			T0s[i] <- sum(sort(rd2[i,],decreasing = TRUE)[1:k])
		}
		pPerm0[j] <- (sum(T0s >= TS_obs[j])+1)/(length(T0s)+1)
	}
	pPerm0
}

# SimpleM
simpleM = function(x, cor_G){
	min_p_obs <- min(x)
	num_of_snps <- length(x)
	cor_r <- cor_G
	
	eigen_values <- eigen(cor_r, only.values = TRUE)$values
	eigen_values_sorted <- sort(eigen_values, decreasing = TRUE)
	
	pca_cut_perc <- 0.99
	sum_eigen_values <- sum(eigen_values_sorted)
	
	M_eff_G <- 1
	for(k in 1:num_of_snps){
		temp <- sum(eigen_values_sorted[1:k])/sum_eigen_values
		if(temp >= pca_cut_perc){
			M_eff_G <- k
			break
		}
	}
	p_simpleM <- 1-(1-min_p_obs)^M_eff_G
}

# main function
COMBAT = function(x, snp.ref, vegas.pct = c(0.1,0.2,1), nperm = 75, seed=12345, ncores=1){
	pvalues <- as.numeric(x)
	n_snps <- length(pvalues)
	cor_G <- cor(as.matrix(snp.ref))
	vegas_vec <- c(1, ceiling(vegas.pct*n_snps))
	
	if(is.positive.definite(cor_G)==FALSE){
		cor_G <- make.positive.definite(cor_G)
	}
	if(is.positive.definite(cor_G)==FALSE){
		cor_G <- cor(as.matrix(snp.ref))
		diag(cor_G) <- 1.0001
	}
	if(is.positive.definite(cor_G)==FALSE){
		diag(cor_G) <- 1.001
	}
	if(is.positive.definite(cor_G)==FALSE){
		diag(cor_G) <- 1.01
	}
	
	set.seed(seed)
	# carry out gates, vegas, and SimpleM 
	pval_gates <- gates(pvalues, cor_G)
	pval_vegas <- vegas(pvalues, cor_G, vegas_vec, 1000)
	if(sum(pval_vegas <= 0.005) >= 1){
		pval_vegas <- vegas(pvalues, cor_G, vegas_vec, 10000)
	}
	if(sum(pval_vegas <= 0.0005) >= 1){
		pval_vegas <- vegas(pvalues, cor_G, vegas_vec, 100000)
	}
	if(sum(pval_vegas <= 0.00005) >= 1){
		pval_vegas <- vegas(pvalues, cor_G, vegas_vec, 1000000)
	}
	
	pval_simpleM <- simpleM(pvalues, cor_G)
	gene_pvals <- c(pval_gates, pval_vegas, pval_simpleM)
	
	# compute p-value correlations matrix of different tests 
	rd  <- rmvnorm(nperm, mean=rep(0,n_snps),sigma=cor_G)
	rd2 <- rd^2
	simul_pval_mat <- pchisq(rd2,1,lower.tail=FALSE) 
	func1=function(x,cor_G,vegas_vec,vegas_samples=1000){
		p_gates <- gates(x, cor_G)
		p_vegas <- vegas(x, cor_G, vegas_vec,vegas_samples)
		p_simpleM <- simpleM(x, cor_G)
		c(p_gates, p_vegas, p_simpleM)
	}
	if(ncores>1 && requireNamespace("parallel",quietly = TRUE)){
		cl=parallel::makeCluster(ncores)
		parallel::clusterSetRNGStream(cl, .Random.seed)
		gene_pval_mat=parallel::parApply(cl,simul_pval_mat,1,func1,cor_G=cor_G,vegas_vec=vegas_vec)
		parallel::stopCluster(cl)
	}else{
		gene_pval_mat=apply(simul_pval_mat,1,func1,cor_G=cor_G,vegas_vec=vegas_vec)
	}
	gene_pval_mat=t(gene_pval_mat)
	method_cor <- cor(gene_pval_mat)
	
	# compute the COMBAT P-value from different methods	
	sort_pvals <- sort(gene_pvals)
	order_pvals <- order(gene_pvals)
	method_cor <- method_cor[order_pvals, order_pvals]
	p_combat <- ext_simes(sort_pvals, method_cor)
	
	res <- c(gene_pvals,p_combat)
	names(res) <- c("GATES","VEGAS-max","VEGAS-10%","VEGAS-20%","VEGAS-sum","SimpleM","COMBAT")	
	res
}
