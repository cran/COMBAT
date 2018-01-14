# Function of VEGAS with different proportion tests and Fisher combination test
vegas = function(x, cor_G, vegas.pct=c(0.1,0.2,0.3,0.4,1), max.simulation=1e6){
	if(is.positive.definite(cor_G)==FALSE) stop('cor_G is not positive definite. Please re-calculate with ld.Rsqure function.\n')
	pval_vegas <- vegas.call(x=x, cor_G=cor_G, vegas.pct=vegas.pct, n_simul=1000)
	if(any(pval_vegas <= 0.005)){
		pval_vegas <- vegas.call(x=x, cor_G=cor_G, vegas.pct=vegas.pct, n_simul=10000)
	}
	if(any(pval_vegas <= 0.0005)){
		pval_vegas <- vegas.call(x=x, cor_G=cor_G, vegas.pct=vegas.pct, n_simul=100000)
	}
	if(any(pval_vegas <= 5e-5)){
		pval_vegas <- vegas.call(x=x, cor_G=cor_G, vegas.pct=vegas.pct, n_simul=1000000)
	}
	if(any(pval_vegas <= 5e-6) && max.simulation>1000000){
		pval_vegas <- vegas.call(x=x, cor_G=cor_G, vegas.pct=vegas.pct, n_simul=max.simulation)
	}
	pval_vegas
}
#this is not to be called directly
vegas.call = function(x, cor_G, vegas.pct, n_simul){
	stopifnot(length(x) == ncol(cor_G))
	vegas_vec <- ceiling(vegas.pct*ncol(cor_G))
	vegas_vec <- sort(vegas_vec)
	if(vegas_vec[1]>1){
		vegas.pct <- c(0,vegas.pct)
		vegas_vec <- c(1,vegas_vec)
	}
	chisq_vec <- qchisq(x,1,lower.tail=FALSE)
	chisq_vec[chisq_vec == Inf] <- 60
	n_snps <- length(x)
	n_tests <- length(vegas_vec)

	TS_obs <- rep(NA,n_tests)
	TS_obs[1] <- max(chisq_vec, na.rm=TRUE)
	chisq_vec <- sort(chisq_vec, decreasing = TRUE)
	for (j in 2:n_tests) TS_obs[j] <- sum(chisq_vec[1:vegas_vec[j]])

	rd  <- rmvnorm(n_simul, mean=rep(0,n_snps),sigma=cor_G)
	rd2 <- rd^2
	rd2 <- apply(rd2,1,sort,decreasing=TRUE)

	pPerm0 <- rep(NA,n_tests)
	T0s <- apply(rd2,2,max)
	pPerm0[1]<- (sum(T0s >= TS_obs[1])+1)/(length(T0s)+1)
	for(j in 2:n_tests){
		for (i in 1:n_simul) T0s[i] <- sum(rd2[1:vegas_vec[j],i])
		pPerm0[j] <- (sum(T0s >= TS_obs[j])+1)/(length(T0s)+1)
	}
	v1 <- paste0('VEGAS.p',vegas.pct)
	v1[vegas_vec==ncol(cor_G)] <- 'VEGAS.all'
	v1[vegas_vec==1] <- 'VEGAS.max'
	names(pPerm0) <- v1
	pPerm0
}
