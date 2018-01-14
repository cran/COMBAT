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
