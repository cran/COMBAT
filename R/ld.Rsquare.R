#compute ld matrix using correlation
ld.Rsquare=function(x){
	cor_G <- cor(as.matrix(x),use='p')
	if(is.positive.definite(cor_G)==FALSE){
		cor_G <- make.positive.definite(cor_G)
	}
	if(is.positive.definite(cor_G)==FALSE){
		cor_G <- cor(as.matrix(x),use='p')
		diag(x) <- 1.0001
	}
	if(is.positive.definite(cor_G)==FALSE){
		diag(cor_G) <- 1.001
	}
	if(is.positive.definite(cor_G)==FALSE){
		diag(cor_G) <- 1.01
	}
	cor_G
}
