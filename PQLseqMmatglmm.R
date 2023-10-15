######################################################################
# Authors: Ping Zeng                                                 #
# Contact: zpstat@xzhmu.edu.cn                                       #
#          Department of Biostatistics, XuZhou Medical University    #
######################################################################
## this codes are modified and constructed based the package glmmPQL of Shiquan Sun and Xiang Zhou
#binomial(link = "logit")
#gaussian(link = "identity")
#Gamma(link = "inverse")
#inverse.gaussian(link = "1/mu^2")
#poisson(link = "log")
#quasi(link = "identity", variance = "constant")
#quasibinomial(link = "logit")
#quasipoisson(link = "log")

glmmPQL <- function(
	y,
	X       = NULL, # not include the intercept term, 会自动添加
	N       = NULL,
	offset  = NULL,
	kmat    = NULL, # n by n matrix
	emat    = TRUE, # for only PMM and BMM；即，是否在PMM and BMM中包含残差项
	model   = "BMM",# LMM, PMM and BMM
	method  = "AI", # AI, EI and OI
	maxiter = 100,
	tol     = 1e-5,
	verbose = TRUE, ...) {

	if (!model %in% c("LMM", "BMM", "PMM"))
		{stop ("Error: \"model\" must be \"LMM\" or \"BMM\" or \"PMM\".")}
	if (!method %in% c("AI", "EI", "OI"))
		{stop ("Error: \"method\" must be \"AI\" or \"EI\" or \"OI\".")}

	numIDV <- length(y)
	if (is.null(X)) {numCov <- 0}
	else {X <- as.matrix(X); numCov <- dim(X)[2]}
	if (is.null(kmat)) {stop ("glmmPQL::please input genetic relatedness matrix!")}

	if (class(kmat) == "matrix") {
		#if (numIDV!=dim(kmat)[1]) {
		#stop ("glmmPQL::the dimention of relatedness matrix doesnot match the number of individuals!")
		#}
		#eig          <- eigen(kmat)
		#eigval       <- eig$value
		#if (any(eigval<1e-10)){
		#	warning("glmmPQL::the relatedness matrix is singular, it has been modified!")
		#	kmat <- as.matrix(nearPD(kmat, corr=F)$mat)
		#}
		#rm(eig)
		#rm(eigval)
		if ((emat) & (model!="LMM")) {kmat <- list(kmat, diag(numIDV))}# LMM时不再加入单位矩阵
		else                         {kmat <- list(kmat)}
		numK <- length(kmat)
		}

	else {
	for (ik in 1:length(kmat)) {
		kmat[[ik]] <- as.matrix(kmat[[ik]])
		if (numIDV!=dim(kmat[[ik]])[1]){stop ("glmmPQL::the dimention of relatedness matrix 
						       is not equal to the number of individuals!")}
		#eig               <- eigen(kmat[[ik]])
		#eigval            <- eig$value
		#if (any(eigval<1e-10)){ 
		#	warning("glmmPQL::the relatedness matrix is singular, it has been modified!")
		#	kmat[[ik]] <- as.matrix(nearPD(kmat[[ik]], corr=F)$mat)
		#	}
		#rm(eig)
		#rm(eigval)
		}
		kmat_new <- list()
		for (ik in 1:(length(kmat))) {kmat_new[[ik]]   <- kmat[[ik]]}
		if ((emat) & (model!="LMM")) {kmat_new[[ik+1]] <- diag(numIDV)}
		kmat <- kmat_new
		numK <- length(kmat)
		rm(kmat_new)
	}

	if (method=="AI") {opti ="Average  Information (AI) REML"}
	if (method=="EI") {opti ="Expected Information (EI) REML"}
	if (method=="OI") {opti ="Observed Information (OI) REML"}
	cat(paste("## number of individuals: ", numIDV,"\n"))
	cat(paste("## number of covariates: ", numCov,"\n"))
	cat(paste("## number of relatedness matrix: ", numK,"\n"))
	cat(paste("## maximum number of iterations: ", maxiter,"\n"))
	cat(paste("## convergence tolerance: ", tol,"\n"))
	cat(paste("## optimization method: ", opti,"\n"))
	cat(paste("..................................................","\n"))

	#***********************************#
	#       Linear Mixed Model         #
	#***********************************#
	if (model == "LMM"){
		cat("## now fitting linear mixed model ... \n")
		if (is.null(offset)) {offset <- rep(0, numIDV)}
		if (numCov==0) {
		model0 <- try(glm(formula = y ~ 1 + offset(offset), family = gaussian(link = "identity")))
		 idx   <- match(rownames(model.frame(formula = y ~  1 + offset(offset), na.action = na.omit)),
				rownames(model.frame(formula = y ~  1 + offset(offset), na.action = na.pass)))
			}
		else {
		model0 <- try(glm(formula = y ~ X + offset(offset), family = gaussian(link = "identity")))
		 idx   <- match(rownames(model.frame(formula = y ~ X + offset(offset), na.action = na.omit)),
				rownames(model.frame(formula = y ~ X + offset(offset), na.action = na.pass)))
		}
		for (ik in 1:(length(kmat))) {kmat[[ik]] <- kmat[[ik]][idx, idx]}
		if (class(model0)[1]!="try-error") { 
		if (method == "AI") {
		model1 <- try(glmmPQL_fit(model0, kmat, method = "AI", maxiter, tol, verbose))}
		else {
		if (method == "EI") {
		model1 <- try(glmmPQL_fit(model0, kmat, method = "EI", maxiter, tol, verbose))}
		else {
		if (method == "OI") {
		model1 <- try(glmmPQL_fit(model0, kmat, method = "OI", maxiter, tol, verbose))}}
		}
		}
		#else {model1 <- NULL}
		if (!is.null(model1) & (class(model1)!="try-error")){
			if (verbose) {cat(paste("glmmPQL::LMM::tau = ", model1$theta,"\n"))}}
		else {converged <- FALSE}
		closeAllConnections()
		return(model1)
	}# end LMM 

	#***********************************#
	#       Poisson Mixed Model         #
	#***********************************#
	if (model == "PMM"){
		cat("# now fitting Poisson Mixed Model ... \n")
		if (is.null(N)) {N <- rep(1,numIDV)}
		else {N <- N}
		if (numCov==0) {
		model0 <- try(glm(formula = y ~ 1 + offset(log(N)), family = poisson(link="log")))
		 idx   <- match(rownames(model.frame(formula = y ~  1 + offset(log(N)), na.action = na.omit)),
				rownames(model.frame(formula = y ~  1 + offset(log(N)), na.action = na.pass)))
			}
		else {
		model0 <- try(glm(formula = y ~ X + offset(log(N)), family = poisson(link="log")))
		 idx   <- match(rownames(model.frame(formula = y ~ X + offset(log(N)), na.action = na.omit)),
				rownames(model.frame(formula = y ~ X + offset(log(N)), na.action = na.pass)))
		}
		for (ik in 1:(length(kmat))) {kmat[[ik]] <- kmat[[ik]][idx, idx]}
		if (class(model0)[1]!="try-error") { 
		if (method == "AI") {
		model1 <- try(glmmPQL_fit(model0, kmat, method = "AI", maxiter, tol, verbose))}
		else {
		if (method == "EI") {
		model1 <- try(glmmPQL_fit(model0, kmat, method = "EI", maxiter, tol, verbose))}
		else {
		if (method == "OI") {
		model1 <- try(glmmPQL_fit(model0, kmat, method = "OI", maxiter, tol, verbose))}}
		}
		}
		#else {model1 <- NULL}
		if (!is.null(model1)&(class(model1)!="try-error")){
			if (verbose) {cat(paste("glmmPQL::PMM::tau = ", model1$theta,"\n"))}}
		else {converged <- FALSE}
		closeAllConnections()
		return(model1)
	}# end PMM 


	#***************************************#
	#          Binary Mixed Model           #
	#***************************************#
	if (model == "BMM"){ 
		cat("## fitting Binomial Mixed Model ... \n")
		if (is.null(N)) {N <- rep(1, numIDV)}
		else {N <- N}
		ratio               <- y/N
		ratio[is.na(ratio)] <- 0
		flag                <- ratio>1.0
		sumflag             <- sum(flag)
		idx                 <- which(sumflag>0)
		if (length(idx)>0){y <- y[-idx]; N <- N[-idx]}
		else {y <- y; N <- N}
		numVar <- 1
		numIDV <- length(y)
		if (is.null(offset)) {offset <- rep(0, numIDV)}
		if (numCov == 0){
			model0 <- glm(formula = y/N ~ 1 + offset(offset), family = binomial(link = "logit"), weights = N)
			idx    <- match(rownames(model.frame(formula = y/N ~ 1 + offset(offset), na.action = na.omit)), 
					rownames(model.frame(formula = y/N ~ 1 + offset(offset), na.action = na.pass)))
			}
		else {
			model0 <- glm(formula = y/N ~ X + offset(offset), family = binomial(link = "logit"), weights = N)
			idx    <- match(rownames(model.frame(formula = y/N ~ X + offset(offset), na.action = na.omit)),
					rownames(model.frame(formula = y/N ~ X + offset(offset), na.action = na.pass)))
			}
		model0$numN <- N[idx];
		model0$numS <- y[idx]
		for (ik in 1:(length(kmat))) {kmat[[ik]] <- kmat[[ik]][idx, idx]}
		if (method == "AI") {
		model1 <- try(glmmPQL_fit(model0, kmat, method = "AI", maxiter, tol, verbose))}
		else {
		if (method == "EI") {
		model1 <- try(glmmPQL_fit(model0, kmat, method = "EI", maxiter, tol, verbose))}
		else {
		if (method == "OI") {
		model1 <- try(glmmPQL_fit(model0, kmat, method = "OI", maxiter, tol, verbose))}}
			}
		if (class(model1) != "try-error"&!is.null(model1)){
			if (verbose){cat(paste("glmmPQL::BMM::tau = ", model1$theta,"\n"))}
				}
		else {converged <- FALSE}
		closeAllConnections()
		return(model1)
	}# end BMM	
}# end function glmmPQL



##########################################################
#           	   glmmPQL FIT FUNCTION			 #
##########################################################
glmmPQL_fit <- function(
	model0,
	kmat,
	method = "AI",
	maxiter,
	tol,
	verbose) {

	#names(kmat) <- paste("kins", 1:length(kmat), sep="")
	if (method == "AI") {
		fixtau_old 	<- rep(0, length(kmat)+1)
		model1 		<- glmmPQL_optim(model0, kmat, method = "AI", 
						maxiter = maxiter, tol = tol, verbose = verbose)
		fixtau_new 	<- 1*(model1$theta < 1.01 * tol)
		while (any(fixtau_new != fixtau_old)) {
			fixtau_old <- fixtau_new
			model1 	<- glmmPQL_optim(model0, kmat, method = "AI", 
						 fixtau = fixtau_old, maxiter = maxiter, tol = tol, verbose = verbose)
			fixtau_new <- 1*(model1$theta < 1.01 * tol)
		}
	}

	else {
	if (method == "EI") {
		fixtau_old 	<- rep(0, length(kmat)+1)
		model1 		<- glmmPQL_optim(model0, kmat, method = "EI", 
						 maxiter = maxiter, tol = tol, verbose = verbose)
		fixtau_new 	<- 1*(model1$theta < 1.01 * tol)
		while(any(fixtau_new != fixtau_old)) {
			fixtau_old <- fixtau_new
			model1 	<- glmmPQL_optim(model0, kmat, method = "EI", 
						 fixtau = fixtau_old, maxiter = maxiter, tol = tol, verbose = verbose)
			fixtau_new <- 1*(model1$theta < 1.01 * tol)
		}
	}

	else {
	if (method == "OI") {
		fixtau_old 	<- rep(0, length(kmat)+1)
		model1 		<- glmmPQL_optim(model0, kmat, method = "OI",
						 maxiter = maxiter, tol = tol, verbose = verbose)
		fixtau_new 	<- 1*(model1$theta < 1.01 * tol)
		while (any(fixtau_new != fixtau_old)) {
			fixtau_old <- fixtau_new
			model1 	<- glmmPQL_optim(model0, kmat, method = "OI",
						 fixtau = fixtau_old, maxiter = maxiter, tol = tol, verbose = verbose)
			fixtau_new <- 1*(model1$theta < 1.01 * tol)
		}
	}
	#else {model1 <- NULL}
	}
	}
	return(model1)
}

##########################################################
#       glmmPQL FIT AVERAGE INFORMATION FUNCTION	 #
##########################################################
glmmPQL_optim <- function(
	model0,
	kmat,
	method = "AI",
	tau = rep(0, length(kmat)+1),
	fixtau = rep(0, length(kmat)+1),
	maxiter,
	tol,
	verbose) {

	if (model0$family$family %in% c("binomial")) { y <- model0$numS}
	else {y <- model0$y}
	numIDV <- length(y)
	offset <- model0$offset
	
	family <- model0$family
	eta    <- model0$linear.predictors# Xb
	mu     <- model0$fitted.values# exp(Xb)/(1+exp(Xb))
	mu.eta <- family$mu.eta(eta)  # mu*(1-mu)
	D <- mu.eta/sqrt(model0$family$variance(mu)) #sqrtW or sqrtD: sqrt(mu.eta)
	#sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0))*family$variance(mu))
	#print(cbind(D, sqrt(mu*(1-mu)))[1:10,]); D = sqrt(mu*(1-mu))
	if (family$family %in% c("binomial")){
	  mu.eta <- model0$numN*mu.eta
	  D <- mu.eta/sqrt(model0$numN*model0$family$variance(mu))
	  mu <- model0$numN*mu
	}

	Y <- eta - offset + (y - mu)/mu.eta	
	X <- model.matrix(model0)
	alpha <- model0$coef
	
	if (family$family %in% c("poisson", "binomial")) {
		tau[1] <- 1
		fixtau[1] <- 1
	}

	numK   <- length(kmat)
	idxtau <- which(fixtau == 0)
	numK2  <- sum(fixtau == 0)

	if (numK2 > 0) {
		tau[fixtau == 0] <- rep(min(0.9,var(Y)/(numK+1)), numK2)
		H <- tau[1]*diag(1/D^2)
		for (ik in 1:numK) {H <- H + tau[ik+1]*kmat[[ik]]}	
		Hinv 	<- chol2inv(chol(H))
		HinvX 	<- crossprod(Hinv, X)
		XHinvX 	<- crossprod(X, HinvX)
		P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))	
		if (class(P) == "try-error"){
			stop("Error in P matrix calculation!")
		}

		PY <- crossprod(P, Y)
		tau0 <- tau
		for (ik in 1:numK2) {
			if (ik == 1 && fixtau[1] == 0) 
			{tau[1] <- max(0, tau0[1] + tau0[1]^2 * (sum((PY/D)^2) - sum(diag(P)/D^2))/numIDV)}
			else {
				PKPY <- crossprod(P, crossprod(kmat[[idxtau[ik]-1]], PY))
				tau[idxtau[ik]] <- max(0, tau0[idxtau[ik]] + 
				tau0[idxtau[ik]]^2 * (crossprod(Y, PKPY) - sum(P*kmat[[idxtau[ik]-1]]))/numIDV)
			}
		}
	}
	
	for (iter in 1:(maxiter)) {
		if (verbose) {print (iter)}
		alpha0 	<- alpha
		tau0 	<- tau
		if (method == "AI") {model1 <- glmmAI(Y, X, length(kmat), kmat, D^2, tau, fixtau, tol)}
		else {
		if (method == "OI") {model1 <- glmmOI(Y, X, length(kmat), kmat, D^2, tau, fixtau, tol)}
		else {
		if (method == "EI") {model1 <- glmmEI(Y, X, length(kmat), kmat, D^2, tau, fixtau, tol)}}
		}

		tau    <- as.numeric(model1$tau)
		cov    <- as.matrix(model1$cov)
		alpha  <- as.numeric(model1$alpha)
		eta    <- as.numeric(model1$eta) + offset
		mu     <- family$linkinv(eta)
		mu.eta <- family$mu.eta(eta)
		D      <- mu.eta/sqrt(family$variance(mu))
		if (family$family %in% c("binomial")){
		mu.eta <- model0$numN*mu.eta
		D      <- mu.eta/sqrt(model0$numN*family$variance(mu))
		mu     <- model0$numN*mu
		}
		Y      <- eta - offset + (y - mu)/mu.eta
		if (2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol),
			abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) {break}
		if (max(tau) > tol^(-2) | any(is.infinite(D)) | any(is.infinite(mu)) | any(is.infinite(eta))) {
			iter <- maxiter
			break
		}
	}

	converged <- ifelse(iter < maxiter, TRUE, FALSE)
	res <- y - mu
	P <- model1$P
	H <- model1$H
	iter <- iter
	FisherMat <- model1$FisherMat

	return(list(
	theta = tau,#第一个表示残差；第二个表示遗传方差
	FisherMat = FisherMat,
	coefficients = alpha,
	cov = cov,
	linear.predictors = eta,
	fitted.values = mu,
	Y = Y,
	P = P,
	H = H,#sigma
	#u = u,
	residuals = res,
	iter = iter,
	converged = converged))
}# end function
#########################################
#             CODE END                  #
#########################################
