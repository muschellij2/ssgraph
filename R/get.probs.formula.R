
get.probs.formula <- function(graph, x, y, family=NULL, weight, ngroups, levs, nkeep, nsubj, corr=FALSE, nadj, usepvals=FALSE, pvals=NULL){
	# if is.null(family)
	
	if (usepvals & is.null(pvals)) stop("Using p-values for probs, but none given")
	if (usepvals & (length(pvals) != nkeep | length(pvals) != ncol(graph))) stop("Something off about p-values")
 	eta <- 1/(10*nsubj)

	if (is.null(family) & nadj > 2) family <- list(family="multinom")
	if (is.null(family) & nadj == 2) family <- binomial()
	#if (!weight & is.null(x)) family <- binomial()
	#if (weight & is.null(x)) family <- gaussian()
	# if (weight & is.null(x) & corr=TRUE) family <- beta()
		
	strfam <- as.character(family$family)
	#print(strfam)
	#print(family)
	if (strfam == "gaussian") dens <- "norm"
	if (strfam %in% c("poisson", "quasipoisson")) dens <- "pois"
	if (strfam == "binomial") dens <- "binom"
	if (strfam == "Gamma") dens <- "gamma"
	if (strfam == "inverse.gaussian") stop("IG not implemented yet")	
	if (strfam == "multinom") dens <- "multinom"
	if (strfam == "beta") dens <- "beta"
	
	# colnames(all.probs) <- levs
	all.probs <- all.params <- all.mods <- Ghat <- runs <- vector(length=ngroups, mode="list")

ig <- 1
icol <- 353
	pi_y <- sapply(1:ngroups, function(x) return(mean(y == levs[x], na.rm=TRUE)))
	if (sum(pi_y) != 1) stop("Something off about pi_y")
	### looping over groups
	for (ig in 1:ngroups){
		
		## grab the level of ylabels
		ilev <- levs[ig]
		ingroup <- y == ilev
		### may not need line below - may actually be bad
		ingroup[is.na(ingroup)] <- FALSE
		
		### get each column probability under certain model
		col.probs <- vector(mode="list", length=ncol(graph))
		# G <- G[ingroup, ]
		# z <- x[ingroup, ]	
		# run.data <- data.frame(z=z)
		# run.data$G <- G					
		# col.probs <- apply(graph, 2, function(x) glm(x ~ ., family=family))
		col.probs <- apply(graph, 2, function(x) glm(x ~ 1, family=family))
		
		### loop over the signal subgraph
		for (icol in 1:ncol(graph)){
			glevs <- NA
			G <- graph[, icol]
			if (is.factor(G)) glevs <- levels(G)
			z <- x[ingroup]				

			### Need to make it a factor 
			if (dens == "multinom") {
				if (!is.factor(G)) G <- factor(G)
					glevs <- levels(G)					
					G <- G[ingroup]
					z <- x[ingroup]
					run.data <- data.frame(z=z)
					run.data$G <- G
					### if all have same value - ASK BRIAN - like if it was poisson and all 3's then how can we say that like 3 is probability 1 and the others
					if (length(unique(G)) == 1){
						mod <- glm(G ~ ., family=binomial(), data=run.data)
					} else {
						mod <- multinom(G ~ ., data=run.data, trace=FALSE)
					}
			} else {
				G <- G[ingroup]
				# run.data <- data.frame(cbind(G=G, z))
				mod <- glm(G ~ z, family=family)
			}
			# print(icol)
			all.mods[[ig]][[icol]] <- mod
		}	 ## over icol	
	
		# all.preds <- pred.mods(all.mods, ig, ngroups, graph, glevs, icol)
		
			### Get Ghat and get parameter estimates from it
			### this is going to be an length(ingroup) by ngroup vector
			if (dens != "multinom") {
				if (inherits(mod, 'glm')){
					# Ghat <- predict(mod, newdata=run.data, type="response")
					# Ghat[[ig]] <- lapply(all.mods[[ig]], function(mod) predict(mod, newdata=run.data, type="response"))
					Ghat[[ig]] <- lapply(all.mods[[ig]], function(mod) predict(mod, type="response"))
				} else {
					# Ghat <- predict(mod, newdata=run.data, type="probs")
	# 				Ghat[[ig]] <- lapply(all.mods[[ig]], function(mod) predict(mod, newdata=run.data, type="probs"))
					Ghat[[ig]] <- lapply(all.mods[[ig]], function(mod) predict(mod, newdata=run.data, type="probs"))
				}
			}
			### need to do something different multinom and beta
			if (dens == "multinom") {
				stop("Multinom hasn't been done")
				if (inherits(mod, 'glm')){
					Ghat <- predict(mod, newdata=run.data, type="response")
					Ghat <- matrix(Ghat, ncol=length(unique(G)) )
#					colnames(Ghat) <- levels(G)					
				} else {
					Ghat <- predict(mod, newdata=run.data, type="probs")
				}
			}
			print(icol)
	}
	
	
			
			for (ig in 1:ngroups){
				for (icol in 1:ncol(graph)) {
					G <- graph[, icol]
					gGhat <- Ghat[[ig]][[icol]]
					runs[[ig]][[icol]] <- dens.probs(G, gGhat, dens, glevs)
				}
			}
			
			all.params <- lapply(runs, function(z) lapply(z, function(x)return(x$params)))
			all.probs <- lapply(runs, function(z) sapply(z, function(x) return(x$probs)))
			
			preds <- lapply(all.probs, function(prob) apply(prob, 1, function(z) sum(log(z))))
			preds <- lapply(1:ngroups, function(ig) preds[[ig]] + log(pi_y[ig]))
			### can we work on log scale for sum? I don't think so
			preds <- do.call(cbind, preds)
			epreds <- exp(preds)
			### Normalize to add to 1 
			preds <- epreds/rowSums(epreds)
			# print(run)
			# col.probs[[icol]] <- run$probs
			# all.params[[ig]][icol] <- run$params

		# #print(paste("Dim col.probs", dim(col.probs)))
		# #print(paste("Length pvals", length(pvals)))
		# print(head(col.probs))	
		# if (usepvals) {
			# ## use weighted sum on the log scale
			# for (iicol in 1:length(col.probs)) col.probs[[iicol]] <- log(col.probs[[iicol]])*pvals[iicol]
			# # col.probs <- sapply(col.probs, 1, function(edge) exp(sum(edge)))
			# col.probs <- sapply(col.probs, function(edge) exp(sum(edge)))
		# } else {
			# # col.probs <- apply(col.probs, 1, prod)
			# col.probs <- sapply(col.probs, prod)
		# }
# #		print(head(col.probs), 20)
# #		stop("here")
		# col.probs <- ifelse(1-col.probs < 1e-9, 1-eta, ifelse(col.probs < 1e-9, eta, col.probs))
		# all.probs[, ig] <- col.probs
	# }
	# names(all.mods) <- levs
	
	# colnames(all.probs) <-levs
	return(list(probs=preds, mods = all.mods, dens=dens, priors=pi_y))
}




#### Dens.probs - find the probabilities from respective densities
dens.probs <- function(G, Ghat, dens, glevs, z = NULL){
		if (dens == "multinom") {
			
			#### if it's binary - then make it 2 columns
			if (length(glevs) == 2) Ghat <- cbind(Ghat, 1-Ghat)
			
			### Need to make it ncol = nsubj, and nrow = number of levels for dmultinom
			tmp <- matrix(rep(0, length=length(glevs)*length(G)), ncol=length(G), nrow=length(glevs), byrow=TRUE)
			G2 <- as.numeric(G)
			for (iicol in 1:ncol(tmp)) tmp[G2[iicol], iicol] <- 1
			G <- tmp
			rownames(G) <- glevs
			mprobs <- colMeans(Ghat)
			params <- mprobs
			probs <- matrix(apply(G, 2, function(edge) dmultinom(edge, size=1, prob=mprobs)), ncol=length(G), nrow=length(glevs))
			rownames(probs) <- glevs
		}			
		
		if (dens == "pois") {
			lambda <- mean(Ghat)
			params <- lambda
			probs <- dpois(G, lambda=lambda)
		}
		if (dens == "binom") {
			phat <- mean(Ghat)
			params <- phat
			# probs <- Ghat
			probs <- dbinom(G, size=1,  prob=phat)
		}
		if (dens == "norm") {
			m <- mean(Ghat)
			sds <- sd(Ghat)
			params <- c(m, sds)
			probs <- dnorm(G, m, sd=sds)
		}
		if (dens == "gamma") {
			### need to scale by 1000  - (https://stat.ethz.ch/pipermail/r-help/2011-August/286125.html)
			tester <- rep(Ghat[1], length=length(Ghat))
			if ( all(Ghat == tester) ) {
				stop("Gamma parameters cannot be estimated - all predictions are equal")
			}
			ests <- fitdistr(Ghat/1000, "gamma")
			ests <- ests$estimate
			params <- ests
			probs <- dgamma(G, shape = ests["shape"], rate = ests["rate"]/1000)
		}
		###need to implement for beta
		if (dens %in% c("beta")) {
			stop("Not Implemented yet for this family")
		}
		return(list(probs=probs, params=params))
}


pred.mods <- function(all.mods, ig, ngroups, graph, glevs, icol){
			### Get Ghat and get parameter estimates from it
			### this is going to be an length(ingroup) by ngroup vector
			if (dens != "multinom") {
				if (inherits(mod, 'glm')){
					# Ghat <- predict(mod, newdata=run.data, type="response")
					# Ghat[[ig]] <- lapply(all.mods[[ig]], function(mod) predict(mod, newdata=run.data, type="response"))
					Ghat[[ig]] <- lapply(all.mods[[ig]], function(mod) predict(mod, type="response"))
				} else {
					# Ghat <- predict(mod, newdata=run.data, type="probs")
	# 				Ghat[[ig]] <- lapply(all.mods[[ig]], function(mod) predict(mod, newdata=run.data, type="probs"))
					Ghat[[ig]] <- lapply(all.mods[[ig]], function(mod) predict(mod, newdata=run.data, type="probs"))
				}
			}
			### need to do something different multinom and beta
			if (dens == "multinom") {
				stop("Multinom hasn't been done")
				if (inherits(mod, 'glm')){
					Ghat <- predict(mod, newdata=run.data, type="response")
					Ghat <- matrix(Ghat, ncol=length(unique(G)) )
#					colnames(Ghat) <- levels(G)					
				} else {
					Ghat <- predict(mod, newdata=run.data, type="probs")
				}
			}
			print(icol)		
				
	return(list(probs=preds, mods = all.mods, dens=dens, priors=pi_y))
}