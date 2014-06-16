getPckg <- function(pckg) install.packages(pckg, repos = "http://cran.r-project.org")
checkPckg <- function(pckg){
	pckg = eval(parse(text=paste("try(require(", pckg, "))")))
	if(!pckg) {
		cat(paste("Installing '", pckg, "' from CRAN\n"))
		getPckg(pckg)
		eval(parse(text=paste("require(", pckg, ")")))
	}
}

# function(wps) -log(wps)

family.sig.sgraph <- function(object, ...) object$family
#binary signal subgraph
sig.sgraph <- sig.sgraph.default <- function(graph, y=NULL, x=NULL, corr=FALSE, family=NULL, na.action = NULL, upper=FALSE, 
	alternative = "two.sided", keep=0.1, weight=FALSE, workspace=400000, formula=NULL, data=NULL,  ...){

	checkPckg("MASS")
	checkPckg("gamlss.dist")
	checkPckg("nnet")
	print(alternative)
	## get levels of adjacencies
	adj <- sort(unique(c(graph)))
	nadj <- length(unique(c(graph)))
	#print("Code block 1")
	
	if (nadj > 10 & !weight) {
		cat("Over 10 labels - running weighted version?\n")
		weight <- TRUE
	}
	
	if (nadj == 2 & weight) {
		cat("Only two labels - using binary classifier\n")
		weight <- FALSE
		family <- binomial()
	}
	
	if ((min(graph) < -1 | max(graph) > 1) && corr) {
			corr <- FALSE
			cat("corr specified as TRUE, but outside {-1, 1}, not transforming (corr+1)/2\n")
	}	

	####Using weighted version of classifier or not
	if (is.null(family) & !corr){
		cat("No Family Specified, assuming Multinomial\n")
		family <- list(family="multinom")
	}
	
	if (weight) {
		## getting the family used for the graph model (multinom is used for P(Y|X))
		if (is.character(family)) {
			### make it so that family$family is the true family
			if (!(family %in% c("multinom", "beta")))	{
				family <- get(family, mode = "function", envir = parent.frame())
			} else {
			   	family <- list(family=family)
	    	}
		}
	    if (is.function(family)) family <- family()
	    if (is.null(family$family) & corr) {
	    		cat("No Family Specified and Corr=TRUE, binomial chosen")
	    		family <- binomial()
		}
	    if (is.null(family$family)) {
	        # print(family)
	        stop(paste("'family'", "family", "not recognized"))
	    }
    }
    
	#print("IN SGRAPH")
	#print(head(y))
    print(family)
	## taken from fisher.test <- check if alternative is in acceptable range
	alternative <- char.expand(alternative, c("two.sided", "less", "greater"))
	if (length(alternative) > 1L || is.na(alternative)) stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
	
	## get dimensions of the graph and the class - matrix/data.frame/array
	dag <- dim(graph)
	cg <- class(graph)
	## see if data is stored as array
	## if ylabels are a column in data.frame, then pull it out
	if (cg == "data.frame" & ("y" %in% colnames(graph) & is.null(y) )) {
		y <- graph$y
		which.y <- which(colnames(graph) == "y")
		graph <- graph[,-which.y]
	}
	if (cg == "matrix" && ("y" %in% colnames(graph) & is.null(y) )) {
		y <- graph[,"y"]
		which.y <- which(colnames(graph) == "y")
		graph <- graph[,-which.y]
	}
	
	# convert into correct format - rows are units/subjects, columns are nodes
	graph <- convert.graph(graph, upper=upper)
	
	nullx <- FALSE
	form <- !is.null(formula)

	if (is.null(y)) stop("No Y input")
	if (is.null(x)) {
		x <- matrix(rep(1, length(y)), ncol=1)
		nullx <- TRUE
		## if formula is specified - should always have something here, even if just intercept
		if (form) stop("No X specified or only intercept")
	}
	
	## check to see Number of graphs equals N
	if (any(is.na(x))) 
        stop("NA not permitted in predictors")
    if (any(is.na(graph))) 
        stop("NA not permitted in graph")
    if (any(is.na(y))) 
        stop("NA not permitted in response")
    if (is.null(dim(x))) 
            dim(x) <- c(length(x), 1)
    if (nrow(graph) != nrow(x) | nrow(x) != length(y)){
		stop("Number of graphs does not match N Subjects from X or Y")
	}
		
	## y is outcome 
	#### change if this fails for dim on a vector
	dy <- dim(y)
	cy <- class(y)
	if (cy == "matrix")
		y <- c(y)
	if (!(cy %in% c("matrix", "numeric", "character", "factor", "integer", "logical")))
		stop("Class of Y unrecognized/too big - array")	
		
	## Get number of groups
	ngroups <- length(unique(y))
	if (ngroups < 2) stop("Need at least 2 groups")
	if (alternative != "two.sided" & ngroups > 2) {
		cat("For groups/labels > 2, alternative is two sided\n")
		alternative <- "two.sided"
	}
	
	if (class(y) != "factor") y <- factor(y)
	# print("Code block 2")
	levs <- levels(y)
	if (ngroups == 2) {
		## make 0/1 variable
		if (class(levs) == "character") y <- as.numeric(factor(y, levels=levs))-1
		if (class(levs) == "factor") y <- as.numeric(y, levels=levs)-1
		if (!all(as.numeric(levs)==c(0, 1))) stop("Label problems")
	}
	#print("Head X")
	#print(head(x))
	# if (form) {
		# mod.priors <- multinom(formula=formula, data=data, trace=FALSE)
	# } else mod.priors <- multinom(formula=y ~ 1, trace=FALSE)
	
	# priors <- predict(mod.priors, type="prob")
	#print("Priors")
	#print(head(priors))
	
	
	dg <- dim(graph)
	## get number of subjects/vertices
	nvert <- dg[2]
	nsubj <- dg[1]
	## rows are units/subjects, columns are vertices
	if (length(dg) > 2) stop("Problem with dimensions")
	
	## get p.value for the graph vertices
	print("in getpvals")
	
	## try rowFtests
	pvals <- apply(graph, 2, function(edge) getpvals(edge, y=y, workspace=workspace, weight=weight, alternative=alternative))
	
	#print(head(sort(pvals)))
	#print("Code block 3")
	
	## keep percentage of the graph (need to CV it probably) or a function to weight with p-values
	# normalize it
	if (is.function(keep)) {
		nkeep <- nvert
		pvals <- keep(pvals)/sum(keep(pvals))
	} else nkeep <- floor(nvert*keep)
	
	## ssgraph: indices of the signal subgraph
	ssgraph <- graph.ind <- sort(order(pvals)[1:nkeep])
	
	if (cg == "array") {
		sig.graph <- matrix(0, nrow=dag[1], ncol=dag[2])
		## make sure to put the ssgraph back in the right spot if upper is on
		if (upper){
			up.tri <- upper.tri(sig.graph, diag=TRUE)
			sig.graph[up.tri][ssgraph] <- 1
		} else sig.graph[ssgraph] <- 1
	} else {
		sig.graph <- rep(0, nvert)
		sig.graph[ssgraph] <- 1
	}
	
#	print("Code block 4")
	
	## only need product over signal subgraph
	## tgraph is total graph
	tgraph <- graph
	graph <- graph[, ssgraph]

	if (corr) {
		graph <- (graph+1)/2
	}

	# print("Code block 5")
	#print("HEAD X IT")
	#print(nkeep)
	print("at get.probs.formula")
	print(list(family=family, weight=weight, ngroups=ngroups, levs=levs, nkeep=nkeep, nsubj= nsubj, corr=corr, nadj= nadj))
	res <- get.probs.formula(graph=graph, x=x, y=y, family=family, weight=weight, ngroups=ngroups, levs=levs, nkeep=nkeep, nsubj= nsubj, corr=corr, nadj= nadj, usepvals=is.function(keep), pvals=pvals)
	
	probs <- res$probs
	priors <- res$priors
	### probs are the probability
	# probs <- res$probs
	# if ( length( unique(probs) ) == 1 ) {
		# # print(head(probs))
		# # print(unique(probs))
		# stop("P(G| Y, X) are all the same")
	# }
	# #print("Head Probs")
	# #print(head(probs))
	# for (icol in 1:ncol(probs)) {
		# if (colnames(probs)[icol] != colnames(priors)[icol]) stop("Something wrong with predictions")	
		# probs[, icol] <- probs[, icol]*priors[, icol]
	# }
	#priors <- res$priors
	
	#print("Code block 6")
	## divide by total to get actual probability and not just numerator P(G = g | Y = y)
	## apply(probs, 1, sum) = P(G=g|Y=1)P(Y=1)+ ... + P(G=g|Y=y_n)P(Y=y_n)
	# probs <- probs / apply(probs, 1, sum)
	#print("Head Probs")
	#print(head(probs))
	
	colnames(probs) <- levs
	preds <- apply(probs, 1, function(x) which(x == max(x)))
	if (class(preds) != "integer") stop("Problem with Prediction")
	preds <- factor(colnames(probs)[preds], levels=levs)
	#preds <- factor(colnames(probs)[preds], levels=levs)
	#preds <- levs[preds]
	# print(length(preds))
	# print(head(preds))
	print(table(preds))
	print(table(preds, y))
	accuracy <- mean(preds == y)
	cat(sprintf("Accuracy is %3.2f%% \n", accuracy*100))
	
	xlevels <- apply(x, 2, levels)
	
	# mod.priors=mod.priors,
	ret <- list(ssgraph=sig.graph, priors=priors, prob=probs, predicted=preds, y=y, accuracy=accuracy, pvals=pvals, usepvals=is.function(keep), graph.models = res$mods,  graph.class = cg, graph.indices = graph.ind, formula=formula, has.formula=form, xlevels=xlevels, terms=NULL, x=x, family=family, dens=res$dens, upper=upper, dg = dg)
    class(ret) <- c("sig.sgraph", class(ret))
	return(ret)
}




## Rpart 
#fit <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis)
# unique(rownames(fit$splits))

#randomForest
#iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE, proximity=TRUE)
# iris.rf$importance


#Possible function for formula
# sig.sgraph.formula <- function (formula, data, subset, na.action, ...) 
# {
    # if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
        # "term.labels")) != 1L)) 
        # stop("'formula' missing or incorrect")
    # m <- match.call(expand.dots = FALSE)
    # if (is.matrix(eval(m$data, parent.frame()))) 
        # m$data <- as.data.frame(data)
    # m[[1L]] <- as.name("model.frame")
    # m$... <- NULL
    # mf <- eval(m, parent.frame())
    # DNAME <- paste(names(mf), collapse = " by ")
    # names(mf) <- NULL
    # response <- attr(attr(mf, "terms"), "response")
    # g <- factor(mf[[-response]])
    # if (nlevels(g) != 2L) 
        # stop("grouping factor must have exactly 2 levels")
    # DATA <- split(mf[[response]], g)
    # names(DATA) <- c("x", "y")
    # y <- do.call("t.test", c(DATA, list(...)))
    # y$data.name <- DNAME
    # if (length(y$estimate) == 2L) 
        # names(y$estimate) <- paste("mean in group", levels(g))
    # y
# }


# function (formula, data = NULL, ..., subset, na.action = na.fail) 
# {
    # if (!inherits(formula, "formula")) 
        # stop("method is only for formula objects")
    # m <- match.call(expand = FALSE)
    # if (any(c("xtest", "ytest") %in% names(m))) 
        # stop("xtest/ytest not supported through the formula interface")
    # names(m)[2] <- "formula"
    # if (is.matrix(eval(m$data, parent.frame()))) 
        # m$data <- as.data.frame(data)
    # m$... <- NULL
    # m$na.action <- na.action
    # m[[1]] <- as.name("model.frame")
    # m <- eval(m, parent.frame())
    # y <- model.response(m)
    # Terms <- attr(m, "terms")
    # attr(Terms, "intercept") <- 0
    # m <- model.frame(terms(reformulate(attributes(Terms)$term.labels)), 
        # data.frame(m))
    # for (i in seq(along = ncol(m))) {
        # if (is.ordered(m[[i]])) 
            # m[[i]] <- as.numeric(m[[i]])
    # }
    # ret <- randomForest(m, y, ...)
    # cl <- match.call()
    # cl[[1]] <- as.name("randomForest")
    # ret$call <- cl
    # ret$terms <- Terms
    # if (!is.null(attr(m, "na.action"))) 
        # ret$na.action <- attr(m, "na.action")
    # class(ret) <- c("randomForest.formula", "randomForest")
    # return(ret)
# }


# function (formula, data = NULL, ..., subset, na.action = na.omit, 
    # scale = TRUE) 
# {
    # call <- match.call()
    # if (!inherits(formula, "formula")) 
        # stop("method is only for formula objects")
    # m <- match.call(expand.dots = FALSE)
    # if (identical(class(eval.parent(m$data)), "matrix")) 
        # m$data <- as.data.frame(eval.parent(m$data))
    # m$... <- NULL
    # m$scale <- NULL
    # m[[1]] <- as.name("model.frame")
    # m$na.action <- na.action
    # m <- eval(m, parent.frame())
    # Terms <- attr(m, "terms")
    # attr(Terms, "intercept") <- 0
    # x <- model.matrix(Terms, m)
    # y <- model.extract(m, "response")
    # attr(x, "na.action") <- attr(y, "na.action") <- attr(m, "na.action")
    # if (length(scale) == 1) 
        # scale <- rep(scale, ncol(x))
    # if (any(scale)) {
        # remove <- unique(c(which(labels(Terms) %in% names(attr(x, 
            # "contrasts"))), which(!scale)))
        # scale <- !attr(x, "assign") %in% remove
    # }
    # ret <- svm.default(x, y, scale = scale, ..., na.action = na.action)
    # ret$call <- call
    # ret$call[[1]] <- as.name("svm")
    # ret$terms <- Terms
    # if (!is.null(attr(m, "na.action"))) 
        # ret$na.action <- attr(m, "na.action")
    # class(ret) <- c("svm.formula", class(ret))
    # return(ret)
# }


# e1071:::predict.svm
# function (object, newdata, decision.values = FALSE, probability = FALSE, 
    # ..., na.action = na.omit) 
# {
    # if (missing(newdata)) 
        # return(fitted(object))
    # if (object$tot.nSV < 1) 
        # stop("Model is empty!")
    # if (inherits(newdata, "Matrix")) {
        # library("SparseM")
        # library("Matrix")
        # newdata <- as(newdata, "matrix.csr")
    # }
    # if (inherits(newdata, "simple_triplet_matrix")) {
        # library("SparseM")
        # ind <- order(newdata$i, newdata$j)
        # newdata <- new("matrix.csr", ra = newdata$v[ind], ja = newdata$j[ind], 
            # ia = as.integer(cumsum(c(1, tabulate(newdata$i[ind])))), 
            # dimension = c(newdata$nrow, newdata$ncol))
    # }
    # sparse <- inherits(newdata, "matrix.csr")
    # if (object$sparse || sparse) 
        # library("SparseM")
    # act <- NULL
    # if ((is.vector(newdata) && is.atomic(newdata))) 
        # newdata <- t(t(newdata))
    # if (sparse) 
        # newdata <- SparseM::t(SparseM::t(newdata))
    # preprocessed <- !is.null(attr(newdata, "na.action"))
    # rowns <- if (!is.null(rownames(newdata))) 
        # rownames(newdata)
    # else 1:nrow(newdata)
    # if (!object$sparse) {
        # if (inherits(object, "svm.formula")) {
            # if (is.null(colnames(newdata))) 
                # colnames(newdata) <- colnames(object$SV)
            # newdata <- na.action(newdata)
            # act <- attr(newdata, "na.action")
            # newdata <- model.matrix(delete.response(terms(object)), 
                # as.data.frame(newdata))
        # }
        # else {
            # newdata <- na.action(as.matrix(newdata))
            # act <- attr(newdata, "na.action")
        # }
    # }
    # if (!is.null(act) && !preprocessed) 
        # rowns <- rowns[-act]
    # if (any(object$scaled)) 
        # newdata[, object$scaled] <- scale(newdata[, object$scaled, 
            # drop = FALSE], center = object$x.scale$"scaled:center", 
            # scale = object$x.scale$"scaled:scale")
    # if (ncol(object$SV) != ncol(newdata)) 
        # stop("test data does not match model !")
    # ret <- .C("svmpredict", as.integer(decision.values), as.integer(probability), 
        # as.double(if (object$sparse) object$SV@ra else t(object$SV)), 
        # as.integer(nrow(object$SV)), as.integer(ncol(object$SV)), 
        # as.integer(if (object$sparse) object$SV@ia else 0), as.integer(if (object$sparse) object$SV@ja else 0), 
        # as.double(as.vector(object$coefs)), as.double(object$rho), 
        # as.integer(object$compprob), as.double(object$probA), 
        # as.double(object$probB), as.integer(object$nclasses), 
        # as.integer(object$tot.nSV), as.integer(object$labels), 
        # as.integer(object$nSV), as.integer(object$sparse), as.integer(object$type), 
        # as.integer(object$kernel), as.integer(object$degree), 
        # as.double(object$gamma), as.double(object$coef0), as.double(if (sparse) newdata@ra else t(newdata)), 
        # as.integer(nrow(newdata)), as.integer(if (sparse) newdata@ia else 0), 
        # as.integer(if (sparse) newdata@ja else 0), as.integer(sparse), 
        # ret = double(nrow(newdata)), dec = double(nrow(newdata) * 
            # object$nclasses * (object$nclasses - 1)/2), prob = double(nrow(newdata) * 
            # object$nclasses), PACKAGE = "e1071")
    # ret2 <- if (is.character(object$levels)) 
        # factor(object$levels[ret$ret], levels = object$levels)
    # else if (object$type == 2) 
        # ret$ret == 1
    # else if (any(object$scaled) && !is.null(object$y.scale)) 
        # ret$ret * object$y.scale$"scaled:scale" + object$y.scale$"scaled:center"
    # else ret$ret
    # names(ret2) <- rowns
    # ret2 <- napredict(act, ret2)
    # if (decision.values) {
        # colns = c()
        # for (i in 1:(object$nclasses - 1)) for (j in (i + 1):object$nclasses) colns <- c(colns, 
            # paste(object$levels[object$labels[i]], "/", object$levels[object$labels[j]], 
                # sep = ""))
        # attr(ret2, "decision.values") <- napredict(act, matrix(ret$dec, 
            # nrow = nrow(newdata), byrow = TRUE, dimnames = list(rowns, 
                # colns)))
    # }
    # if (probability && object$type < 2) 
        # attr(ret2, "probabilities") <- napredict(act, matrix(ret$prob, 
            # nrow = nrow(newdata), byrow = TRUE, dimnames = list(rowns, 
                # object$levels[object$labels])))
    # ret2
# }
# <environment: namespace:e1071>
