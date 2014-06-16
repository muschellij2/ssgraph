predict.sig.sgraph <- function (object, newdata, newgraph, type = "response", ...) 
{

    if (!inherits(object, "sig.sgraph")) 
        stop("object not of class sig.sgraph")

	# graph <- newgraph
#	Need to check for na 
	if (sum(is.na(newgraph) > 0)){
		stop("NA's in newgraph, may impute data for sensitivity analysis")
	}
	
	dag <- dim(newgraph)
	cg <- class(newgraph)
	if (object$graph.class != cg){
		warning("Classes for new newgraph different than original")
	}

	
	if (cg == "data.frame" & ("y" %in% colnames(newgraph)) ) {
		which.y <- which(colnames(newgraph) == "y")
		newgraph <- newgraph[,-which.y]
	}
	if (cg == "matrix" && ("y" %in% colnames(newgraph) ) ) {
		which.y <- which(colnames(newgraph) == "y")
		newgraph <- newgraph[,-which.y]
	}
	
	tnewgraph <- newgraph
	newgraph <- newgraph[, object$ssgraph]
		
	# convert into correct format - rows are units/subjects, columns are nodes
	newgraph <- convert.graph(newgraph, upper=object$upper)
	dg <- dim(newgraph)
	nvert <- dg[2]
	print(dg)
	print(object$dg)
	if (nvert != object$dg[2]) 
		stop("Number of vertices not the same as in training data")
	
		
	
	if (object$usepvals & is.null(object$pvals)) stop("Using p-values for probs, but none given")
	# if (object$usepvals & (length(object$pvals) != nkeep | length(object$pvals) != ncol(object$graph))) stop("Something off about p-values")
 	# eta <- 1/(10*nsubj)
	
    out.type <- charmatch(tolower(type), c("response", "prob", "class"))
    if (is.na(out.type)) 
        stop("type must be one of 'response', 'prob'")
    if (out.type == 3) 
        out.type <- 1
    # if (out.type != 1 && object$type == "regression") 
        # stop("'prob' or 'vote' not meaningful for regression")
        # print(missing(newdata))
        # print(out.type)
    if (missing(newdata) && missing(newgraph)) {
        if (out.type == 1) return(object$predicted)
        if (out.type == 2) return(object$prob)
    }
   	if (missing(newdata) || missing(newgraph)) {
   		warning("Both newdata and newgraph not specified, ignoring arguments")
        if (out.type == 1) return(object$predicted)
        if (out.type == 2) return(object$prob)
    }
    
    nsubj <- dim(newgraph)[1]
    if (is.null(newdata)) 	newdata <- data.frame(z=rep(1, nsubj ))

    
    ### if it's a formula - then put the data into a model frame
    if (inherits(object, "sig.sgraph.formula")) {
        # newdata <- data.frame(newdata, row.names=rownames(newdata))
        #rn <- row.names(newdata)
        Terms <- delete.response(object$terms)
        print(is(newdata))
        print(class(newdata))
        m <- model.frame.default(Terms, data=as.data.frame(newdata), na.action = na.omit, xlev = object$xlevels)
        x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        ### taking out intercept term
        # x <- X[, -which(colnames(X) == "(Intercept)")]
    } else {
        if (is.null(dim(newdata))) 
            dim(newdata) <- c(1, length(newdata))
        if (class(newdata) != "data.frame")
        	stop("newdata must be a data.frame or NULL")
         #newdata <- as.data.frame(newdata)
        x <- newdata
        if (nrow(x) == 0) 
            stop("newdata has 0 rows")
        if (any(is.na(x))) 
            stop("missing values in newdata")
        # x <- data.frame(rep(1, length(newdata)), newdata)
    }
    ## check if number of cols are same
    #####TOOK THIS OUT 19Apr12
    # train.ncol <- ncol(object$x)
    # vname <- colnames(object$x)
    # if (is.null(colnames(x))) {
        # if (ncol(x) != train.ncol) {
            # stop("number of variables in newdata does not match that in the training data")
        # }
    # } else {
        # if (any(!vname %in% colnames(x))) 
            # stop("variables in the training data missing in newdata")
        # x <- x[, vname, drop = FALSE]
    # }
    
    # took out 19Apr12
    # # # if (is.data.frame(x)) {
        # # # xfactor <- which(sapply(x, is.factor))
        # # # if (length(xfactor) > 0 && "xlevels" %in% names(object)) {
            # # # for (i in xfactor) {
                # # # if (any(!levels(x[[i]]) %in% object$xlevels[[i]])) 
                  # # # stop("New factor levels not present in the training data")
                # # # x[[i]] <- factor(x[[i]], levels = levels(x[[i]])[match(levels(x[[i]]), 
                  # # # object$xlevels[[i]])])
            # # # }
        # # # }
        # # # # cat.new <- sapply(x, function(x) if (is.factor(x) && 
            # # # # !is.ordered(x)) 
            # # # # length(levels(x))
        # # # # else 1)
        # # # # if (!all(object$forest$ncat == cat.new)) 
            # # # # stop("Type of predictors in new data do not match that of the training data.")
    # # # }
    
	##### Still need to work out data.frame for newdata

	ngroups <- length(object$graph.models)
	probs <- lapply(1:ngroups, function(ig) sapply(object$graph.models[[ig]], function(mod) predict(mod, newdata=newdata, type="response")))
    # for (ig in 1:ngroups){
    	# for (icol in 1:length(object$graph.models[[igrp]])){
    		# mod <- object$graph.models[[igrp]][[iedge]]
    		# if (object$dens != "multinom") predict(mod, newdata=x)
    		# else predict(mod, newdata=x, type="prob")
    	# }
    # }
    	runs <- vector(length=ngroups, mode="list")
    	dens <- object$dens
    	glevs <- NULL
			for (ig in 1:ngroups){
				for (icol in 1:length(object$graph.models[[ig]])) {
					G <- newgraph[, icol]
					# mod <- object$graph.models[[ig]][[icol]]
					Ghat <- probs[[ig]][[icol]]
		    		# if (object$dens != "multinom") Ghat <- predict(mod, newdata=x)
    				# else Ghat <- predict(mod, newdata=x, type="prob")
    				print(G)
					runs[[ig]][[icol]] <- dens.probs(G, Ghat, dens, glevs)
				}
			}
			    			
	# probs <- lapply(object$graph.models, function(mod) predict(mod, newdata=x, type="prob"))
	print(probs)
	stop("here")
#### STOPPED HERE
   	if (object$usepvals) {
		## use weighted sum on the log scale
		for (iicol in 1:ncol(probs)) col.probs[, iicol] <- log(probs[, iicol])*pvals[iicol]
		col.probs <- apply(col.probs, 1, function(edge) exp(sum(edge)))
	} else {
		col.probs <- apply(col.probs, 1, prod)
	}
	col.probs <- ifelse(col.probs==1, 1-eta, ifelse(col.probs < .Machine$double.eps, eta, col.probs))
#### STOPPED HERE
	
    
    priors <- predict(object$mod.priors, newdata=x, type='prob')
    
	for (icol in 1:ncol(probs)) {
		if (colnames(probs)[icol] != colnames(priors)[icol]) stop("Something wrong with predictions")	
		probs[, icol] <- probs[, icol]*priors[, icol]
	}
	#priors <- res$priors
	
	#print("Code block 6")
	## divide by total to get actual probability and not just numerator P(G = g | Y = y)
	## apply(probs, 1, sum) = P(G=g|Y=1)P(Y=1)+ ... + P(G=g|Y=y_n)P(Y=y_n)
	probs <- probs / apply(probs, 1, sum)
	#print("Head Probs")
	#print(head(probs))
	
	preds <- apply(probs, 1, function(x) which(x == max(x)))
	if (class(preds) != "integer") stop("Problem with Prediction")
	preds <- factor(colnames(probs)[preds], levels=levs)
	#preds <- levs[preds]
	print(length(preds))
	print(head(preds))
	print(table(preds, y))
	accuracy <- mean(preds == y)
	cat(sprintf("Accuracy is %3.2f%% \n", accuracy*100))



    # return("BLAH")
    
    # mdim <- ncol(x)
    # ntest <- nrow(x)
    #res

	
}