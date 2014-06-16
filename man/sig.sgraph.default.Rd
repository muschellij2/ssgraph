\name{sig.sgraph.default}
\alias{sig.sgraph.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sig.sgraph.default(graph, y = NULL, x = NULL, corr = FALSE, family = NULL, na.action = NULL, upper = FALSE, alternative = "two.sided", keep = 0.1, weight = FALSE, workspace = 4e+05, formula = NULL, data = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
%%     ~~Describe \code{graph} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{corr}{
%%     ~~Describe \code{corr} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{na.action}{
%%     ~~Describe \code{na.action} here~~
}
  \item{upper}{
%%     ~~Describe \code{upper} here~~
}
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
}
  \item{keep}{
%%     ~~Describe \code{keep} here~~
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{workspace}{
%%     ~~Describe \code{workspace} here~~
}
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (graph, y = NULL, x = NULL, corr = FALSE, family = NULL, 
    na.action = NULL, upper = FALSE, alternative = "two.sided", 
    keep = 0.1, weight = FALSE, workspace = 4e+05, formula = NULL, 
    data = NULL, ...) 
{
    checkPckg("MASS")
    checkPckg("gamlss.dist")
    checkPckg("nnet")
    print(alternative)
    adj <- sort(unique(c(graph)))
    nadj <- length(unique(c(graph)))
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
    if (is.null(family) & !corr) {
        cat("No Family Specified, assuming Multinomial\n")
        family <- list(family = "multinom")
    }
    if (weight) {
        if (is.character(family)) {
            if (!(family \%in\% c("multinom", "beta"))) {
                family <- get(family, mode = "function", envir = parent.frame())
            }
            else {
                family <- list(family = family)
            }
        }
        if (is.function(family)) 
            family <- family()
        if (is.null(family$family) & corr) {
            cat("No Family Specified and Corr=TRUE, binomial chosen")
            family <- binomial()
        }
        if (is.null(family$family)) {
            stop(paste("'family'", "family", "not recognized"))
        }
    }
    print(family)
    alternative <- char.expand(alternative, c("two.sided", "less", 
        "greater"))
    if (length(alternative) > 1L || is.na(alternative)) 
        stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
    dag <- dim(graph)
    cg <- class(graph)
    if (cg == "data.frame" & ("y" \%in\% colnames(graph) & is.null(y))) {
        y <- graph$y
        which.y <- which(colnames(graph) == "y")
        graph <- graph[, -which.y]
    }
    if (cg == "matrix" && ("y" \%in\% colnames(graph) & is.null(y))) {
        y <- graph[, "y"]
        which.y <- which(colnames(graph) == "y")
        graph <- graph[, -which.y]
    }
    graph <- convert.graph(graph, upper = upper)
    nullx <- FALSE
    form <- !is.null(formula)
    if (is.null(y)) 
        stop("No Y input")
    if (is.null(x)) {
        x <- matrix(rep(1, length(y)), ncol = 1)
        nullx <- TRUE
        if (form) 
            stop("No X specified or only intercept")
    }
    if (any(is.na(x))) 
        stop("NA not permitted in predictors")
    if (any(is.na(graph))) 
        stop("NA not permitted in graph")
    if (any(is.na(y))) 
        stop("NA not permitted in response")
    if (is.null(dim(x))) 
        dim(x) <- c(length(x), 1)
    if (nrow(graph) != nrow(x) | nrow(x) != length(y)) {
        stop("Number of graphs does not match N Subjects from X or Y")
    }
    dy <- dim(y)
    cy <- class(y)
    if (cy == "matrix") 
        y <- c(y)
    if (!(cy \%in\% c("matrix", "numeric", "character", "factor", 
        "integer", "logical"))) 
        stop("Class of Y unrecognized/too big - array")
    ngroups <- length(unique(y))
    if (ngroups < 2) 
        stop("Need at least 2 groups")
    if (alternative != "two.sided" & ngroups > 2) {
        cat("For groups/labels > 2, alternative is two sided\n")
        alternative <- "two.sided"
    }
    if (class(y) != "factor") 
        y <- factor(y)
    levs <- levels(y)
    if (ngroups == 2) {
        if (class(levs) == "character") 
            y <- as.numeric(factor(y, levels = levs)) - 1
        if (class(levs) == "factor") 
            y <- as.numeric(y, levels = levs) - 1
        if (!all(as.numeric(levs) == c(0, 1))) 
            stop("Label problems")
    }
    dg <- dim(graph)
    nvert <- dg[2]
    nsubj <- dg[1]
    if (length(dg) > 2) 
        stop("Problem with dimensions")
    print("in getpvals")
    pvals <- apply(graph, 2, function(edge) getpvals(edge, y = y, 
        workspace = workspace, weight = weight, alternative = alternative))
    if (is.function(keep)) {
        nkeep <- nvert
        pvals <- keep(pvals)/sum(keep(pvals))
    }
    else nkeep <- floor(nvert * keep)
    ssgraph <- graph.ind <- sort(order(pvals)[1:nkeep])
    if (cg == "array") {
        sig.graph <- matrix(0, nrow = dag[1], ncol = dag[2])
        if (upper) {
            up.tri <- upper.tri(sig.graph, diag = TRUE)
            sig.graph[up.tri][ssgraph] <- 1
        }
        else sig.graph[ssgraph] <- 1
    }
    else {
        sig.graph <- rep(0, nvert)
        sig.graph[ssgraph] <- 1
    }
    tgraph <- graph
    graph <- graph[, ssgraph]
    if (corr) {
        graph <- (graph + 1)/2
    }
    print("at get.probs.formula")
    print(list(family = family, weight = weight, ngroups = ngroups, 
        levs = levs, nkeep = nkeep, nsubj = nsubj, corr = corr, 
        nadj = nadj))
    res <- get.probs.formula(graph = graph, x = x, y = y, family = family, 
        weight = weight, ngroups = ngroups, levs = levs, nkeep = nkeep, 
        nsubj = nsubj, corr = corr, nadj = nadj, usepvals = is.function(keep), 
        pvals = pvals)
    probs <- res$probs
    priors <- res$priors
    colnames(probs) <- levs
    preds <- apply(probs, 1, function(x) which(x == max(x)))
    if (class(preds) != "integer") 
        stop("Problem with Prediction")
    preds <- factor(colnames(probs)[preds], levels = levs)
    print(table(preds))
    print(table(preds, y))
    accuracy <- mean(preds == y)
    cat(sprintf("Accuracy is \%3.2f\%\% \n", accuracy * 100))
    xlevels <- apply(x, 2, levels)
    ret <- list(ssgraph = sig.graph, priors = priors, prob = probs, 
        predicted = preds, y = y, accuracy = accuracy, pvals = pvals, 
        usepvals = is.function(keep), graph.models = res$mods, 
        graph.class = cg, graph.indices = graph.ind, formula = formula, 
        has.formula = form, xlevels = xlevels, terms = NULL, 
        x = x, family = family, dens = res$dens, upper = upper, 
        dg = dg)
    class(ret) <- c("sig.sgraph", class(ret))
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
