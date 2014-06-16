rm(list=ls())
source("convert.graph.R")
source("get.probs.R")
source("get.probs.formula.R")
source("sig.sgraph.formula.R")
source("sig.sgraph.R")
source("predict.sig.sgraph.R")
library(nnet)
#### Test data
graph <- agraph <- array(rbinom(1000 * 100, 1, prob=0.5), dim=c(10, 10, 1000))
ylabels <- rbinom(1000, 5, 0.5)
#graph <- t(apply(agraph, 3, function(x) return(x[upper.tri(x)])))
#graph <- data.frame(cbind(graph, ylabels=ylabels))
upper <- FALSE
alternative = "two.sided"
keep.pct=0.1
weight <- TRUE
corr <- FALSE
workspace=400000



ylabels <- rbinom(1000, 5, 0.5)
ylab <- rep("Control", 1000)
ylab[ylabels == 1] <- "ADHD1"
ylab[ylabels == 2] <- "ADHD2"
ylab[ylabels == 3] <- "ADHD3"
ylab[ylabels == 4] <- "ADHD4"

ylabels2 <- floor(runif(1000, 0, 5))
ylab2 <- rep("Control", 1000)
ylab2[ylabels2 == 1] <- "ADHD1"
ylab2[ylabels2 == 2] <- "ADHD2"
ylab2[ylabels2 == 3] <- "ADHD3"
ylab2[ylabels2 == 4] <- "ADHD4"

ylab <- factor(ylab, levels= c("Control", "ADHD1", "ADHD2", "ADHD3", "ADHD4"))
X <- matrix(rnorm(10*1000), nrow=1000)

df <- data.frame(y=ylab2, X)
graph <- agraph
# Q <- sig.sgraph.formula(formula=y ~ 1  , data=df, graph=agraph, keep=.1)
###Q <- sig.sgraph.formula(formula=y ~ 1 , data=df, graph=matrix(rgamma(1000*100, 2, 3), nrow=1000, ncol=100), family=Gamma(), keep=.1)
object <- Q <- sig.sgraph.formula(formula=y ~ 1, data=df, graph=agraph, family=poisson(), keep=.1)
Q$ssgraph
predict(Q, newdata = df, newgraph = agraph)

object <- Q <- sig.sgraph.formula(formula=y ~ 1, data=df, graph=agraph, family=poisson(), keep=.1, upper=TRUE)
Q$ssgraph

object <- Q <- sig.sgraph.formula(formula=y ~ 1, data=df, graph=agraph, family=poisson(), upper=TRUE)
Q$ssgraph
predict(Q, newdata = df, newgraph = agraph)


graph <- t(apply(agraph, 3, function(x) return(x[upper.tri(x)])))
object <- Q <- sig.sgraph.formula(formula=y ~ 1, data=df, graph=graph, family=poisson(), keep=.1)

predict(Q, newdata = df, newgraphs = agraph)
Q <- sig.sgraph(y=df$y, x=NULL, data=df, graph=agraph, family=poisson(), keep=.1)
age <- rnorm(1000, mean=25)
Q <- sig.sgraph(y=df$y, x=age, data=df, graph=agraph, family=poisson(), keep=.1)


Q <- sig.sgraph.formula(formula=y ~ 1, data=df, graph=agraph, family=poisson(), keep=.1, corr=TRUE)
Q <- sig.sgraph.formula(formula=y ~ 1  , data=df, graph=agraph, family=quasipoisson(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ 1  , data=df, graph=agraph, family=gaussian(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ 1  , data=df, graph=agraph, family=binomial(), keep=.1)

Q <- sig.sgraph.formula(formula=y ~ X , data=df, graph=matrix(rgamma(1000*100, 2, 3), nrow=1000, ncol=100), family=Gamma(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ X , data=df, graph=matrix(rgamma(1000*100, 2, 3), nrow=1000, ncol=100), family=Gamma(), keep=.1, corr=TRUE)
Q <- sig.sgraph.formula(formula=y ~ X  , data=df, graph=agraph, family=poisson(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ X  , data=df, graph=agraph, family=quasipoisson(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ X  , data=df, graph=agraph, family=gaussian(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ X  , data=df, graph=agraph, family=binomial(), keep=.1)

graph <- agraph <- array(rpois(1000 * 100, lambda=0.5), dim=c(10, 10, 1000))
Q <- sig.sgraph.formula(formula=y ~ 1, data=df, graph=agraph, family=poisson(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ 1, data=df, graph=agraph, family=poisson(), keep=.1, corr=TRUE)
Q <- sig.sgraph.formula(formula=y ~ 1  , data=df, graph=agraph, family=quasipoisson(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ 1  , data=df, graph=agraph, family=gaussian(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ 1  , data=df, graph=agraph, family=binomial(), keep=.1)

Q <- sig.sgraph.formula(formula=y ~ X , data=df, graph=matrix(rgamma(1000*100, 2, 3), nrow=1000, ncol=100), family=Gamma(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ X  , data=df, graph=agraph, family=poisson(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ X  , data=df, graph=agraph, family=quasipoisson(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ X  , data=df, graph=agraph, family=gaussian(), keep=.1)
Q <- sig.sgraph.formula(formula=y ~ X  , data=df, graph=agraph, family=binomial(), keep=.1)

Q <- sig.sgraph(y=ylab, x=NULL, data=df, graph=agraph)
ylabels[523] <- NA


ylabels[523] <- NA
ylabels2 <- rbinom(1000, 1, 0.5)
graph <- agraph <- array(rbinom(1000 * 100, 1, prob=0.5), dim=c(10, 10, 1000))
x <- sig.sgraph(graph=agraph, y=ylabels)
x <- sig.sgraph(agraph, y=ylabels, x=X)


df$y[523] <- 1
Q <- sig.sgraph.formula(formula=y ~ 1, data=df, graph=agraph)

graph <- agraph <- array(runif(1000 * 100, min=-1, max=1), dim=c(10, 10, 1000))
x <- sig.sgraph(agraph, ylabels=ylabels, corr=TRUE, weight=TRUE)
x <- sig.sgraph(agraph, ylabels=ylabels2, corr=TRUE, weight=TRUE)

x <- sig.sgraph(agraph, ylabels=ylabels, corr=TRUE, weight=FALSE)
x <- sig.sgraph(agraph, ylabels=ylabels2, corr=TRUE, weight=FALSE)

x <- sig.sgraph(agraph, ylabels=ylabels, corr=FALSE)
x <- sig.sgraph(agraph, ylabels=ylabels2, corr=FALSE)

x <- sig.sgraph(agraph, ylabels=ylabels, corr=FALSE, weight=TRUE)
x <- sig.sgraph(agraph, ylabels=ylabels2, corr=FALSE, weight=TRUE)

x <- sig.sgraph(graph, ylabels=ylabels, alternative="less")
x <- sig.sgraph(graph, ylabels=ylabels, alternative="great")

graph <- t(apply(agraph, 3, function(x) return(x[upper.tri(x)])))
graph <- data.frame(cbind(graph, ylabels=ylabels))

x <- sig.sgraph(graph)
x <- sig.sgraph(graph, ylabels=ylabels, alternative="less")
x <- sig.sgraph(graph, ylabels=ylabels, alternative="great")

graph <- agraph <- array(rbinom(1000 * 100, 1, prob=0.5), dim=c(10, 10, 1000))
ylabels <- rbinom(1000, 5, 0.5)
x <- sig.sgraph(agraph, ylabels=ylabels, workspace=10000000)
x <- sig.sgraph(agraph, ylabels=ylabels, alternative="less")
x <- sig.sgraph(agraph, ylabels=ylabels, alternative="great")


ylabels <- c(rep(0, 100), rep(1, 100))
mat0 <- matrix(runif(400, min=0, max=1), nrow=20)
mat1 <- matrix(runif(400, min=0, max=1), nrow=20)

p0 <- c(mat0)
p1 <- c(mat1)

graph0 <- sapply(p0, function(x) rbinom(100, 1, x))
graph1 <- sapply(p1, function(x) rbinom(100, 1, x))

graph <- data.frame(rbind(graph0, graph1))
graph$ylabels <- ylabels
x <- sig.sgraph(graph)