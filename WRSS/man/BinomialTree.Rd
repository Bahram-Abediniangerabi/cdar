\name{BinomialTree}
\alias{BinomialTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BinomialTree(S, I, Time, r, sigma, dt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{I}{
%%     ~~Describe \code{I} here~~
}
  \item{Time}{
%%     ~~Describe \code{Time} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{dt}{
%%     ~~Describe \code{dt} here~~
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
function (S, I, Time, r, sigma, dt) 
{
    options(warn = -1)
    n = (Time/dt) + 1
    u = exp(sigma * sqrt(dt))
    d = 1/u
    k = exp(r * dt)
    p = (k - d)/(u - d)
    Df = exp(-r * dt)
    q = 1 - p
    OptionValue = (S * u^(0:n) * d^(n:0))
    offset = 1
    Tree = OptionValue = (abs(OptionValue) + OptionValue)/2
    {
        for (j in (n - 1):0) {
            Tree <- c(Tree, rep(0, times = n - j))
            for (i in 0:j) {
                OptionValue[i + offset] = max(((S * u^i * d^(abs(i - 
                  j)))), (p * OptionValue[i + 1 + offset] + (1 - 
                  p) * OptionValue[i + offset]) * Df)
                Tree = c(Tree, OptionValue[i + offset])
            }
        }
    }
    Tree = matrix(rev(Tree), byrow = FALSE, ncol = n + 1)
    colnames(Tree) <- paste(0:n, sep = "")
    rownames(Tree) <- paste(0:n, sep = "")
    dx = -0.025
    dy = 0.4
    cex = 1
    Tree_rounded = round(Tree, digits = 2)
    depth = ncol(Tree_rounded) - 1
    lablist <- as.vector(c(0:(depth - 1)))
    plot(x = c(1, depth), y = c(-depth, depth), col = 0, xlab = "Investment Horizon (Years)", 
        ylab = "Up or Down", labels = FALSE)
    axis(1, at = seq(1, n, by = 1), labels = FALSE)
    axis(2)
    text(seq(1, n, by = 1), par("usr")[3] - 0.2, labels = lablist, 
        pos = 1, xpd = TRUE)
    points(x = 1, y = 0)
    text(1 + dx, 0 + dy, deparse(Tree_rounded[1, 1]), cex = cex)
    for (i in 1:(depth - 1)) {
        y = seq(from = -i, by = 2, length = i + 1)
        x = rep(i, times = length(y)) + 1
        points(x, y, col = 1)
        for (j in 1:length(x)) text(x[j] + dx, y[j] + dy, deparse(Tree_rounded[length(x) + 
            1 - j, i + 1]), cex = cex)
        y = (-i):i
        x = rep(c(i + 1, i), times = 2 * i)[1:length(y)]
        lines(x, y, col = 2)
    }
    Cashflow = matrix(0, nrow = n + 1, ncol = n + 1)
    for (i in n:1) {
        for (j in n:1) {
            Cashflow[i, j] = -I + Tree[i, j] + ((p * Cashflow[i, 
                j + 1]) + ((1 - p) * Cashflow[i + 1, j + 1]))/(1 + 
                r)
        }
    }
    for (i in 0:n + 1) {
        for (j in 0:n + 1) {
            if (i > j) {
                Cashflow[i, j] = 0
            }
            else {
                Cashflow[i, j] = Cashflow[i, j]
            }
        }
    }
    colnames(Cashflow) <- paste(0:n, sep = "")
    rownames(Cashflow) <- paste(0:n, sep = "")
    Cashflow = round(Cashflow, digits = 2)
    DecisionTree = as.data.frame(Cashflow)
    for (i in 1:n) {
        for (j in 1:n) {
            if (DecisionTree[i, j] > 0) {
                DecisionTree[i, j] = "Invest"
            }
            else {
                DecisionTree[i, j] = "Wait"
            }
        }
    }
    for (i in 0:n + 1) {
        for (j in 0:n + 1) {
            if (i > j) {
                DecisionTree[i, j] = "_"
            }
            else {
                DecisionTree[i, j] = DecisionTree[i, j]
            }
        }
    }
    colnames(DecisionTree) <- paste("Year ", 0:n, sep = "")
    rownames(DecisionTree) <- paste(0:n, sep = "")
    param = c(S, Time, 1 + r, sigma, n - 1, u, 1/u, p, 1 - p)
    (Parameters <- structure(param, names = c("S", "Time", "Rf", 
        "sigma", "n", "Up", "Down", "Pi_Up", "Pi_Down")))
    print("Model Parameters:")
    print(Parameters, digits = 2)
    print("Binomial Tree:")
    print(Tree_rounded, digits = 2)
    print("Cashflow:")
    print(Cashflow)
    print("Decision Tree:")
    print(DecisionTree[0:n, 0:n])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
