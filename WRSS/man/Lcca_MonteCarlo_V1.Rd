\name{Lcca_MonteCarlo_V1}
\alias{Lcca_MonteCarlo_V1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lcca_MonteCarlo_V1(comp1 = NA, comp2 = NA, comp3 = NA, comp4 = NA, comp5 = NA, recurring_cost = NA, r, n, n_loop)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comp1}{
%%     ~~Describe \code{comp1} here~~
}
  \item{comp2}{
%%     ~~Describe \code{comp2} here~~
}
  \item{comp3}{
%%     ~~Describe \code{comp3} here~~
}
  \item{comp4}{
%%     ~~Describe \code{comp4} here~~
}
  \item{comp5}{
%%     ~~Describe \code{comp5} here~~
}
  \item{recurring_cost}{
%%     ~~Describe \code{recurring_cost} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{n_loop}{
%%     ~~Describe \code{n_loop} here~~
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
function (comp1 = NA, comp2 = NA, comp3 = NA, comp4 = NA, comp5 = NA, 
    recurring_cost = NA, r, n, n_loop) 
{
    all = list(comp1, comp2, comp3, comp4, comp5, recurring_cost, 
        r, n_loop, n)
    if (any(lapply(all, is.null) == T)) 
        stop("Cannot input any variables as NULL.")
    if (missing(comp1)) {
        comp1 = rep(0, length.out = n_loop)
    }
    else {
        comp1 = comp1
    }
    if (missing(comp2)) {
        comp2 = rep(0, length.out = n_loop)
    }
    else {
        comp2 = comp2
    }
    if (missing(comp3)) {
        comp3 = rep(0, length.out = n_loop)
    }
    else {
        comp3 = comp3
    }
    if (missing(comp4)) {
        comp4 = rep(0, length.out = n_loop)
    }
    else {
        comp4 = comp4
    }
    if (missing(comp5)) {
        comp5 = rep(0, length.out = n_loop)
    }
    else {
        comp5 = comp5
    }
    if (missing(recurring_cost)) {
        recurring_cost_list = rep(0, length.out = n)
    }
    else {
        recurring_cost_list = rep(recurring_cost * (((1 + r)^n) - 
            1)/(r * (1 + r)^n), length.out = n)
    }
    npv_total_cost = comp1 + comp2 + comp3 + comp4 + comp5 + 
        recurring_cost_list
    euv_total_cost = npv_total_cost * (r * (1 + r)^n)/(((1 + 
        r)^n) - 1)
    hist(comp1, col = "Gray", xlab = "Comp1", main = "Histogram of Comp1")
    hist(comp2, col = "Gray", xlab = "Comp1", main = "Histogram of Comp2")
    hist(comp3, col = "Gray", xlab = "Comp1", main = "Histogram of Comp3")
    hist(comp4, col = "Gray", xlab = "Comp1", main = "Histogram of Comp4")
    hist(comp5, col = "Gray", xlab = "Comp1", main = "Histogram of Comp5")
    hist(npv_total_cost, col = "Gray", xlab = "Net Present Values", 
        main = "Histogram of Net Present Values")
    plot(density(npv_total_cost), main = "Kernel Density of Net Present Values", 
        col = "blue")
    plot(ecdf(npv_total_cost), main = "Cumulative Density Function of Net Present Values", 
        xlab = "Net Present Values", col = "blue")
    hist(euv_total_cost, col = "Gray", xlab = "Equivalent Uniform Values", 
        main = "Histogram of Equivalent Uniform Values")
    plot(density(euv_total_cost), main = "Kernel Density of Equivalent Uniform Values", 
        col = "blue")
    plot(ecdf(euv_total_cost), main = "Cumulative Density Function of Equivalent Uniform Values", 
        xlab = "Equivalent Uniform Values", col = "blue")
    hist(comp3)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
