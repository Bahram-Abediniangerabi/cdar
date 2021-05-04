\name{byd.arch}
\alias{byd.arch}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("byd.arch")}
\format{
  The format is:
List of 10
 $ order        : Named num [1:2] 0 1
  ..- attr(*, "names")= chr [1:2] "p" "q"
 $ coef         : Named num [1:2] 499.21 0.993
  ..- attr(*, "names")= chr [1:2] "a0" "a1"
 $ n.likeli     : num 921
 $ n.used       : int 198
 $ residuals    : num [1:198] NA 0.984 1.063 0.728 0.725 ...
 $ fitted.values: num [1:198, 1:2] NA 114.8 114.8 123.6 92.4 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "sigt" "-sigt"
 $ series       : chr "hcs"
 $ frequency    : num 1
 $ call         : language garch(x = hcs, order = c(0, 1))
 $ vcov         : num [1:2, 1:2] 4.66e+05 -1.68e+02 -1.68e+02 1.05e-01
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:2] "a0" "a1"
  .. ..$ : chr [1:2] "a0" "a1"
 - attr(*, "class")= chr "garch"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(byd.arch)
## maybe str(byd.arch) ; plot(byd.arch) ...
}
\keyword{datasets}
