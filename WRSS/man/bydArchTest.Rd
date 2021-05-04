\name{bydArchTest}
\alias{bydArchTest}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("bydArchTest")}
\format{
  The format is:
List of 5
 $ statistic: Named num 22.1
  ..- attr(*, "names")= chr "Chi-squared"
 $ parameter: Named num 1
  ..- attr(*, "names")= chr "df"
 $ p.value  : Named num 2.63e-06
  ..- attr(*, "names")= chr "Chi-squared"
 $ method   : chr "ARCH LM-test;  Null hypothesis:  no ARCH effects"
 $ data.name: chr "hcs"
 - attr(*, "class")= chr "htest"
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
data(bydArchTest)
## maybe str(bydArchTest) ; plot(bydArchTest) ...
}
\keyword{datasets}
