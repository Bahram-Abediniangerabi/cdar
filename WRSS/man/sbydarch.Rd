\name{sbydarch}
\alias{sbydarch}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("sbydarch")}
\format{
  The format is:
List of 12
 $ call         : language dynlm(formula = ehatsq ~ L(ehatsq))
 $ terms        :Classes 'terms', 'formula'  language ehatsq ~ L(ehatsq)
  .. ..- attr(*, "variables")= language list(ehatsq, L(ehatsq))
  .. ..- attr(*, "factors")= int [1:2, 1] 0 1
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:2] "ehatsq" "L(ehatsq)"
  .. .. .. ..$ : chr "L(ehatsq)"
  .. ..- attr(*, "term.labels")= chr "L(ehatsq)"
  .. ..- attr(*, "order")= int 1
  .. ..- attr(*, "intercept")= int 1
  .. ..- attr(*, "response")= int 1
  .. ..- attr(*, ".Environment")=<environment: 0x0000021a79891228> 
 $ residuals    : Time-Series [1:197] from 2 to 198: 822 1749 -1123 -842 -635 ...
  ..- attr(*, "names")= chr [1:197] "2" "3" "4" "5" ...
 $ coefficients : num [1:2, 1:4] 650.4518 0.3344 142.7995 0.0674 4.555 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:2] "(Intercept)" "L(ehatsq)"
  .. ..$ : chr [1:4] "Estimate" "Std. Error" "t value" "Pr(>|t|)"
 $ aliased      : Named logi [1:2] FALSE FALSE
  ..- attr(*, "names")= chr [1:2] "(Intercept)" "L(ehatsq)"
 $ sigma        : num 1775
 $ df           : int [1:3] 2 195 2
 $ r.squared    : num 0.112
 $ adj.r.squared: num 0.107
 $ fstatistic   : Named num [1:3] 24.6 1 195
  ..- attr(*, "names")= chr [1:3] "value" "numdf" "dendf"
 $ cov.unscaled : num [1:2, 1:2] 6.48e-03 -1.42e-06 -1.42e-06 1.44e-09
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:2] "(Intercept)" "L(ehatsq)"
  .. ..$ : chr [1:2] "(Intercept)" "L(ehatsq)"
 $ frequency    : num 1
 - attr(*, "class")= chr [1:2] "summary.dynlm" "summary.lm"
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
data(sbydarch)
## maybe str(sbydarch) ; plot(sbydarch) ...
}
\keyword{datasets}
