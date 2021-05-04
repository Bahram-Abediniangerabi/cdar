\name{hcs.ARCH}
\alias{hcs.ARCH}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("hcs.ARCH")}
\format{
  The format is:
List of 15
 $ coefficients : Named num [1:2] 650.452 0.334
  ..- attr(*, "names")= chr [1:2] "(Intercept)" "L(ehatsq)"
 $ residuals    : Time-Series [1:197] from 2 to 198: 822 1749 -1123 -842 -635 ...
  ..- attr(*, "names")= chr [1:197] "2" "3" "4" "5" ...
 $ effects      : Time-Series [1:197] from 2 to 198: -13751 -8802 -1305 -867 -624 ...
  ..- attr(*, "names")= chr [1:197] "(Intercept)" "L(ehatsq)" "" "" ...
 $ rank         : int 2
 $ fitted.values: Time-Series [1:197] from 2 to 198: 1390 1390 1700 843 651 ...
 $ assign       : int [1:2] 0 1
 $ qr           :List of 5
  ..$ qr   : num [1:197, 1:2] -14.0357 0.0712 0.0712 0.0712 0.0712 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:197] "2" "3" "4" "5" ...
  .. .. ..$ : chr [1:2] "(Intercept)" "L(ehatsq)"
  .. ..- attr(*, "assign")= int [1:2] 0 1
  ..$ qraux: num [1:2] 1.07 1.04
  ..$ pivot: int [1:2] 1 2
  ..$ tol  : num 1e-07
  ..$ rank : int 2
  ..- attr(*, "class")= chr "qr"
 $ df.residual  : int 195
 $ xlevels      : Named list()
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
 $ model        :'data.frame':	197 obs. of  2 variables:
  ..$ ehatsq   : Time-Series [1:197] from 2 to 198: 2211.37 3138.83 577.21 1.05 16.2 ...
  ..$ L(ehatsq): Time-Series [1:197] from 2 to 198: 2211.37 2211.37 3138.83 577.21 1.05 ...
  ..- attr(*, "terms")=Classes 'terms', 'formula'  language ehatsq ~ L(ehatsq)
  .. .. ..- attr(*, "variables")= language list(ehatsq, L(ehatsq))
  .. .. ..- attr(*, "factors")= int [1:2, 1] 0 1
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:2] "ehatsq" "L(ehatsq)"
  .. .. .. .. ..$ : chr "L(ehatsq)"
  .. .. ..- attr(*, "term.labels")= chr "L(ehatsq)"
  .. .. ..- attr(*, "order")= int 1
  .. .. ..- attr(*, "intercept")= int 1
  .. .. ..- attr(*, "response")= int 1
  .. .. ..- attr(*, ".Environment")=<environment: 0x0000021a79891228> 
  .. .. ..- attr(*, "predvars")= language merge.zoo(ehatsq, L(ehatsq), retclass = "list", all = FALSE)
  .. .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "numeric"
  .. .. .. ..- attr(*, "names")= chr [1:2] "ehatsq" "L(ehatsq)"
 $ index        : num [1:197] 2 3 4 5 6 7 8 9 10 11 ...
 $ frequency    : num 1
 $ twostage     : logi FALSE
 - attr(*, "class")= chr [1:2] "dynlm" "lm"
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
data(hcs.ARCH)
## maybe str(hcs.ARCH) ; plot(hcs.ARCH) ...
}
\keyword{datasets}
