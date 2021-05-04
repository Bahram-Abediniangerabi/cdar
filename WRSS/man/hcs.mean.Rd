\name{hcs.mean}
\alias{hcs.mean}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("hcs.mean")}
\format{
  The format is:
List of 14
 $ coefficients : Named num 66
  ..- attr(*, "names")= chr "(Intercept)"
 $ residuals    :'zoo' series from 1 to 198
  Data: num [1:198] 47.03 47.03 56.03 24.03 1.03 ...
  Index:  int [1:198] 1 2 3 4 5 6 7 8 9 10 ...
 $ effects      :'zoo' series from 1 to 198
  Data: num [1:198] -928.35 43.91 52.91 20.91 -2.09 ...
  Index:  int [1:198] 1 2 3 4 5 6 7 8 9 10 ...
 $ rank         : int 1
 $ fitted.values:'zoo' series from 1 to 198
  Data: num [1:198] 66 66 66 66 66 ...
  Index:  int [1:198] 1 2 3 4 5 6 7 8 9 10 ...
 $ assign       : int 0
 $ qr           :List of 5
  ..$ qr   : num [1:198, 1] -14.0712 0.0711 0.0711 0.0711 0.0711 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:198] "1" "2" "3" "4" ...
  .. .. ..$ : chr "(Intercept)"
  .. ..- attr(*, "assign")= int 0
  ..$ qraux: num 1.07
  ..$ pivot: int 1
  ..$ tol  : num 1e-07
  ..$ rank : int 1
  ..- attr(*, "class")= chr "qr"
 $ df.residual  : int 197
 $ call         : language dynlm(formula = hcs ~ 1)
 $ terms        :Classes 'terms', 'formula'  language hcs ~ 1
  .. ..- attr(*, "variables")= language list(hcs)
  .. ..- attr(*, "factors")= int(0) 
  .. ..- attr(*, "term.labels")= chr(0) 
  .. ..- attr(*, "order")= int(0) 
  .. ..- attr(*, "intercept")= int 1
  .. ..- attr(*, "response")= int 1
  .. ..- attr(*, ".Environment")=<environment: 0x0000021a76388640> 
 $ model        :'data.frame':	198 obs. of  1 variable:
  ..$ hcs:'zoo' series from 1 to 198
  Data: int [1:198, 1] 113 113 122 90 67 70 73 71 85 90 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr "hcs"
  Index:  int [1:198] 1 2 3 4 5 6 7 8 9 10 ...
  ..- attr(*, "terms")=Classes 'terms', 'formula'  language hcs ~ 1
  .. .. ..- attr(*, "variables")= language list(hcs)
  .. .. ..- attr(*, "factors")= int(0) 
  .. .. ..- attr(*, "term.labels")= chr(0) 
  .. .. ..- attr(*, "order")= int(0) 
  .. .. ..- attr(*, "intercept")= int 1
  .. .. ..- attr(*, "response")= int 1
  .. .. ..- attr(*, ".Environment")=<environment: 0x0000021a76388640> 
  .. .. ..- attr(*, "predvars")= language merge.zoo(hcs, retclass = "list", all = FALSE)
  .. .. ..- attr(*, "dataClasses")= Named chr "nmatrix.1"
  .. .. .. ..- attr(*, "names")= chr "hcs"
 $ index        : int [1:198] 1 2 3 4 5 6 7 8 9 10 ...
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
data(hcs.mean)
## maybe str(hcs.mean) ; plot(hcs.mean) ...
}
\keyword{datasets}
