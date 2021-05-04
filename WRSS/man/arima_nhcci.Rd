\name{arima_nhcci}
\alias{arima_nhcci}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("arima_nhcci")}
\format{
  The format is:
List of 14
 $ coef     : Named num -0.345
  ..- attr(*, "names")= chr "ar1"
 $ sigma2   : num 0.00387
 $ var.coef : num [1, 1] 0.0158
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr "ar1"
  .. ..$ : chr "ar1"
 $ mask     : logi TRUE
 $ loglik   : num 78.7
 $ aic      : num -153
 $ arma     : int [1:7] 1 0 0 0 4 2 0
 $ residuals: Time-Series [1:60] from 2003 to 2018: 0.000447 -0.00132 0.004449 -0.015084 0.020892 ...
 $ call     : language arima(x = nhcci_train.ts, order = c(1, 2, 0))
 $ series   : chr "nhcci_train.ts"
 $ code     : int 0
 $ n.cond   : int 0
 $ nobs     : int 58
 $ model    :List of 10
  ..$ phi  : num -0.345
  ..$ theta: num(0) 
  ..$ Delta: num [1:2] 2 -1
  ..$ Z    : num [1:3] 1 2 -1
  ..$ a    : num [1:3] -0.122 1.734 1.685
  ..$ P    : num [1:3, 1:3] -2.22e-16 1.25e-16 7.67e-17 1.25e-16 -7.18e-17 ...
  ..$ T    : num [1:3, 1:3] -0.345 1 0 0 2 ...
  ..$ V    : num [1:3, 1:3] 1 0 0 0 0 0 0 0 0
  ..$ h    : num 0
  ..$ Pn   : num [1:3, 1:3] 1.00 2.48e-17 5.95e-17 2.48e-17 -7.18e-17 ...
 - attr(*, "class")= chr "Arima"
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
data(arima_nhcci)
## maybe str(arima_nhcci) ; plot(arima_nhcci) ...
}
\keyword{datasets}
