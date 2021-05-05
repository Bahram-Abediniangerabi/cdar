#' @title Probabilistic Life Cycle Cost Analysis
#' @description This R function will provide probabilistic LCCA for a projects with several uncertain initail cost components.
#' @param comp1 A vector for the first initial cost component  #Example for comp1 = rtriangle(n_loop, a =10, b = 200, c = 30)
#' @param comp2 A vector for the second initial cost component #Example for comp2 = rnorm(n_loop, mean=10, sd=3)
#' @param comp3 A vector for the third initial cost component  #Example for comp3 = rlnorm(n_loop, meanlog = 10, sdlog = 1)
#' @param comp4 A vector for the forth initial cost component  #Example for comp4 = runif(n_loop, min = 0, max = 25)
#' @param comp5 A vector for the fifth initial cost component  #Example for comp5 = runif(n_loop, min = 0, max = 25)
#' @param recurring_cost The annual reccuring cost of the project
#' @param r Rate of Return
#' @param n Investment Horizon (yearly)
#' @param n_loop Number of Monte Carlo Iterations
#' @return Returns the histogram plot of all recurring cost components (inputs) and histogram, kernel density, and cumulative density function plots of net present value and equivalent uniform value (outputs) for the project.
#' @examples Lcca_MonteCarlo(comp1 = rnorm(1000, mean=34, sd=8), comp2 = rnorm(1000, mean=10, sd=3),
#' comp3 = rlnorm(1000, meanlog = 10, sdlog = 1), r =0.05, recurring_cost = 100,  n=10, n_loop = 1000)
#' @export
Lcca_MonteCarlo=function(comp1 = NA, comp2 = NA,comp3 = NA, comp4 = NA, comp5 = NA, recurring_cost = NA, r, n, n_loop){
  all=list(comp1,comp2,comp3,comp4,comp5,recurring_cost,r,n_loop, n)

  #n_loop=10000
  #comp1 = rtriangle(n_loop, a =10, b = 200, c = 180)
  #comp2 = rnorm(n_loop, mean=10, sd=10)
  #comp3 = rlnorm(n_loop, meanlog = 10, sdlog = 1)
  #comp4 = rnorm(n_loop, mean=100, sd=20)
  #comp5 = runif(n_loop, min = 100, max = 300)
  #Lcca_MonteCarlo(comp1 = comp1, comp2 = comp2, comp3 = comp3, comp4 = comp4, comp5 = comp5,  r=0.05, n_loop=n_loop, n=50, recurring_cost = 100)

  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  if (missing(comp1)){comp1 = rep(0, length.out=n_loop)} else {comp1 = comp1}
  if (missing(comp2)){comp2 = rep(0, length.out=n_loop)} else {comp2 = comp2}
  if (missing(comp3)){comp3 = rep(0, length.out=n_loop)} else {comp3 = comp3}
  if (missing(comp4)){comp4 = rep(0, length.out=n_loop)} else {comp4 = comp4}
  if (missing(comp5)){comp5 = rep(0, length.out=n_loop)} else {comp5 = comp5}
  if (missing(recurring_cost)){recurring_cost_list = rep(0, length.out=n)} else {
    recurring_cost_list = rep(recurring_cost*(((1+r)^n)-1)/(r*(1+r)^n), length.out = n)}

  npv_total_cost =  comp1 + comp2 + comp3 + comp4 + comp5 + recurring_cost_list
  euv_total_cost = npv_total_cost*(r*(1+r)^n)/(((1+r)^n)-1)

  #print(npv_recurring_cost)
  hist(comp1, col = "Gray",xlab="Comp1", main="Histogram of Comp1")
  hist(comp2, col = "Gray",xlab="Comp1", main="Histogram of Comp2")
  hist(comp3, col = "Gray",xlab="Comp1", main="Histogram of Comp3")
  hist(comp4, col = "Gray",xlab="Comp1", main="Histogram of Comp4")
  hist(comp5, col = "Gray",xlab="Comp1", main="Histogram of Comp5")

  hist(npv_total_cost, col = "Gray",xlab="Net Present Values", main="Histogram of Net Present Values")
  plot(density(npv_total_cost),main="Kernel Density of Net Present Values",col="blue")
  plot(ecdf(npv_total_cost),main="Cumulative Density Function of Net Present Values",xlab="Net Present Values",col="blue")

  #print(euv_total_cost)
  hist(euv_total_cost, col = "Gray",xlab="Equivalent Uniform Values", main="Histogram of Equivalent Uniform Values")
  plot(density(euv_total_cost),main="Kernel Density of Equivalent Uniform Values",col="blue")
  plot(ecdf(euv_total_cost),main="Cumulative Density Function of Equivalent Uniform Values",xlab="Equivalent Uniform Values",col="blue")
}

#' @title Real Options Pricing using Binomial Lattice
#' @description This R function will provide real option prices for underlying assets (Infrastructure Projects, Buildings, etc) using binomial lattice.
#' @param S State Variable
#' @param I Investment
#' @param dt Time Intervals within a Year
#' @param r Rate of Return
#' @param sigma Fluctuations in the Price of State Variable
#' @param Time Investment Horizon (yearly)
#' @return Returns a binomial tree for the state variable "S", cashflow matrix calculated from the binomial tree and the investment cost, decision matrix for investment for different situations through the investment horizon, and a binomial tree plot.
#' @examples BinomialTree(S=50, I=50, Time=5, r=0.2, sigma=0.4, dt=1)
#' @export
BinomialTree <- function(S, I, Time, r, sigma, dt)
  {
  options(warn=-1)
    # Parameters:
    n = (Time / dt)+1
    u  = exp(sigma*sqrt(dt))
    d  = 1 / u
    k  = exp(r*dt)
    p  = (k - d) / (u - d)
    Df = exp(-r*dt)
    q  = 1-p
    # Algorithm:
    OptionValue = (S*u^(0:n)*d^(n:0))
    offset = 1
    Tree = OptionValue = (abs(OptionValue)+OptionValue)/2
    {for (j in (n-1):0) {
      Tree <-c(Tree, rep(0, times=n-j))
      for (i in 0:j) {
        OptionValue[i+offset] =
          max(((S*u^i*d^(abs(i-j)))),
              (p*OptionValue[i+1+offset] +
                 (1-p)*OptionValue[i+offset]) * Df )
        Tree = c(Tree, OptionValue[i+offset]) } } }
    Tree = matrix(rev(Tree), byrow = FALSE, ncol = n+1)
    colnames(Tree) <- paste(0:n, sep = "")
    rownames(Tree) <- paste( 0:n, sep = "")

    # Binomial Lattice for State Variable:
    dx = -0.025
    dy = 0.4
    cex = 1
    Tree_rounded = round(Tree, digits = 2)
    depth = ncol(Tree_rounded)-1
    lablist<-as.vector(c(0:(depth-1)))
    plot(x = c(1,depth), y = c(-depth, depth), col = 0, xlab = "Investment Horizon (Years)", ylab = "Up or Down", labels = FALSE)
    axis(1, at=seq(1, n, by=1),labels = FALSE)
    axis(2)
    text(seq(1, n, by=1), par("usr")[3] - 0.2, labels = lablist, pos = 1, xpd = TRUE)
    points(x = 1, y = 0)
    text(1+dx, 0+dy, deparse(Tree_rounded[1, 1]), cex = cex)
    for (i in 1:(depth-1) ) {
      y = seq(from = -i, by = 2, length = i+1)
      x = rep(i, times = length(y))+1
      points(x, y, col = 1)
      for (j in 1:length(x))
        text(x[j]+dx, y[j]+dy, deparse(Tree_rounded[length(x)+1-j,i+1]), cex = cex)
      y = (-i):i
      x = rep(c(i+1,i), times = 2*i)[1:length(y)]
      lines(x, y, col = 2)
    }
    # Cashflow Matrix
    Cashflow=matrix(0, nrow = n+1, ncol = n+1)
    for(i in n:1){
      for(j in n:1){
        Cashflow[i,j]=-I+Tree[i,j]+((p*Cashflow[i,j+1])+((1-p)*Cashflow[i+1,j+1]))/(1+r)
      }}
    # Triangle(Cashflow)
    for(i in 0:n+1){
      for(j in 0:n+1){
        if(i>j){Cashflow[i,j]=0
        }else{
          Cashflow[i,j]=Cashflow[i,j]
        }
      }
    }
    colnames(Cashflow) <- paste(0:n, sep = "")
    rownames(Cashflow) <- paste( 0:n, sep = "")
    Cashflow = round(Cashflow, digits = 2)
    # DecisionTree Matrix
    DecisionTree = as.data.frame(Cashflow)
    for(i in 1:n){
      for(j in 1:n){
        if(DecisionTree[i,j] > 0){
          DecisionTree[i,j]='Invest'
        } else {
          DecisionTree[i,j]='Wait'
        }
      }
    }
    # Triangle(DecisionTree)
    for(i in 0:n+1){
      for(j in 0:n+1){
        if(i>j){DecisionTree[i,j]='_'
        }else{
          DecisionTree[i,j]=DecisionTree[i,j]
        }
      }
    }
    colnames(DecisionTree) <- paste("Year ", 0:n, sep = "")
    rownames(DecisionTree) <- paste( 0:n, sep = "")

    ### Function Outputs
    # Printing Model Parameters
    param = c(S,Time,1+r,sigma,n-1,u,1/u,p,1-p)
    (Parameters <- structure(param,names=c("S","Time","Rf","sigma","n","Up","Down","Pi_Up","Pi_Down")))
    print("Model Parameters:")
    print(Parameters,digits=2)
    print("Binomial Tree:")
    print(Tree_rounded,digits=2)
    print("Cashflow:")
    print(Cashflow)
    print("Decision Tree:")
    print(DecisionTree[0:n,0:n])
    }

#' @title Binomial Real Options Pricing with Monte Carlo Simulation
#' @description This R function will provide real option prices and the probability of investment within the investment horizon using binomial lattice and monte carlo simulations.
#' @param S State Variable
#' @param I Investment
#' @param dt Time Intervals within a Year
#' @param r Rate of Return
#' @param sigma Fluctuations in the Price of State Variable
#' @param Time Investment Horizon (yearly)
#' @param MC_loops Number of Monte Carlo Simulations
#' @return Returns a binomial tree for the state variable "S", cashflow matrix calculated from the binomial tree and the investment cost, decision matrix for investment for different situations through the investment horizon, a binomial tree plot, and the likelihood of implementation plot.
#' @examples BinomialTree_MC(S=10, I=100, Time=10, r=.01, sigma=0.6, dt=1, MC_loops = 1000)
#' @export
BinomialTree_MC <- function(S, I, Time, r, sigma, dt, MC_loops)
  {
  options(warn=-1)
  # Parameters:
  n = (Time / dt)+1
  u  = exp(sigma*sqrt(dt))
  d  = 1 / u
  k  = exp(r*dt)
  p  = (k - d) / (u - d)
  Df = exp(-r*dt)
  # Algorithm:
  OptionValue = (S*u^(0:n)*d^(n:0))
  offset = 1
  Tree = OptionValue = (abs(OptionValue)+OptionValue)/2
  {for (j in (n-1):0) {
    Tree <-c(Tree, rep(0, times=n-j))
    for (i in 0:j) {
      OptionValue[i+offset] =
        max(((S*u^i*d^(abs(i-j)))),
            (p*OptionValue[i+1+offset] +
               (1-p)*OptionValue[i+offset]) * Df )
      Tree = c(Tree, OptionValue[i+offset]) } } }
  Tree = matrix(rev(Tree), byrow = FALSE, ncol = n+1)
  colnames(Tree) <- paste(0:n, sep = "")
  rownames(Tree) <- paste( 0:n, sep = "")

  # Tree Output
  # Binomial Lattice for State Variable:
  dx = -0.025
  dy = 0.4
  cex = 1
  Tree_rounded = round(Tree, digits = 2)
  depth = ncol(Tree_rounded)-1
  lablist<-as.vector(c(0:(depth-1)))
  plot(x = c(1,depth), y = c(-depth, depth), col = 0, xlab = "Investment Horizon (Years)", ylab = "Up or Down", labels = FALSE)
  axis(1, at=seq(1, n, by=1),labels = FALSE)
  axis(2)
  text(seq(1, n, by=1), par("usr")[3] - 0.2, labels = lablist, pos = 1, xpd = TRUE)
  points(x = 1, y = 0)
  text(1+dx, 0+dy, deparse(Tree_rounded[1, 1]), cex = cex)
  for (i in 1:(depth-1) ) {
    y = seq(from = -i, by = 2, length = i+1)
    x = rep(i, times = length(y))+1
    points(x, y, col = 1)
    for (j in 1:length(x))
      text(x[j]+dx, y[j]+dy, deparse(Tree_rounded[length(x)+1-j,i+1]), cex = cex)
    y = (-i):i
    x = rep(c(i+1,i), times = 2*i)[1:length(y)]
    lines(x, y, col = 2)
  }
  #Cashflow Matrix
  Cashflow=matrix(0, nrow = n+1, ncol = n+1)
  for(i in n:1){
    for(j in n:1){
      Cashflow[i,j]=-I+Tree[i,j]+((p*Cashflow[i,j+1])+((1-p)*Cashflow[i+1,j+1]))/(1+r)
    }}
  #Triangle(Cashflow)
  for(i in 0:n+1){
    for(j in 0:n+1){
      if(i>j){Cashflow[i,j]=0
      }else{
        Cashflow[i,j]=Cashflow[i,j]
      }
    }
  }
  colnames(Cashflow) <- paste(0:n, sep = "")
  rownames(Cashflow) <- paste( 0:n, sep = "")
  Cashflow = round(Cashflow, digits = 2)
  #DecisionTree Matrix
  DecisionTree = as.data.frame(Cashflow)
  for(i in 1:n){
    for(j in 1:n){
      if(DecisionTree[i,j] > 0){
        DecisionTree[i,j]='Invest'
      } else {
        DecisionTree[i,j]='Wait'
      }
    }
  }
  #Triangle(DecisionTree)
  for(i in 0:n+1){
    for(j in 0:n+1){
      if(i>j){DecisionTree[i,j]='_'
      }else{
        DecisionTree[i,j]=DecisionTree[i,j]
      }
    }
  }
  colnames(DecisionTree) <- paste("Year ", 0:n, sep = "")
  rownames(DecisionTree) <- paste( 0:n, sep = "")

  Decision_Mat= matrix(0, 1, n-1)
  for (i in 0:MC_loops){
    Sample_path = rbinom(n-1,1,0.5)
    Decision_Path = vector()
    x_position = 1
    y_position = 1
    for (x in Sample_path) {
      if (x==1) {
        x_position = x_position
        y_position = y_position+1
        Decision_Path = append(Decision_Path,DecisionTree[x_position,y_position])}
      else if (x==0) {
        x_position = x_position+1
        y_position = y_position+1
        Decision_Path = append(Decision_Path,DecisionTree[x_position,y_position])}
    }
    Decision_Mat = rbind.data.frame(Decision_Path,Decision_Mat)
  }
  Decision_Mat$Never <- apply(Decision_Mat, 1, function(x) length(which(x=="Wait")))
  #Counting loops with no Investment
  Full_Wait <- count(Decision_Mat$Never==(n-1))
  Investment_Probability_Table <- ldply(Decision_Mat, function(c) sum(c=="Invest"))
  Investment_Probability_Table["V1"] = Investment_Probability_Table["V1"]/(MC_loops)
  Investment_Probability_Table[n,"V1"] = Full_Wait[2,2]/(MC_loops)
  #Plotting Investment Probabilities
  X_axis = paste("Year ", 1:n, sep = "")
  X_axis[n] = "Never"
  X_axis<-factor(X_axis, levels = X_axis)
  Label = Investment_Probability_Table[,2]
  Invest.per.year.item <-
    ggplot(Investment_Probability_Table, aes(x=as.array(X_axis), y=Label)) +
    geom_bar(stat="identity", colour="Navy") + xlab("Year") + ylab("Likelihood of Implementation") +
    geom_text(aes(label=Label), position=position_dodge(width=0.9), vjust=-0.25)
  param = c(S,Time,1+r,sigma,n-1,u,1/u,p,1-p)
  (Parameters <- structure(param,names=c("S","Time","Rf","sigma","n","Up","Down","Pi_Up","Pi_Down")))

  ######Function Outputs######
  print("Model Parameters:")
  print(Parameters,digits=2)
  #State Variable Binomial Tree
  print("Binomial Tree:")
  print(Tree_rounded,digits=2)
  #Cashflow Matrix
  print("Cashflow:")
  print(Cashflow)
  #Decision Tree
  print("Decision Tree:")
  print(DecisionTree[0:n,0:n])
  #print(Investment_Probability_Table)
  print(Invest.per.year.item)
}

#' @title Cashflow Diagram
#' @description This R function will provide a cashflow diagram with the given cashflow streams and their times.
#' @param cf_t0 Cashflows at time 0
#' @param cf Cashflow vector
#' @param times Cashflow times
#' @return Returns the Cashflow digram for the given cashflow streams
#' @examples cfd(cf_t0 = -1000, cf=c(200,200,1000), times = c(1,1,2))
#' @export
cfd=function(cf_t0, cf,times){
  all=list(cf_t0, cf,times)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Numeric
  if(!is.vector(cf) | !is.numeric(cf)) stop("cf must be a numeric vector.")
  if(!is.vector(times) | !is.numeric(times)) stop("times must be a numeric vector.")
  #NA
  if(any(is.na(cf)) | any(is.na(times))) stop("Cannot input any variables as NA.")
  #Infinite
  if(any(cf==Inf)) stop("Cannot have infinite in cf.")
  if(any(times==Inf)) stop("Cannot have infinite in times.")
  #Positive
  if(any(times<0)) stop("Cannot have negative values in times.")
  if(length(cf)==0 | length(times)==0 ) stop("Missing cf or times information.")
  if(length(cf) != length(times)) stop("Number of payments not equal to number of period values.")

  d=unique(times[which(duplicated(times)==T)])
  if(length(d)>0){ cfd=rep(0,length(d))
  for(r in 1:length(d)) cfd[r]=sum(cf[which(times==d[r])])}
  p=rep(0,max(times))
  p[times]=cf
  if(length(d)>0) {for(r in 1:length(d)) p[d[r]]=cfd[r]}

  x.pv=p[which(p!=0)]
  if(length(x.pv)==0) x.pv=0
  x.t=unique(times)[order(unique(times))]
  plot(0,0,type="n",axes=F,ann=F,xlim = c(-1,max(x.t)+1),ylim=c(0, 10))
  axis(1,at=c(0,x.t),labels=c(0,x.t),line=-5)
  axis(1,at=0,labels=cf_t0,line=-7,tck=0)
  par(mgp=c(3,2,0))
  par(col.axis="blue")
  par(col.axis="black")
  text(x.t,rep(5.5,max(x.t)),labels=x.pv,cex=.85)
  text((max(x.t))/2,10,labels="Time Diagram",cex=1.2)
  legend(0,1,legend="Period",lty=0,bty="n")
  par(mgp=c(3,1,0))
}

#' @title Net Present Value
#' @description This R function will calculate NPV for the given cashflow streams.
#' @param cf_t0 Cashflows at time 0
#' @param cf Cashflow vector
#' @param times Cashflow times
#' @param r rate of return
#' @return Returns the NPV for the given cashflow streams
#' @examples NPV(cf_t0 = -1000 , cf = c(100,300)  , times = c(3,2)  , r = 0.05)
#' @export
NPV=function(cf_t0,cf,times,r){
  all=list(cf_t0,cf,times,r)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(list(cf_t0,r),length) != 1)==T) stop("cf_t0 and r   must be of length 1.")
  #Numeric
  if(!is.vector(cf) | !is.numeric(cf)) stop("cf must be a numeric vector.")
  if(!is.vector(times) | !is.numeric(times)) stop("times must be a numeric vector.")
  if(!is.numeric(cf_t0) | !is.numeric(r)) stop("cf_t0 and r must be numeric.")
  #NA
  if(any(is.na(cf)) | any(is.na(times)) | any(is.na(c(cf_t0,r)))) stop("Cannot input NA for any variables.")
  #Infinite
  if(any(cf==Inf) | any(times==Inf) | cf_t0==Inf | r==Inf) stop("Cannot input infinite for any variables.")
  #Positive
  if(any(times<=0)) stop("Cannot have negative values in times.")
  if(r<0) stop("r cannot be negative.")

  if(length(cf)==0 | length(times)==0 ) stop("Not enough cash flow information.")
  if(length(cf) != length(times)) stop("Amount of cash flows not equal to amount of time values.")

  cf_t0=cf_t0
  pv=sum(cf/(1+r)^times)
  npv=pv+cf_t0

  return(npv)
}

#' @title Time Value of Money
#' @description This R function can calculate present value and future value the given information.
#' @param pv present value
#' @param fv future value
#' @param n number of periods
#' @param r nominal interest rate convertible ic times per period
#' @param ic interest conversion frequency per period
#' @return Returns the NPV or FV for the given information
#' @examples tvm(pv=50, n=5, r=.04)
#' @export
tvm=function(pv=NA,fv=NA,n=NA,r=NA,ic=1) {
  ###Example: TVM(pv=50,n=5,r=.04)
  all=list(pv,fv,n,r,ic)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(pv,fv,n,r,ic)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("pv, fv, n, r, and ic must be numeric.")
  #NA
  nalist=list(ic)
  if(any(lapply(nalist,is.na)==T)) stop("Cannot input ic as NA.")

  NA.Neg=array(c(pv,fv,n,r,ic))
  NA.Neg.Str=c("pv","fv","n","r","ic")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}

  test=c(pv,fv,n,r)
  if(length(test[which(is.na(test))])>=2) stop("Too many unknowns.")
  if(length(test[which(is.na(test))])==0) stop("Must have one unknown.")
  if(!is.na(fv) & pv>fv & !is.na(pv)) stop("PV cannot be greater than FV.")

  nom1=NA
  if(!is.na(r)) {int=(1+r/ic)^ic-1
  if(ic!=1) nom1=r}

  if(is.na(pv)){pv=fv/(1+int)^(n);solve="PV"}
  if(is.na(fv)){fv=pv*(1+int)^(n);solve="FV"}
  if(is.na(n)){n=log(fv/pv)/log(1+int);solve="n"}
  if(is.na(r)){int=(fv/pv)^(1/n)-1;solve="r"
  if(ic!=1) nom1=((1+int)^(1/ic)-1)*ic}

  out=c(pv,fv,n,int,nom1)
  m.out=matrix(out,nrow=length(out))
  rownames(m.out)=c("PV","FV","Periods",
                    "Eff Rate",paste("r^(",round(ic,2),")",sep=""))
  na=apply(m.out, 1, function(x) all(is.na(x)))
  m.out=as.matrix(m.out[!na,])
  colnames(m.out)=c("TVM")
  return(m.out)
}

#' @title Equalnet Uniform Value or Annuity
#' @description This R function calculate annuity for the present value, future value, number of payments/periods, amount of the first payment, the payment increment amount per period, and/or the interest rate for an arithmetically growing annuity.
#' @param pv present value
#' @param fv future value
#' @param n number of periods
#' @param r nominal interest rate convertible ic times per period
#' @param ic interest conversion frequency per period
#' @param a amount of the first payment
#' @param q payment increment amount per period
#' @param pf Payment frequency - number of payments per year
#' @param imm option for annuity immediate or annuity due, default is immediate (TRUE)
#' @param plot option to display a time diagram of the payments
#' @return Returns Annuity
#' @examples euv(pv=1000, fv=NA, n=10, a=NA, q=0, r=.05, ic=1, pf=1, imm=FALSE)
#' @export
euv=function(pv=NA,fv=NA,n=NA,a=NA,q=NA,r=NA,ic=1,pf=1,imm=TRUE,plot=FALSE){

  all=list(pv,fv,n,a,q,r,ic,pf,imm,plot)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(pv,fv,n,a,r,ic,pf,q)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("pv, fv, n, a, q, r, ic, and pf must be numeric.")
  #NA
  nalist=list(ic,pf,imm,plot)
  if(any(lapply(nalist,is.na)==T)) stop("Cannot input ic, pf, imm, and plot as NA.")
  #Logical
  stopifnot(is.logical(imm),is.logical(plot))

  num=c(pv,fv,n,a,r,ic,pf)
  NA.Neg=array(num)
  NA.Neg.Str=c("pv","fv","n","a","r","ic","pf")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}

  test=c(n,r);t1=length(test[which(is.na(test))])
  test2=c(pv,fv);t2=length(test2[which(is.na(test2))])
  test3=c(a,q);t3=length(test3[which(is.na(test3))])
  test4=c(n,r,a,q);t4=length(test4[which(is.na(test4))])
  if(t2 == 0) stop("Cannot specify both pv and fv.")
  if(t1 >= 2) stop("Too many unknowns.")
  if(t3 >= 2) stop("Cannot have both a and q unknown")
  if(t4==0 & t2==0) stop("No unknowns specified.")
  if(t4==1 & t2==2) stop("Too many unknowns.")
  if(t1==0 & t3==0 & t2==1) stop("one of n, a, q, or r must be unknown if either pv or fv is known.")

  if(!is.na(n) & n !=round(n)) stop("n must be a positive integer")
  if(!is.na(a) & !is.na(q) & a+q<=0) stop("q is too small.")
  if(is.na(n)) sn=T else sn=F
  if(!is.na(a) & !is.na(q) & !is.na(n) & a+q*(n-1)<0) stop("Payments become negative.")

  nom1=NA;nom2=NA

  if(is.na(r)){
    if(imm==T) {
      if(is.na(pv)) {r=round(Re(polyroot(c(-fv+a+q*(n-1),rev(seq(from=a,by=q,length.out=n-1)))))-1,8)
      r=r[which(round((a*(1-(1+r)^-n)/r+q*((1-(1+r)^-n)/r-n*(1+r)^-n)/r)*(1+r)^n,0)==round(fv,0))]}
      if(is.na(fv)) {r=round(1/Re(polyroot(c(-pv,(seq(from=a,by=q,length.out=n)))))-1,8)
      r=r[which(round((a*(1-(1+r)^-n)/r+q*((1-(1+r)^-n)/r-n*(1+r)^-n)/r),2)==round(pv,2))]}
    }
    if(imm==F) {
      if(is.na(pv)) {r=round(Re(polyroot(c(-fv,rev(seq(from=a,by=q,length.out=n)))))-1,8)
      r=r[which(round((a*(1-(1+r)^-n)/r+q*((1-(1+r)^-n)/r-n*(1+r)^-n)/r)*(1+r)^(n+1),0)==round(fv,0))]}
      if(is.na(fv)) {r=round(1/Re(polyroot(c(-pv+a,(seq(from=a+q,by=q,length.out=n-1)))))-1,8)
      r=r[which(round((a*(1-(1+r)^-n)/r+q*((1-(1+r)^-n)/r-n*(1+r)^-n)/r)*(1+r),0)==round(pv,0))]}
    }
    if(length(r)>1) {if(any(r>0)) r=r[which(r>0)] else r=r[1];r=r[1]}
    eff.r=(1+r)^pf-1
    n.r=(1+eff.r)^(1/ic)-1
    if(ic != 1) nom1=n.r*ic
    if(pf != 1 & pf != ic) nom2=r*pf
    r=n.r*ic
  }

  eff.r=(1+r/ic)^(ic)-1
  int=(1+eff.r)^(1/pf)-1
  if(imm==T) imm_r=1 else imm_r=1+int
  if(ic !=1) nom1=((1+eff.r)^(1/ic)-1)*ic
  if(pf != ic & pf !=1) nom2=((1+eff.r)^(1/pf)-1)*pf

  if(is.na(a)){
    if(is.na(fv)) a=(pv/imm_r-q*((1-(1+int)^-n)/int-n*(1+int)^-n)/int)/((1-(1+int)^-n)/int)
    if(is.na(pv)) a=(fv/imm_r/(1+int)^n-q*((1-(1+int)^-n)/int-n*(1+int)^-n)/int)/((1-(1+int)^-n)/int)
  }

  if(is.na(q)){
    if(is.na(fv)) q=(pv/imm_r-a*(1-(1+int)^-n)/int)/(((1-(1+int)^-n)/int-n*(1+int)^-n)/int)
    if(is.na(pv)) q=(fv/imm_r/(1+int)^n-a*(1-(1+int)^-n)/int)/(((1-(1+int)^-n)/int-n*(1+int)^-n)/int)
  }

  if(is.na(n)){
    if(is.na(pv)) h=fv else h=pv
    for(g in 0:1000){
      if(is.na(fv)) h=h-imm_r*(a+q*g)/(1+int)^g else h=h-(a+q*g)
      if(h<=0) {g=g+1;break}
    }
    n=seq(from=g/4,to=3*g,by=.0001)
    if(is.na(fv)) n=n[which(round(imm_r*(a*(1-(1+int)^-n)/int+q*((1-(1+int)^-n)/int-n*(1+int)^-n)/int),0)==round(pv,0))]
    if(is.na(pv)) n=n[which(round(imm_r*(a*(1-(1+int)^-n)/int+q*((1-(1+int)^-n)/int-n*(1+int)^-n)/int)*(1+int)^n,1)==round(fv,1))]
    n=median(n)
  }

  if(is.na(pv)) {
    pv=(a*(1-(1+int)^-n)/int+q*((1-(1+int)^-n)/int-n*(1+int)^-n)/int)*imm_r}
  if(is.na(fv)) {
    fv=(a*(1-(1+int)^-n)/int+q*((1-(1+int)^-n)/int-n*(1+int)^-n)/int)*imm_r*(1+int)^n}

  if(plot==T){
    if(sn==F){
      plot.new()
      plot(0,0,type="n",axes=F,ann=F,xlim = c(0,ceiling(n)+1),ylim=c(0, 10))
      axis(1,at=seq(0,ceiling(n)),labels=seq(0,ceiling(n)),line=-8)
      if(imm==T){
        text(seq(1,ceiling(n)),rep(5.5,ceiling(n)),labels=seq(round(a,2),by=round(q,2),length.out=ceiling(n)),cex=.75)
      }
      if(imm==F){
        text(seq(0,ceiling(n)-1),rep(5.5,ceiling(n)),labels=seq(round(a,2),by=round(q,2),length.out=ceiling(n)),cex=.75)
      }
      if(pf != 1) axis(1,at=seq(0,ceiling(n),by=pf),labels=seq(0,ceiling(n),by=pf)/pf,line=-5,col="blue")
      text(.05*n,7.2,labels=round(pv,2),cex=.85)
      text(.05*n,7.7,labels="PV")
      text(n-.05*n,7.2,labels=round(fv,2),cex=.85)
      text(n-.05*n,7.7,labels="FV")
      text(n/2,10,labels="Time Diagram",cex=1.2)
      text(n/2,9.4,labels=if(q != 0) "Arithmetic Annuity" else "Level Annuity",cex=1)
      text(n*.4,8.5,bquote("Eff Rate "== .(round(eff.r,4))),cex=.85)
      if(pf != 1) text(n*.4,7.8,bquote( r^(.(pf))== .(round(int*pf,4))),cex=.85)
      r.ic=((1+eff.r)^(1/ic)-1)*ic
      if(ic != 1 & ic != pf) text(n*.4,7.1,bquote( r^(.(ic))== .(round(r.ic,4))),cex=.85)
      if(pf==1) legend(0,1,legend="Years",lty=1,bty="n") else legend(0,1,legend=c("Periods","Years"),lty=1,col=c("black","blue"),bty="n",ncol=2)
    } else warning("No time diagram is provided when solving for n.\n") }

  if(pf==1) {years=n;n=NA} else years=n/pf
  out=c(pv,fv,a,q,eff.r,nom1,nom2,n,years)
  m.out=matrix(out,nrow=length(out))
  rownames(m.out)=c("PV","FV","A","Q","Eff Rate",paste("r^(",round(ic,2),")",sep=""),
                    paste("r^(",round(pf,2),")",sep=""),"Periods","Years")
  na=apply(m.out, 1, function(x) all(is.na(x)))
  m.out=as.matrix(m.out[!na,])
  if(round(q,2)==0) m="Level Annuity" else m="Arithmetic Annuity"
  colnames(m.out)=m
  return(m.out)
}

