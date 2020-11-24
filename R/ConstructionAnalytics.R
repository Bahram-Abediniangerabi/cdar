#' @title Real Options Pricing using Binomial Lattice
#' @description This R function will provide real option prices for underlying assets (Infrastructure Projects, Buildings, etc) using binomial lattice.
#' @param S State Variable
#' @param I Investment
#' @param Time Time Intervals within a Year
#' @param r Rate of Return
#' @param sigma Fluctuations in the Price of State Variable
#' @param n Investment Horizon (yearly)
#' @return NULL
#' @examples BinomialTree(S=60, I=50,Time=1, r=.1,sigma=0.2,n=5)
#' @export
BinomialTree <- function(S, I, Time, r, sigma, n)
  { # A function implemented by Diethelm Wuertz

    ### Note:
    # Calculates option prices from the Cox-Ross-Rubinstein Binomial tree model.
    #   The model described here is a version of the CRR Binomial Tree model.

    ### FUNCTION:
    # Parameters:
    n  = n+1
    dt = Time / (n-1)
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

    # Binomial Lattice for State Variable:
    dx = -0.025
    dy = 0.4
    cex = 1
    Tree_rounded = round(Tree, digits = 2)
    depth = ncol(Tree_rounded)
    plot(x = c(1,depth), y = c(-depth+1, depth-1), col = 0)
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

    Sample_path = rbinom(n,1,0.5)
    Decision = c()
    Decision_Path = vector()
    x_position = 0
    y_position = 0
    for (i in Sample_path) {
      if (i==1) {
        x_position = x_position
        y_position = y_position+1}
      if (i==0) {
        x_position = x_position+1
        y_position = y_position+1}
      Decision_Path=append(Decision_Path,DecisionTree[x_position,y_position])
    }
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
    print(DecisionTree[0:n,0:n])}

#' @title Binomial Real Options Pricing with Monte Carlo Simulation
#' @description This R function will provide real option prices and the probability of investment within the investment horizon using binomial lattice and monte carlo simulations.
#' @param S State Variable
#' @param I Investment
#' @param Time Time Intervals within a Year
#' @param r Rate of Return
#' @param sigma Fluctuations in the Price of State Variable
#' @param n Investment Horizon (yearly)
#' @param MC_loops Number of Monte Carlo Simulations
#' @return NULL
#' @examples BinomialTree_MC(S=70, I=100, Time=1, r=.05, sigma=0.3, n=10, MC_loops = 1000)
#' @export
BinomialTree_MC <- function(S, I, Time, r, sigma, n, MC_loops)
  {   # A function implemented by Diethelm Wuertz

  # Description:
  #   Calculates option prices from the Cox-Ross-Rubinstein
  #   Binomial tree model.

  # Note:
  #   The model described here is a version of the CRR Binomial
  #   Tree model.

  # Example:
  # BinomialTree_MC(S=70, I=100, Time=1, r=.05, sigma=0.3, n=10, MC_loops = 1000)

  # FUNCTION:

  # Parameters:
  n  = n+1
  dt = Time / (n-1)
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
  # Tree Output:

  # Binomial Lattice for State Variable:
  dx = -0.025
  dy = 0.4
  cex = 1
  Tree_rounded = round(Tree, digits = 2)
  depth = ncol(Tree_rounded)
  plot(x = c(1,depth), y = c(-depth+1, depth-1), col = 0)
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
  require(plyr)
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
  require(ggplot2)
  Invest.per.year.item <-
    ggplot(Investment_Probability_Table, aes(x=as.array(X_axis), y=V1)) +
    geom_bar(stat="identity", colour="Navy") + xlab("Year") + ylab("Probability of Investment") +
    geom_text(aes(label=V1), position=position_dodge(width=0.9), vjust=-0.25)


  ######Function Outputs######
  #Printing Model Parameters
  param = c(S,Time,1+r,sigma,n-1,u,1/u,p,1-p)
  (Parameters <- structure(param,names=c("S","Time","Rf","sigma","n","Up","Down","Pi_Up","Pi_Down")))
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
  #print(Sample_path)
  #print(Decision_Path)
  #print(Decision_Mat)
  print(Investment_Probability_Table)
  print(Invest.per.year.item)}
