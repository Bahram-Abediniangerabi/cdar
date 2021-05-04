Lcca_MonteCarlo_V1 <-
function(comp1 = NA, comp2 = NA,comp3 = NA, comp4 = NA, comp5 = NA, recurring_cost = NA, r, n, n_loop){
  all=list(comp1,comp2,comp3,comp4,comp5,recurring_cost,r,n_loop, n)

  #n_loop=10000
  #comp1 = rtriangle(n_loop, a =10, b = 200, c = 180)
  #comp2 = rnorm(n_loop, mean=10, sd=10)
  #comp3 = rlnorm(n_loop, meanlog = 10, sdlog = 1)
  #comp4 = rnorm(n_loop, mean=100, sd=20)
  #comp5 = runif(n_loop, min = 100, max = 300)
  #Lcca_MonteCarlo_V1(comp1 = comp1, comp2 = comp2, comp3 = comp3, comp4 = comp4, comp5 = comp5,  r=0.05, n_loop=n_loop, n=50, recurring_cost = 100)

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
  hist(comp1, col = "Gray",xlab="Comp1", main="Histogram of Comp1")
  hist(comp2, col = "Gray",xlab="Comp1", main="Histogram of Comp2")
  hist(comp3, col = "Gray",xlab="Comp1", main="Histogram of Comp3")
  hist(comp4, col = "Gray",xlab="Comp1", main="Histogram of Comp4")
  hist(comp5, col = "Gray",xlab="Comp1", main="Histogram of Comp5")
  #print(npv_recurring_cost)
  hist(npv_total_cost, col = "Gray",xlab="Net Present Values", main="Histogram of Net Present Values")
  plot(density(npv_total_cost),main="Kernel Density of Net Present Values",col="blue")
  plot(ecdf(npv_total_cost),main="Cumulative Density Function of Net Present Values",xlab="Net Present Values",col="blue")

  #print(euv_total_cost)
  hist(euv_total_cost, col = "Gray",xlab="Equivalent Uniform Values", main="Histogram of Equivalent Uniform Values")
  plot(density(euv_total_cost),main="Kernel Density of Equivalent Uniform Values",col="blue")
  plot(ecdf(euv_total_cost),main="Cumulative Density Function of Equivalent Uniform Values",xlab="Equivalent Uniform Values",col="blue")
  hist(comp3)
  }
