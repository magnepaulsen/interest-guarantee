#Part 1
#Black & Scholes Price of Interest Rate Guarantee with Continuous Rebalancing of Portfolio

InterestGuarantee <- function(t, T, g, r, sigma, lambda){
  #Function to calculate the interest rate guarantee price  on a portfolio consisting of a stock and a bond (standard Black & Scholes assumption)
  #Parameters:
  #t: Start time
  #T: End time
  #g: guaranteed interest rate
  #r: risk free interest rate
  #sigma: volatility 
  #lambda: proportion of stocks in the portfolio
  
  K = exp(g*(T-t)) - (1-lambda)*exp(r*(T-t))
  d1 = (log(lambda/K) + (r+0.5*(sigma^2))*(T-t))/
    (sigma*sqrt(T-t))
  d2 = d1 - sigma*sqrt(T-t)
  P = K*exp(-r*(T-t))*pnorm(-d2)-lambda*pnorm(-d1)
  P
}


#Initalize values for the interest rate guarantee calculation:
g = 0.0275
r = 0.025
mu = 0.1
sigma = 0.2
lambda = 0.3

lambdas = seq(0.05, 1, 0.05)
sigmas = seq(0.05, 0.5, 0.05)

#Matrix of interest rate guarantees
InterestGuarantees = matrix(nrow=length(lambdas),
                           ncol=length(sigmas))


#Calculates the interest rate gurantee for different values of sigma and lambda:

for (i in 1:length(lambdas))
{
  for (j in 1:length(sigmas))
  {
    InterestGuarantees[i,j] = InterestGuarantee(0, 1, g,
                                         r, sigmas[j], lambdas[i])
  }
}


#Graph of the interest rate guarantee surface

postscript("avkGarFlate3.eps", horizontal=FALSE, height=6, width=8)
persp(lambdas, sigmas, InterestGuarantees, theta=-30,
      phi=30, xlab="Stock Proportion", ylab="Volatility",
      zlab="", expand=0.5, ticktype="detailed", shade=0.5,
      ltheta=120)
dev.off()

sigma = 0.2
lambda = 0.3

#Calculate how the Interest Rate Guarantee price varies with the guaranteed interest and the risk free interest rate
rs = seq(0.005, 0.1, 0.005)
gs = seq(0.005, 0.1, 0.005)
InterestGuaranteesRG = matrix(nrow=length(rs), ncol=length(gs))

for (i in 1:length(rs))
{
  for (j in 1:length(gs))
  {
    InterestGuaranteesRG[i,j] = InterestGuarantee(0, 1, gs[j],
                                rs[i], sigma, lambda)
  }
}


postscript("avkGarFlate3.eps", horizontal=FALSE, height=6, width=8)
persp(rs, gs, InterestGuaranteesRG, theta=-30, phi=30,
      xlab="Risk Free Interest Rate", ylab="Guaranteed Interest Rate",
      zlab="", expand=0.5, ticktype="detailed", shade=0.5,
      ltheta=120)
dev.off()


#Part 2
#Defined Contribution pension with and without interes rate guarantee

g = 0.0275 #Guaranteed interest rate
r = 0.025 #Risk free interest rate
sigma = 0.2 #Yearly stock return volatility
mu = 0.1 #Yearly stock return mean (drift)
lambda = 0.3 #Proportion of stocks in the portfolio

n = 100000 #Number of Monte Carlo simulations
T = 20 #Contract duration (years)
Pi = 1 #Yearly contribution

FU = matrix(nrow=n, ncol=T+1)
FG = matrix(nrow=n, ncol=T+1)
Psi = matrix(nrow=n, ncol=T+1)
S = matrix(nrow=n, ncol=T+1)

S[,] = exp(rnorm(n*(T+1), mu -
                   (sigma^2)/2, sigma))
Pg = InterestGuarantee(0, 1, g, r, sigma, lambda)
FU[,1] = Pi + Pg
FG[,1] = Pi
Psi[,1] = ((FG[,1]/FU[,1]) - 1)*100

for (i in 1:T+1){
  FU[,i] = (lambda*S[,i] + (1-lambda)*exp(r))*(Pi+Pg+FU[,i-1])
  FG[,i] = (pmax(exp(g),lambda*S[,i]+(1-lambda)*exp(r)))*(Pi+FG[,i-1])
  Psi[,i] = ((FG[,i]/FU[,i])-1)*100
}

postscript("psi1.eps", horizontal=FALSE, height=6, width=8)
par(mfrow=c(2,2))
plot(density(Psi[, 2+1]), ylab="Density",
     xlab="", sub="t=2", main="")
plot(density(Psi[, 5+1]), ylab="Density",
     xlab="", sub="t=5", main="")
plot(density(Psi[, 10+1]), ylab="Density",
     xlab="", sub="t=10", main="")
plot(density(Psi[, T+1]), ylab="Density",
     xlab="", sub="t=20", main="")
dev.off()

postscript("psi2.eps", horizontal=FALSE, height=6, width=8)
par(mfrow=c(2,2))
plot(density(FU[,2+1]), type="l",
     ylab="Density", xlab="",
     sub="t=2", main="", ylim=c(0,5), lty=2)
lines(density(FG[,2+1]), type="l", lty=1)
#legend("topright", c("With Guarantee", "Without Guarantee"), lty=c(1,2))
plot(density(FU[,5+1]), type="l", ylab="Density",
     xlab="", sub="t=5", main="", ylim=c(0,0.8),
     lty=2)
lines(density(FG[,5+1]), type="l", lty=1)
#legend("topright", c("With Guarantee", "Without Guarantee"), lty=c(1,2))
plot(density(FU[,10+1]), type="l", ylab="Density",
     xlab="", sub="t=10", main="", lty=2,
     ylim=c(0,0.3))
lines(density(FG[,10+1]), type="l", lty=1)
#legend("topright", c("With Guarantee", "Without Guarantee"), lty=c(1,2))
plot(density(FU[,20+1]), type="l", ylab="Density",
     xlab="", sub="t=20", main="", ylim=c(0,0.07),
      lty=2)
lines(density(FG[,20+1]), type="l", lty=1)
#legend("topright", c("With Guarantee", "Without Guarantee"), lty=c(1,2))
dev.off()


#Calculates the extreme values (VaR and conditional VaR) for the different scenarios

years = c(3, 6, 11, 21)
varPsi = matrix(nrow=4, ncol=2)
cvarPsi = matrix(nrow=4, ncol=2)
varFU = matrix(nrow=4, ncol=2)
cvarFU = matrix(nrow=4, ncol=2)
varFG = matrix(nrow=4, ncol=2)
cvarFG = matrix(nrow=4, ncol=2)
for (i in 1:4){
  varPsi[i,1] = quantile(Psi[,years[i]], 0.05)
  varPsi[i,2] = quantile(Psi[,years[i]], 0.01)
  varFU[i,1] = quantile(FU[,years[i]], 0.05)
  varFU[i,2] = quantile(FU[,years[i]], 0.01)
  varFG[i,1] = quantile(FG[,years[i]], 0.05)
  varFG[i,2] = quantile(FG[,years[i]], 0.01)
  cvarPsi[i,1] = mean(sort(Psi[, years[i]])
                      [seq(1, 0.05*n)])
  cvarPsi[i,2] = mean(sort(Psi[, years[i]])
                      [seq(1,0.01*n)])
  cvarFU[i,1] = mean(sort(FU[, years[i]])
                     [seq(1,0.05*n)])
  cvarFU[i,2] = mean(sort(FU[, years[i]])
                     [seq(1,0.01*n)])
  cvarFG[i,1] = mean(sort(FG[, years[i]])
                     [seq(1,0.05*n)])
  cvarFG[i,2] = mean(sort(FG[, years[i]])
                     [seq(1,0.01*n)])
}



#Finds P(Psi>0)
PsiG = vector(length=4)
plass = vector(length=4)
sortPsi = matrix(nrow=n, ncol=4)
sortPsi[,1] = sort(Psi[,2+1])
sortPsi[,2] = sort(Psi[,5+1])
sortPsi[,3] = sort(Psi[,10+1])
sortPsi[,4] = sort(Psi[,20+1])

for(i in 1:4){
  for (j in 1:n){
    if(sortPsi[j,i]>0){
      plass[i] = j-1
      break
    }
  }
  PsiG[i] = (n-plass[i])/n
}



#Part 3
#Defined Contribution Pension with Portfolio Rebalancing on Fixed Time Steps


n = 100000 #Number of Monte Carlo simulations
h = c(6,12,50,250) #Number of trading days (rebalancing)

Pi = 1 #Yearly contribution
g = 0.0275 #Guaranteed interest rate
r = 0.025 #Risk free interest rate
sigma = 0.2 #Yearly stock return volatility
mu = 0.1 #Yearly stock return mean (drift)
lambda = 0.3 #Proportion of stocks in the portfolio

stock = matrix(nrow=n, ncol=h)
hedge_error = matrix(0, n, length(h))
years = vector(length=n)

for(i in 1:length(h)){
  Returns = matrix(exp(rnorm(h[i]*n, (mu-sigma^2/2)*
                           (1/h[i]),sigma*sqrt(1/h[i]))), n, h[i])
  S = lambda
  
  for(k in 1:n){
    years[k] = prod(Returns[k,1:h[i]])
  }
  
  K = exp(g) - (1-lambda)*exp(r)
  a = -pnorm(-(log(S/K)+(r+(sigma^2)/2))/sigma)
  b = K*exp(-r)*pnorm(-(log(S/K)+(r-(sigma^2)/2))/sigma)
  
  for (t in 2:h[i])
  {
    S = S*Returns[,t]
    Vf = a*S + b*exp(r/h[i])
    a = -pnorm(-(log(S/K)+(r+(sigma^2)/2)*(1-t/h[i]))
               /(sigma*sqrt(1-t/h[i])))
    b = K*exp(-r*(1-t/h[i]))*pnorm(-(log(S/K) +
                                       (r-(sigma^2)/2)*(1-t/h[i]))/(sigma*sqrt(1-t/h[i])))
    Ve = a*S + b
    hedge_error[,i] = hedge_error[,i] + (Ve-Vf)*
      exp(-r*(t/h[i]))
  }
  
  yearly_return = lambda*years + (1-lambda)*exp(r)
  gar = Pi*pmax((exp(g)-yearly_return),0)
  Vs = a*S*Returns[,h[i]] + b*exp(r/h)
  hedge_error[,i] = hedge_error[,i] + (gar-Vs)*exp(-r)
}

var = matrix(nrow=length(h), ncol=2)
cvar = matrix(nrow=length(h), ncol=2)
PG = vector(length=length(h))

meanHF = vector(length=length(h))
sdHF = vector(length=length(h))
for (i in 1:length(h)){
  var[i,1] = quantile(hedge_error[,i], 0.95)
  var[i,2] = quantile(hedge_error[,i], 0.99)
  cvar[i,1] = mean(sort(hedge_error[,i])[seq(0.95*n,n)])
  cvar[i,2] = mean(sort(hedge_error[,i])[seq(0.99*n,n)])
  meanHF[i] = mean(hedge_error[,i])
  sdHF[i] = sd(hedge_error[,i])
}

plass = vector(length=length(h))
sortHF = matrix(nrow=n, ncol=length(h))
sortHF[,1] = sort(hedge_error[,1])
sortHF[,2] = sort(hedge_error[,2])
sortHF[,3] = sort(hedge_error[,3])
sortHF[,4] = sort(hedge_error[,4])

for(i in 1:length(h)){
  for (j in 1:n){
    if(sortHF[j,i]>0){
      plass[i] = j-1
      break
    }
  }
  PG[i] = (n-plass[i])/n
}

postscript("hf1.eps", horizontal=FALSE, height=6, width=8)
par(mfrow=c(2,2))
plot(density(hedge_error[,1]), type="l",
     ylab="Density", xlab="",
     sub="6 trading days", main="")
plot(density(hedge_error[,2]), type="l",
     ylab="Density", xlab="",
     sub="12 trading days", main="")
plot(density(hedge_error[,3]), type="l",
     ylab="Density", xlab="",
     sub="50 trading days", main="")
plot(density(hedge_error[,4]), type="l",
     ylab="Density", xlab="",
     sub="250 trading days", main="")
dev.off()

#Part 4 
#Defined Contribution Pension with Portfolio Rebalancing on Threshold Values

n = 100000 #Number of Monte Carlo simulations
h = 250 #Number of trading days

Pi = 1 #Yearly contribution
g = 0.0275 #Guaranteed interest rate
r = 0.025 #Risk free interest rate
sigma = 0.2 #Yearly stock return volatility
mu = 0.1 #Yearly stock return mean (drift)
lambda = 0.3 #Proportion of stocks in the portfolio

max_variance = c(0, 0.025, 0.05, 0.1)
times_under = matrix(0, n, length(max_variance))
times_over = matrix(0, n, length(max_variance))
times_rebalance = matrix(1, n, length(max_variance))
S = matrix(nrow=n, ncol=h)
hedge_error = matrix(0, n, length(max_variance))
years = vector(length=n)

Returns = matrix(exp(rnorm(h*n, (mu-sigma^2/2)*(1/h),
                       sigma*sqrt(1/h))),n,h)
S[,1] = lambda

for(k in 1:n){
  years[k] = prod(Returns[k,1:h])
}

for (i in 1:length(max_variance)){
  for (j in 1:n){
    limit_upper = lambda+max_variance[i]
    limit_lower = lambda-max_variance[i]
    update_time = 1
    K = exp(g) - (1-lambda)*exp(r)
    a = -pnorm(-(log(lambda/K)+(r+(sigma^2)/2))/sigma)
    b = K*exp(-r)*pnorm(-(log(lambda/K)+(r-(sigma^2)/2))
                        /sigma)
    
    for (t in 2:h){
      S[j,t] = S[j,t-1]*Returns[j,t]
      
      if(S[j,t] > limit_upper){
        afor = a
        bfor = b
        a = -pnorm(-(log(S[j,t]/K)+(r+(sigma^2)/2)*
                       (1-t/h))/(sigma*sqrt(1-t/h)))
        b = K*exp(-r*(1-t/h))*pnorm(-(log(S[j,t]/K) +
                                        (r-(sigma^2)/2)*(1-t/h))/(sigma*sqrt(1-t/h)))
        
        times_over[j,i] = times_over[j,i]+1
        limit_upper = S[j,t] + max_variance[i]
        limit_lower = S[j,t] - max_variance[i]
        
        Vf = afor*S[j,t] + bfor*exp(r*(t-update_time)/h)
        Ve = a*S[j,t] + b
        hedge_error[j,i] = hedge_error[j,i] + (Ve-Vf)*exp(-r*(t/h))
        update_time = t
        
      }else if(S[j,t] < limit_lower){
        afor = a
        bfor = b
        
        a = -pnorm(-(log(S[j,t]/K)+(r+(sigma^2)/2)*
                       (1-t/h))/(sigma*sqrt(1-t/h)))
        b = K*exp(-r*(1-t/h))*pnorm(-(log(S[j,t]/K) +
                                        (r-(sigma^2)/2)*(1-t/h))/(sigma*sqrt(1-t/h)))
        times_under[j,i] = times_under[j,i]+1
        grenseOvre = S[j,t]+max_variance[i]
        grenseNedre = S[j,t]-max_variance[i]
        Vf = afor*S[j,t] + bfor*exp(r*(t-update_time)/h)
        Ve = a*S[j,t] + b
        hedge_error[j,i]=hedge_error[j,i]+
          (Ve-Vf)*exp(-r*(t/h))
        update_time = t
        
      }
    }
    yearly_return = lambda*years[j] + (1-lambda)*exp(r)
    gar = Pi*pmax((exp(g)-yearly_return), 0)
    Vs = a*S[j,t]*Returns[j,t] + b*exp(r/h)
    hedge_error[j,i] = hedge_error[j,i] + (gar-Vs)*exp(-r)
  }
}



times_rebalance = times_rebalance + times_under + times_over
#VaR og CVaR etc. for hedge_errorene:
varHF = matrix(nrow=length(max_variance), ncol=2)
cvarHF = matrix(nrow=length(max_variance), ncol=2)
PG = vector(length=length(max_variance))
meanHF = vector(length=length(max_variance))
sdHF = vector(length=length(max_variance))

for (i in 1:length(max_variance)){
  varHF[i,1] = quantile(hedge_error[,i], 0.95)
  varHF[i,2] = quantile(hedge_error[,i], 0.99)
  cvarHF[i,1] = mean(sort(hedge_error[,i])[seq(0.95*n, n)])
  cvarHF[i,2] = mean(sort(hedge_error[,i])[seq(0.99*n, n)])
  meanHF[i] = mean(hedge_error[,i])
  sdHF[i] = sd(hedge_error[,i])
}

plass = vector(length=length(max_variance))
sortHF = matrix(nrow=n, ncol=length(max_variance))
sortHF[,1] = sort(hedge_error[,1])
sortHF[,2] = sort(hedge_error[,2])
sortHF[,3] = sort(hedge_error[,3])
sortHF[,4] = sort(hedge_error[,4])

for(i in 1:length(max_variance)){
  for (j in 1:n){
    if(sortHF[j,i]>0){
      plass[i] = j-1
      break
    }
  }
  PG[i] = (n-plass[i])/n
}

postscript("hf2.eps", horizontal=FALSE, height=6, width=8)
par(mfrow=c(2,2))
plot(density(hedge_error[,1]), ylab="Density", xlab="",
     sub="Maximal variation: 0%", main="")
plot(density(hedge_error[,2]), ylab="Density", xlab="",
     sub="Maximal variation: 2.5%", main="")
plot(density(hedge_error[,3]), ylab="Density", xlab="",
     sub="Maximal variation: 5%", main="")
plot(density(hedge_error[,4]), ylab="Density", xlab="",
     sub="Maximal variation: 10%", main="")
dev.off()

#Part 5
#Pricing of interest rate guarantee as a knock-out option

put <- function(t, T, g, r, sigma, lambda){
  K = exp(g*(T-t)) - (1-lambda)*exp(r*(T-t))
  d1 = (log(lambda/K) + (r+0.5*(sigma^2))*(T-t))/
    (sigma*sqrt(T-t))
  d2 = d1 - sigma*sqrt(T-t)
  P = K*exp(-r*(T-t))*pnorm(-d2)-lambda*pnorm(-d1)
  P
}

knock <- function(t, T, g, r, sigma, lambda, H){
  K = exp(g*(T-t)) - (1-lambda)*exp(r*(T-t))
  d1 = (log(lambda/K) + (r+0.5*(sigma^2))*(T-t))/
    (sigma*sqrt(T-t))
  d2 = d1 - sigma*sqrt(T-t)
  theta = (r + 0.5*(sigma^2))/(sigma^2)
  y = (log((H^2)/(lambda*K)))/(sigma*sqrt(T-t)) +
    theta*sigma*sqrt(T-t)
  x1 = (log(lambda/H))/(sigma*sqrt(T-t)) +
    theta*sigma*sqrt(T-t)
  y1 = (log(H/lambda)/(sigma*sqrt(T-t))) +
    theta*sigma*sqrt(T-t)
  pPut = put(t,T,g,r,sigma,lambda)
  P = pPut + lambda*pnorm(-x1) -
    K*exp(-r*(T-t))*pnorm(-x1 +
                            sigma*sqrt(T-t)) - sigma*((H/lambda)^(2*theta))*
    (pnorm(y)-pnorm(y1)) +
    K*exp(-r*(T-t))*((H/lambda)^(2*theta-2))*
    (pnorm(y-sigma*sqrt(T-t)) -
       pnorm(y1-sigma*sqrt(T-t)))
  P
}

g = 0.0275 #Guaranteed interest rate
r = 0.025 #Risk free interest rate
sigma = 0.2 #Yearly stock return volatility
mu = 0.1 #Yearly stock return mean (drift)
lambda = 0.3 #Proportion of stocks in the portfolio

Hs = seq(0.001, 0.4, 0.001)

knocks = vector(length=length(Hs))
for(i in 1:length(Hs)){
  knocks[i] = knock(0,1,g,r,sigma,lambda,Hs[i])
}

plot(Hs, knocks, ylab="Price", xlab="H",
     type="l", main="", sub="")


H = 0.20
lambdas = seq(0.05, 1, 0.05)
sigmas = seq(0.05, 0.5, 0.05)
knockSigmaLambda = matrix(nrow=length(lambdas),
                          ncol=length(sigmas))
for (i in 1:length(lambdas)){
  for (j in 1:length(sigmas)){
    knockSigmaLambda[i,j] = knock(0,1,g,r,sigmas[j],
                                  lambdas[i],H)
  }
}

persp(lambdas, sigmas, knockSigmaLambda, theta=-30,
      phi=30, xlab="Stock Proportion", ylab="Volatility",
      zlab="", expand=0.5, ticktype="detailed",
      shade=0.5, ltheta=120)

sigma = 0.2
lambda = 0.3
rs = seq(0.005, 0.1, 0.005)
gs = seq(0.005, 0.1, 0.005)
knockRG = matrix(nrow=length(rs), ncol=length(gs))

for (i in 1:length(rs))
{
  for (j in 1:length(gs))
  {
    knockRG[i,j] = knock(0, 1, gs[j],
                         rs[i], sigma, lambda, H)
  }
}

persp(rs, gs, knockRG, theta=-30, phi=30,
      xlab="Risk Free Interest Rate", ylab="Guaranteed Interest Rate",
      zlab="", expand=0.5, ticktype="detailed", shade=0.5,
      ltheta=120)



#Part 6
#The relationship between the knock-out level and the price reduction

findH<-function(knockReduction,n,t,T,g,r,sigma,lambda){
  pPut = put(t, T, g, r, sigma, lambda)
  Hmin = 0
  Hmax = 1
  for(i in 1:n){
    H = (Hmin+Hmax)/2
    pTest = knock(t, T, g, r, sigma, lambda, H)
    if(pTest>pPut*knockReduction){
      Hmax = H
    }else{
      Hmin = H
    }
  }
  H
}

sigma = 0.2
mu = 0.1
r = 0.025
g = 0.0275
lambda = 0.3
reductions = seq(0,1, 0.01)
Hs = vector(length=length(reductions))

for(i in 1:length(redukctions)){
  Hs[i] = findH(reductions[i], 100, 0, 1, g,
                r, sigma, lambda)
}


plot(reductions, Hs, xlab="Price Reduction",
     ylab="Knock-Out Level",type="l", main="", sub="")

