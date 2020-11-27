library(tseries)     # provide get.hist.quote()
library(timeSeries)  # to provide as.timeSeries()
library(zoo)
library(moments)
library(MASS)
library(fPortfolio) # masks kurtosis and moments from moments
library(copula)     # masks getSigma from fPortfolio
library(fCopulae)
library(mvtnorm)
library(forecast)
library(fGarch)

#########
# 02A EDA
#########

# loading data
P1 = get.hist.quote(instrument = 'T',
                    start="2015-01-01",  end="2015-12-31",
                    quote = c("AdjClose"),provider = "yahoo",
                    compress = "d")

P2 = get.hist.quote(instrument = 'IBM',
                    start="2015-01-01",  end="2015-12-31",
                    quote = c("AdjClose"),provider = "yahoo",
                    compress = "d")

# get daily log returns
r1 = diff(log(P1))
r2 = diff(log(P2))

# calculating Sharpe Ratio
sr1 = mean(r1)/sd(r1)
sr2 = mean(r2)/sd(r2)

# calculating standard error of SR
SE_r1 = sqrt((1+0.5*sr1^2)/length(r1))
SE_r2 = sqrt((1+0.5*sr2^2)/length(r2))

# 95% confidence interval is \hat{SR}\pm SE

# VaR and ES
P = get.hist.quote(instrument = 'T',
                   start="2000-01-01",  end="2012-12-31",
                   quote = c("AdjClose"),provider = "yahoo",
                   compress = "d")

r = diff(log(P))

# VaR_{0.5}
VaR05 = -quantile(r, 0.5)
# VaR_{0.1}
VaR01 = -quantile(r, 0.1)
# VaR_{0.01}
VaR001 = -quantile(r, 0.01)
# ES_{0.01}
ES001 = -mean(r[which(r <= -VaR001)])

# Note: VaR is not coherent because it is does not have subadditivity
# ES is coherent.

#########
# 02B EDA
#########

pSP = get.hist.quote(instrument = '^GSPC',
                     start="2001-01-01",  
                     quote = c("AdjClose"),provider = "yahoo",
                     compress = "d")

timeSP = index(pSP)
rSP = diff(log(pSP))
n = length(rSP)

# Skewness
skewness(rSP)
# skewness test (H0: X is symmetric)
agostino.test(rSP, alternative = "greater")
# negative (left) skewness: longer left tail

# kurtosis
kurtosis(rSP)
# anscombe test (H0: kurtosis = 3)
anscombe.test(rSP)
# if kurtosis > 3, it has heavier tails than normal distribution.

# Jarque Bera test for normality (H0: X follows normal distribution)
jarque.bera.test(rSP)

# Fitting t distribution and QQ plot (from Tutorial 2)
(fit_t = fitdistr(as.numeric(rSP), "t"))
m <- fit_t$estimate[1]
s <- fit_t$estimate[2]
df <- fit_t$estimate[3]
qqplot(qt(ppoints(n), df=df), (as.numeric(rSP)-m)/s,
       xlab = "Quantiles of t distribution",
       ylab = "Quantiles of scaled S&P500 data")
lines(qt(ppoints(n), df=df),
      qt(ppoints(n), df=df))

# Fitting normal distribution
(fit_normal = fitdistr(as.numeric(rSP), "normal"))
mu <- fit_normal$estimate[1]
sd <- fit_normal$estimate[2]
qqplot(qnorm(ppoints(n), mean=mu, sd=sd), (as.numeric(rSP)-mu)/s,
       xlab = "Quantiles of normal distribution",
       ylab = "Quantiles of scaled S&P500 data")
qqline(qnorm(ppoints(100), mean=mu, sd=sd))

# Let X be the fitted distribution. X = (rSP-m)/s follows t_{df}
# VaR at p=0.001
-qt(0.001, df)*s + m
# ES at p=0.001
-mean(rSP[which(rSP <= (qt(0.001, df)*s + m))])

# comparing between fitted distributions
# using AIC/BIC, model with smaller value preferred
BIC(fit_t)
BIC(fit_normal)

#############################
# 03A Multivariate, Portfolio
#############################

MSFTdata = get.hist.quote(instrument = "MSFT", start="2000-01-01", end="2012-12-31",
                          quote = c("AdjClose"),provider = "yahoo", compression = "m")
INTCdata = get.hist.quote(instrument = "INTC", start="2000-01-01", end="2012-12-31",
                          quote = c("AdjClose"),provider = "yahoo", compression = "m")
AAPLdata = get.hist.quote(instrument = "AAPL", start="2000-01-01", end="2012-12-31",
                          quote = c("AdjClose"),provider = "yahoo", compression = "m")
GEdata = get.hist.quote(instrument = "GE", start="2000-01-01", end="2012-12-31",
                        quote = c("AdjClose"),provider = "yahoo", compression = "m")
Tdata = get.hist.quote(instrument = "T", start="2000-01-01", end="2012-12-31",
                       quote = c("AdjClose"),provider = "yahoo", compression = "m")
EBAYdata = get.hist.quote(instrument = "EBAY", start="2000-01-01", end="2012-12-31",
                          quote = c("AdjClose"),provider = "yahoo", compression = "m")
Fdata = get.hist.quote(instrument = "F", start="2000-01-01", end="2012-12-31",
                       quote = c("AdjClose"),provider = "yahoo", compression = "m")
x = merge(MSFTdata, INTCdata, AAPLdata, GEdata, Tdata, EBAYdata, Fdata)
R = diff(log(x))
colnames(R)[1] = "MSFT"
colnames(R)[2] = "INTC"
colnames(R)[3] = "AAPL"
colnames(R)[4] = "GE"
colnames(R)[5] = "ATT"
colnames(R)[6] = "EBAY"
colnames(R)[7] = "F"

pairs(R)
cor(R)

R = timeSeries(R)

# other estimators for covariance matrix from fPortfolio
covEstimator(R)
shrinkEstimator(R)
# shrink estimator usually have smaller estimation error
# than sample covariance matrix

#############################
# 03B Multivariate, Portfolio
#############################

# Plot efficient portfolio between 2 assets (Tutorial 2)
# note that the Cov matrix here suggests no correlation between them
cov <- matrix(c(1, 0, 0, 9), 2, 2)
mu <- matrix(c(0.1, 0.2), 2)
r <- matrix(0, 10000, 1)
s <- matrix(0, 10000, 1)

for (i in 1:10000){
  x <- runif(1)
  w <- matrix(c(x, 1-x), 2, 1)
  r[i] <- t(w) %*% mu
  s[i] <- sqrt(t(w) %*% cov %*% w)
}
plot(s, r, cex = 0.5, xlab=expression(sigma(w)),
     ylab="r(w)", main=bquote(rho~"= 0"))

# Specifying portfolios using fPortfolio using 7 stocks above
Spec = portfolioSpec(portfolio=list(targetReturn = 0.003, 
                                    nFrontierPoints = 500, 
                                    riskFreeRate = 0.001)) 
Constraints = "Longonly"

efficientPortfolio(R, Spec, Constraints)
tangencyPortfolio(R, Spec, Constraints)

# Outputs:
# Cov refers to SD / risk
# CVar refers to ES.
# VaR is the value-at-risk at probability 0.05

# weights can be extracted using getPortfolio(portfolio)$weights

# portfolio returns for tangency portfolio
RN = R %*% getPortfolio(tangencyPortfolio(R, Spec, Constraints))$weights

# plotting tangency portfolio returns
plot(RN, type="l")

# expected/average return
mean(RN)

# standard dev of return
sd(RN)

# VaR
-quantile(RN, 0.05)

# min risk portfolio
minriskPortfolio(R)

# visualizing the efficient frontier and portfolios (with selections)
out = portfolioFrontier(R, spec=Spec, constraints = "Longonly")
plot(out)

#############################
# 03C Multivariate, Portfolio
#############################

# Finding efficient portfolios (Tutorial 3 Q2)

# min risk portfolio
Sigma = matrix(c(0.0210,0.0095, 0.0131 , 0.0119,
                  0.0095 , 0.0175 , 0.0091 , 0.0092,
                  0.0131 , 0.0091 , 0.0230 , 0.0095,
                  0.0119 , 0.0092 , 0.0095 , 0.0189), 4, 4)
mu = c(0.0090, 0.0125, 0.0070, 0.0111)
ones = c(1,1,1, 1)
Sigmainv = solve(Sigma)
top = Sigmainv %*% ones
bot = t(ones) %*% Sigmainv %*% ones
(w = top/as.numeric(bot))

# set target return = 0.0111, find efficient portfolio
target.mu = 0.0111
a = as.numeric( t(ones) %*% solve(Sigma) %*% ones )
b = as.numeric( t(ones) %*% solve(Sigma) %*%mu )
c = as.numeric( t(mu) %*% solve(Sigma) %*%mu )
part.1 = (c-b*target.mu)/(a*c-b*b) * solve(Sigma) %*% ones
part.2 = (a*target.mu-b)/(a*c-b*b) * solve(Sigma) %*% mu
(w = part.1+part.2)

# tangency portfolio, assume risk free is 0.001
top.tangency = solve(Sigma) %*% (mu - 0.001 * ones)
bot.tangency = as.numeric( t(ones) %*% top.tangency)
(w.tangency = top.tangency/bot.tangency)

############
# 04A Copula
############

# fitting distribution and copula (Tutorial 4 Q7)
MMM <- get.hist.quote(instrument="MMM", start="2008-01-01", end="2009-12-31",
                      quote="AdjClose", provider="yahoo", compression="d")
MSFT <- get.hist.quote(instrument="MSFT", start="2008-01-01", end="2009-12-31",
                       quote="AdjClose", provider="yahoo", compression="d")
# log return
r_MMM <- as.matrix(diff(log(MMM)))
r_MSFT <- as.matrix(diff(log(MSFT)))
# fitting t-distribution
(fit_t_MMM <- fitdistr(r_MMM, "t"))
(fit_t_MSFT <- fitdistr(r_MSFT, "t"))

# Transforming data using rank
n <- nrow(r_MMM)
U_MMM <- rank(r_MMM)/(n+1)
U_MSFT <- rank(r_MSFT)/(n+1)
data <- data.frame(list(U_MMM=U_MMM, U_MSFT=U_MSFT))
data <- data.matrix(data)
# fit t-copula
tCopula <- tCopula(dim=2)
(fit_cop <- fitCopula(tCopula, data))

# to specify the full (complicated) joint distribution, we need to plug in
# the m, s from t-fit and rho, df from t-copula-fit.

# Monte-Carlo Simulation (Tutorial 4 Q8)
N <- 1000000
r1 <- rnorm(N, mean = 0.01, sd = sqrt(0.1))
r2 <- rnorm(N, mean = 0.02, sd = sqrt(0.2))
joint <- 0.4*r1 + 0.6*r2
(VaR <- -quantile(joint, 0.01))
(ES <- -mean(joint[joint < -VaR]))

# jointly normal with rho=0.6 (use normal copula)
# use copula samples to generate samples u=[0,1]
# then use qnorm to get x values
ncop <- normalCopula(param=0.6, dim=2)
rC <- rCopula(N, ncop)
x1 <- qnorm(rC[,1], mean=0.01, sd=sqrt(0.1))
x2 <- qnorm(rC[,2], mean=0.02, sd=sqrt(0.2))
portf <- 0.4*x1 + 0.6*x2
(Var <- -quantile(portf, 0.01))
(ES <- -mean(portf[portf < -Var]))

# t-copula with df=2, rho=0.6
tcop <- tCopula(param=0.6, df=2)
rCt <- rCopula(N, tcop)
x1 <- qnorm(rCt[,1], mean=0.01, sd=sqrt(0.1))
x2 <- qnorm(rCt[,2], mean=0.02, sd=sqrt(0.2))
portf <- 0.4*x1 + 0.6*x2
(VaR <- -quantile(portf, 0.01))
(ES <- -mean(portf[portf < -VaR]))

#########################
# 05A CAPM, Factor Models
#########################

# Single Index Model

r = read.table('./FS05A/6_Portfolios_2x3_daily.txt')
r = data.frame(r[,2:7])

colnames(r)[1] = "R1"
colnames(r)[2] = "R2"
colnames(r)[3] = "R3"
colnames(r)[4] = "R4"
colnames(r)[5] = "R5"
colnames(r)[6] = "R6"

Factor = read.table('./FS05A/FF_Research_Data_Factors_daily.txt')
F = data.frame(Factor[,2:5])

colnames(F)[1] = "Mkt"
colnames(F)[2] = "SMB"
colnames(F)[3] = "HML"
colnames(F)[4] = "rf"

reg1 = lm(r$R1 ~ F$Mkt ) 
reg2 = lm(r$R2 ~ F$Mkt ) 
reg3 = lm(r$R3 ~ F$Mkt ) 
reg4 = lm(r$R4 ~ F$Mkt ) 
reg5 = lm(r$R5 ~ F$Mkt ) 
reg6 = lm(r$R6 ~ F$Mkt ) 

summary(reg1)
summary(reg2)
summary(reg3)
summary(reg4)
summary(reg5)
summary(reg6)

# CAPM

xDJI = get.hist.quote(instrument = "^DJI", start="2000-01-01", end="2013-02-11",
                      quote = c("AdjClose"),provider = "yahoo", compression = "d")


xi = get.hist.quote(instrument = "T", start="2000-01-01", end="2013-02-11",
                    quote = c("AdjClose"),provider = "yahoo", compression = "d")

x = merge(xi, xDJI)
R = diff(log(x))
FactorS = read.csv('./FS05A/Factors.csv')
date = as.Date(strptime(FactorS$date, "%Y%m%d"))
rF = FactorS$rF/30;    # the risk-free return is monthly
rF = zoo(rF, date)

R = na.omit(merge(R, rF))
# regressing excess return of stock against excess market return
RmRf = R[,2]-R[,3]    # use DJI as market
Ri = R[,1]-R[,3]      # AT&T is our stock
summary(lm(Ri~RmRf))

#########################
# 05B CAPM, Factor Models
#########################

# Fama-French 3 Factor Model

r = read.table('./FS05B/6_Portfolios_2x3_daily.txt')
r = data.frame(r[,2:7])

colnames(r)[1] = "R1"
colnames(r)[2] = "R2"
colnames(r)[3] = "R3"
colnames(r)[4] = "R4"
colnames(r)[5] = "R5"
colnames(r)[6] = "R6"

Factor = read.table('./FS05B/FF_Research_Data_Factors_daily.txt')
F = data.frame(Factor[,2:5])

colnames(F)[1] = "Mkt_rf"
colnames(F)[2] = "SMB"
colnames(F)[3] = "HML"
colnames(F)[4] = "rf"

reg1 = lm(r$R1-F$rf ~ F$Mkt_rf + F$SMB + F$HML) 
reg2 = lm(r$R2-F$rf ~ F$Mkt_rf + F$SMB + F$HML) 
reg3 = lm(r$R3-F$rf ~ F$Mkt_rf + F$SMB + F$HML) 
reg4 = lm(r$R4-F$rf ~ F$Mkt_rf + F$SMB + F$HML) 
reg5 = lm(r$R5-F$rf ~ F$Mkt_rf + F$SMB + F$HML) 
reg6 = lm(r$R6-F$rf ~ F$Mkt_rf + F$SMB + F$HML) 

summary(reg1)
summary(reg2)
summary(reg3)
summary(reg4)
summary(reg5)
summary(reg6)
# regression of return on factors are significant at level 0.01.
# R^2 all above 90%

#########################
# 05C CAPM, Factor Models
#########################

# Principal Component Analysis

MSFTdata = get.hist.quote(instrument = "MSFT", start="2000-01-01", end="2012-12-31",
                          quote = c("Adjusted"), provider = "yahoo", compression = "m")
MSFT = diff(data.matrix(log(MSFTdata$Adjusted)))

INTCdata = get.hist.quote(instrument = "INTC", start="2000-01-01", end="2012-12-31",
                          quote = c("Adjusted"), provider = "yahoo", compression = "m")
INTC = diff(data.matrix(log(INTCdata$Adjusted)))

GEdata = get.hist.quote(instrument = "GE", start="2000-01-01", end="2012-12-31",
                        quote = c("Adjusted"),provider = "yahoo", compression = "m")
GE  = diff(data.matrix(log(GEdata$Adjusted)))


YHOOdata = get.hist.quote(instrument = "BA", start="2000-01-01", end="2012-12-31",
                          quote = c("Adjusted"),provider = "yahoo", compression = "m")
YHOO = diff(data.matrix(log(YHOOdata$Adjusted)))

mydata = data.frame(list(MSFT=MSFT, INTC=INTC, GE=GE, YHOO=YHOO))
mydata = as.timeSeries(mydata)

colnames(mydata)[1] = "MSFT"
colnames(mydata)[2] = "INTC"
colnames(mydata)[3] = "GE"
colnames(mydata)[4] = "YHOO"

(mypca = princomp(mydata, cor = TRUE))

summary(mypca)

loadings(mypca)   # note that blank entries are small but not zero
plot(mypca)       # shows a screeplot.

mydataPCA = merge(mydata, mypca$scores)
pairs( mydataPCA[,1:6])

# to find the eigenvectors and eigenvalues
(S = cov(mydata)) # or S = cor(mydata)
(eigen(S)$vector)
(eigen(S)$value)

# using 6-portfolio data

r = read.table('./FS05C/6_Portfolios_2x3_daily.txt')
r = data.matrix(r[,2:7])

mypca = princomp(r, cor = TRUE)
summary(mypca)

loadings(mypca)  
plot(mypca)

PCA = data.matrix(mypca$scores[,1:3]) # first 3 components

mydataPCA = data.frame(matrix(c(r, PCA), length(r[,1]), 9))

colnames(mydataPCA)[1] = "R1"
colnames(mydataPCA)[2] = "R2"
colnames(mydataPCA)[3] = "R3"
colnames(mydataPCA)[4] = "R4"
colnames(mydataPCA)[5] = "R5"
colnames(mydataPCA)[6] = "R6"
colnames(mydataPCA)[7] = "PC1"
colnames(mydataPCA)[8] = "PC2"
colnames(mydataPCA)[9] = "PC3"

pairs(mydataPCA)

## Tutorial 5 Q6
# regression model for 30 stocks

name0 = ""
name0[ 1 ]="MMM"
name0[ 2 ]="AXP"
name0[ 3 ]="AAPL"
name0[ 4 ]="BA"
name0[ 5 ]="CAT"
name0[ 6 ]="CVX"
name0[ 7 ]="CSCO"
name0[ 8 ]="KO"
name0[ 9 ]="DIS"
name0[ 10 ]="DOW"
name0[ 11 ]="XOM"
name0[ 12 ]="GS"
name0[ 13 ]="HD"
name0[ 14 ]="IBM"
name0[ 15 ]="INTC"
name0[ 16 ]="JNJ"
name0[ 17 ]="JPM"
name0[ 18 ]="MCD"
name0[ 19 ]="MRK"
name0[ 20 ]="MSFT"
name0[ 21 ]="NKE"
name0[ 22 ]="PFE"
name0[ 23 ]="PG"
name0[ 24 ]="TRV"
name0[ 25 ]="UTX"
name0[ 26 ]="UNH"
name0[ 27 ]="VZ"
name0[ 28 ]="V"
name0[ 29 ]="WMT"
name0[ 30 ]="WBA"

x = get.hist.quote(instrument = "^DJI", start="2000-01-01",
                   quote = c("AdjClose"),provider = "yahoo", compression = "m")
for (i in 1:30){
  xi = get.hist.quote(instrument = name0[i], start="2000-01-01",
                      quote = c("AdjClose"),provider = "yahoo", compression = "m")
  x = merge(x, xi)
}
R = diff(log(x))
R = na.omit(R)

# matrices to store alpha, beta and returns
a = matrix(0, 30, 1)
b = matrix(0, 30, 1)
r = matrix(0, 30, 1)

for (i in 1:30)
{
  abi = lm(R[,i+1]~R[,1])
  r[i] = mean(R[,i+1])
  a[i] = abi$coefficients[1]
  b[i] = abi$coefficients[2]
  print(summary(abi))
}

# regressing alpha against beta
ab = lm(a~b)
summary(ab)

# fitting security market line
SML = lm(r~b)
summary(SML)

# plotting
plot(b, r, ylab="return")
abline(SML)
text(b+0.03, r, labels = name0, cex=0.7, offset=0.5)
# recommend the stocks above the SML (underpriced)

#################
# 06A Time Series
#################

P = get.hist.quote(instrument = "MSFT", start="2010-01-01", end="2013-4-1",
                   quote = c("AdjClose"),provider = "yahoo", compression = "d")
r = as.timeSeries(diff(log(P)))

# to see autocorrelation (ACF)
acf(r)
# the region bounded by 2 horizontal blue lines is the acceptance region
# for the hypothesis H0: ACF(h) = 0
# for any lag h > 0 outside of the region, we can say that ACF is not 0
# a white noise has all ACF(h) = 0, where h > 0

# using Ljung-Box test to test for white noise (output p-value)
Box.test(r, type = "Ljung")

# test for stationarity: ADF test
# H0: time series has unit root - not stationary
# H1: time series does not have unit root - stationary
adf.test(r, alternative = "stationary")

#################
# 06B Time Series
#################

P = get.hist.quote(instrument = "DJIA", start="2008-01-01", end="2013-01-1",
                   quote = c("AdjClose"),provider = "yahoo", compression = "d")
plot(P)
r = diff(data.matrix(log(P)))

# model fitting
mARIMA = arima(r, order = c(2,0,0))

# check diagnostics
tsdiag(mARIMA)
# one of the output is the p-values for Ljung-Box statistic
# if all the p-values are high, the residuals are white noise,
# which means the model fits the data well.
# this is necessary to check for fitting models, but when comparing models,
# also compare AIC - output when calling arima(), or use AIC(model).

# make predictions
predict(mARIMA, n.ahead = 5)

############
# ARCH/GARCH
############

# Fit ARMA(p,q) + GARCH(m,s) model
# cond.dist = distribution of standardized residuals
# shape = if t-dist, this is the df

# t-dist for standardized residuals
mGARCHt = garchFit(~arma(1,1)+garch(2,1),data=r*100, 
                   cond.dist="std", shape=7, include.shape=FALSE, trace=FALSE)
summary(mGARCHt)
plot(mGARCHt)

# t-dist for standardized residuals but auto select df
mGARCHtA = garchFit(~arma(1,1)+garch(2,1),data=r*100, 
                    cond.dist="std", trace=FALSE)

# fitting APARCH with normal standardized residuals
mAPARCHn = garchFit(~aparch(1,1),data=r*100,trace=F, 
                    delta=2, include.delta=F,
                    cond.dist="norm")
summary(mGARCHtA)
plot(mGARCHt)

# a good model is when both residuals and standardized residuals are white noise

# make predictions for n days ahead
predict(mGARCHt, n.ahead = 5)
# meanForecast: predicted returns
# meanError: predicted volatility
# standardDeviation: SD of predicted returns

# Appendix:
# [object]@residuals: (raw, unstandardized) residual values 
# [object]@fitted: the fitted values
# [object]@h.t: the conditional variances (h.t = sigma.t^2)
# [object]@sigma.t: the conditional standard deviation.
# [object]@fit$coef: estimated parameters in the model
# [object]@fit$se.coef: SE of estimated parameters
































