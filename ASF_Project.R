rm(list = ls())
graphics.off()

setwd("/Users/manenti_paolo/Desktop/FINASS/Advanced Statistic/PROJECT/")

library(car)
library(dplyr)
library(sf)
library(magrittr)

data <- read.table("datBI.txt", sep = ";")
dataset <- data.frame(data$HB0900, data$HB0100, data$HB0200, data$HB4400, data$HD1110, data$HG0410, data$HI0100, data$HNB0920)
summary(dataset) 

dataset1 <- na.omit(dataset)
rm(data, dataset)
Dataset <- dataset1
rm(dataset1)

rownames(Dataset) <- 1:nrow(Dataset)
dataset <- Dataset
dataset$Observations <- 1:nrow(dataset)

rm(Dataset)
colnames(dataset) <- c("Y, Current $ main residence", "X1, SQ House", "X2, YRS in household",
                      "X3, Tot. Val. Cars", "X4, Sight Account", "X5, Gross Income from fin. invest.", 
                      "X6, Food expences", "X7, HMR/Imputed rent", "Observations")

# VARIABLES OF INTEREST
Y <- dataset$`Y, Current $ main residence` #current price main residence (in the original dataset corresponds to FHB0900)
X1 <- dataset$`X1, SQ House` #SQ houses (in the original dataset corresponds to FHB0100) 
X2 <- dataset$`X2, YRS in household` #years spent in the household (in the original dataset corresponds to FHB0200)
X3 <- dataset$`X3, Tot. Val. Cars` #total value of cars (in the original dataset corresponds to FHB4400)
X4 <- dataset$`X4, Sight Account` #sight accounts (in the original dataset corresponds to FHD1110)
X5 <- dataset$`X5, Gross Income from fin. invest.` #gross income from financial investments (in the original dataset corresponds to FHG0410)
X6 <- dataset$`X6, Food expences` #food expences (in the original dataset corresponds to FHI0100)
X7 <- dataset$`X7, HMR/Imputed rent` #HMR/imputed rent (in the original dataset corresponds to FHNB0920)

summary(dataset)

####################################
###### MODEL SPECIFICATIONS ########
####################################

pairs(dataset, upper.panel = panel.smooth)
GGally::ggpairs(dataset[, c(1:8)], axisLabels = "none") #Another way to represent the relation among regressors and the respective correlation

model1 <- lm(log(Y) ~ log(X1) + X2 + log(X3) + X4 + log(X5) + log(X6) + log(X7), data = dataset, x = T, y = T)
summary(model1)

# we eliminated X2 and log(X6) because of their high p-values 
model2 <- update(model1, .~. -X2 -log(X6))
summary(model2)

MASS::stepAIC(model2) # As a confirmation, the model4 seems to be the best and we eliminate the others from the environment

res <- model2$residuals

# Plots of residuals
par(mfrow = c(1,3))
plot(res, ylab = "Raw residuals", xlab = "Observation index", main = "Raw residuals")
abline(h = mean(res), col = "red")
plot(rstandard(model2), ylab = "Standardized residuals", xlab = "Observation index", main = "Plot of standardized residuals")
abline(h = mean(res), col = "red")
plot(rstudent(model2), ylab = "Studentized residuals", xlab = "Observation index", main = "Plot of studentized residuals")
abline(h = mean(res), col = "red")
# The distribution of the point clouds seems to be near to 0, thus the linearity assumption is valid
# There is a plausible presence of heteroskedasticity

par(mfrow = c(1,1))

# Residuals vs each regressors
head(model2$x)
par(mfrow=c(3,2))
for (i in c(2:6)){
  plot(model2$x[,i], rstandard(model2), ylab='Standardized Residuals', xlab=colnames(model2$x)[i],main=paste('Standardized residuals vs',colnames(model2$x)[i], sep=' '))
}
# As a further confirmation, also the plots of residuals vs each regressors show the confirmation of the linearity hypothesis, and the plausible presence of heteroskedasticity and outliers 

# Fitted values and scatterplots
par(mfrow = c(1,2))
plot(model2$fitted.values, model2$residuals, ylab = "Residuals", xlab = "Fitted values", main = "Plot of fitted values and residuals")
abline(h = mean(res), col = 2)
plot(model2$fitted.values, rstandard(model2), ylab = "Standardized residuals", xlab = "Fitted values", main = "Plot of standardized fitted values and residuals")
abline(h = mean(res), col = 2)

par(mfrow = c(1,1))
# Squared residuals
sq_res <- res^2 
# y hat
plot(model2$fitted.values, sq_res, xlab = expression(paste("Estimated", sep = " ", hat("y"))), ylab = "Residuals")
#As further confirmation of the plausible presence of heteroskedascitity we compute the squared residuals

# Squared residuals vs each regressors 
par(mfrow = c(3,2))
for (i in c(2:6)){
  plot(model2$x[, i], sq_res, ylab = expression(paste("Residuals"^2)), xlab = colnames(model2$x)[i])
}
par(mfrow = c(1,1))

plot(model2) #The scale location and the Cook's distance plots help us to identify the presence of outliers

# We can obtain Leverages by the function hat (it stands for hat values) applied to the model.matrix:
lev <- hat(model.matrix(model2))
plot(lev, ylab = "Leverages", main = "Index plot of Leverages") 
lev.t <- 2*ncol(model.matrix(model2))/nrow(model.matrix(model2)) # threshold leverage, twice the average level
abline(h = lev.t, col = "red")
# Observations with high leverage:
h.l <- cbind(which(lev > lev.t), lev[c(which(lev > lev.t))]) # there are 269 cases

#Cook's distance
cooksD <- cooks.distance(model2)
plot(cooksD)
abline(h = 4/length(cooksD), col = "green")
d.inf <- cooksD <= 4/length(cooksD)
table(d.inf) # there are 231 influential observations

# Check for NORMALITY ASSUMPTION
qqnorm(model2$residuals, main = "Normal Q-Q Plot of residuals")
qqline(model2$residuals, col = "red")
# The normality assumption is not satisfied given the heavy tails in the representation of the distribution
# We need to compare the normal qqplot with the student-t one to point out the presence of outliers
car::qqPlot(model2, ylab = "", main = "Studentized residuals Q-Q plot")

# Tests for the normality assumption
tseries::jarque.bera.test(resid(model2)) # Ho: Normality
# since the p-value is less than alpha, we have to reject Ho
# non-parametric test
ks.test(model2$residuals, "pnorm") # Ho: Normality (REJECTED)

# The Shapiro - Wilk test is not performed due to the dimension of the dataset
# As can be seen form the plots and also from the tests, the residuals are NOT normally distributed, 
# thus we proceed with the elimination of the influential observations

# let's try to see how many of the influential observations are outliers:
print(outl <- car::outlierTest(model2, labels = dataset$Observations)) # Bonferroni outlier Test
# For the interpretation of the results, the observations characterized by a p-value < 0.05 
# are highlighted by the test.

d.outl <- dataset$Observations == names(outl$rstudent)[1] | dataset$Observations == names(outl$rstudent)[2]|
  dataset$Observations == names(outl$rstudent)[3] | dataset$Observations == names(outl$rstudent)[4]|
  dataset$Observations == names(outl$rstudent)[5] | dataset$Observations == names(outl$rstudent)[6]|
  dataset$Observations == names(outl$rstudent)[7] | dataset$Observations == names(outl$rstudent)[8]|
  dataset$Observations == names(outl$rstudent)[9] | dataset$Observations == names(outl$rstudent)[10]

# Before excluding the outliers found, we elaborate a different version of the model 
# that includes the dummy variable d.outl that might be a solution for this issue

model3 <- update(model2, .~. + d.outl)
tseries::jarque.bera.test(resid(model3))
ks.test(model3$residuals, "pnorm") 

model4 <- update(model2, .~., subset = d.outl == 0)
summary(model4) 
summary(model3) # best model which contains the outliers
summary(model2)

# there is not the normality yet
tseries::jarque.bera.test(resid(model4))
shapiro.test(resid(model4))
ks.test(model4$residuals, "pnorm")

# we can try to exclude all the influential observations
model5 <- update(model4, subset = d.inf)
summary(model5)

tseries::jarque.bera.test(resid(model5))

# After obtaining the normality of the residuals we can proceed with the Mis-specifications analysis 

# UNDERFITTING

lmtest:: dwtest(model5, order.by = dataset$`X2, YRS in household`[cooksD<(4/length(cooksD))])
lmtest:: dwtest(model5, order.by = dataset$`X6, Food expences`[cooksD<(4/length(cooksD))])

#  One out of two regressors that has been excluded from the analysis seem not to be relevant in the analysis
# the Durbin - Watson test, in fact, shows that the residuals are not correlated. Only X6 needs to be readmitted 
# H0: autocorrelation of the disturbances is 0 (ACCEPTED)

model6 <- update(model5, .~. + I(X6))
summary(model6)
tseries::jarque.bera.test(model6$residuals)

# OVERFITTING

vif(model6) #Considering that values under 2 are still considered acceptable, we can be confident that the problem of multicollinearity is not present

View(cor(model.matrix(model6)[, c(2:7)]))
# Considering the correlation matrix, the absence of multicollinearity is confirmed because of the low values of the correlation among the regressors. 

# LINEARITY

lmtest::resettest(model6, power = 2:3, type = "fitted")
# The RESET test shows that none of the relations has to be changed (p-value > 0.05) 

res2 <- model6$residuals

par(mfrow = c(1,3))
plot(res2, ylab = "Raw residuals", xlab = "Observation index", main = "Raw residuals")
abline(h = mean(res2), col = "red")
plot(rstandard(model6), ylab = "Standardized residuals", xlab = "Observation index", main = "Plot of standardized residuals")
abline(h = mean(res2), col = "red")
plot(rstudent(model6), ylab = "Studentized residuals", xlab = "Observation index", main = "Plot of studentized residuals")
abline(h = mean(res2), col = "red")
# The mean of the residuals is concentrated in 0 and there may be presence of heteroskedasticity 

# Residuals vs each regressors
head(model6$x)
par(mfrow=c(3,2))
for (i in c(2:7)){
  plot(model6$x[,i], rstandard(model6), ylab='Standardized Residuals', xlab=colnames(model6$x)[i],main=paste('Standardized residuals vs',colnames(model6$x)[i], sep=' '))
}

# Fitted values and scatterplots
par(mfrow = c(1,2))
plot(model6$fitted.values, model6$residuals, ylab = "Residuals", xlab = "Fitted values", main = "Plot of fitted values and residuals")
abline(h = mean(res2), col = 2)
plot(model6$fitted.values, rstandard(model6), ylab = "Standardized residuals", xlab = "Fitted values", main = "Plot of standardized fitted values and residuals")
abline(h = mean(res2), col = 2)

par(mfrow = c(1,1))

# Squared residuals
sq_res2 <- res2^2 

# Squared residuals vs each regressors 
par(mfrow = c(3,2))
for (i in c(2:7)){
  plot(model6$x[, i], sq_res2, ylab = expression(paste("Residuals"^2)), xlab = colnames(model6$x)[i])
}

par(mfrow = c(1,1))

plot(model6) 
plot(model6, which = 4, col = "red") #We can see a smaller number of influential observations

#Check for normality
par(mfrow=c(1,2))

qqnorm(model6$residuals, main = "Normal Q-Q Plot of residuals")
qqline(model6$residuals, col = "red")

hist(model6$residuals, main = "Histogram of residuals", xlab = "Residuals", freq = F, breaks = 10)
xfit <- seq(min(model6$residuals), max(model6$residuals), length = 1000)
yfit <- dnorm(xfit, mean = mean(model6$residuals), sd = sd(model6$residuals)) #normal distribution
lines(xfit, yfit, type = "l", col = "blue", lwd = 1.5) 

car::qqPlot(model6, ylab = "", main = "Studentized residuals Q-Q plot")

# The normality is now obtained just by looking at the qqplot, that shows 
# correspondence among the theoretical distribution and the realized one 

# The normality is also confirmed by the tests
tseries::jarque.bera.test(resid(model6)) 
ks.test(model6$residuals, "pnorm")

olsrr::ols_plot_cooksd_bar(model6)

# HOMOSKEDASTICITY TEST

# Towards predicted y values
plot(sq_res2 ~ fitted(model6), xlab = expression(paste("fitted values", sep = " ", hat("y"))), ylab = expression(paste("squared residuals", sep = " ", hat("e"))))

# Towards each regressor
x.m6 <- as.data.frame(model.matrix(model6))

par(mfrow = c(3, 2))
for (i in 2:7)  {
  plot(sq_res2 ~ x.m6[, i], ylab = expression(paste("Residuals", sep = " ", hat("e")^2)), xlab = colnames(x.m6)[i])
}

# Breusch-Pagan test

lmtest::bptest(model6, studentize = T) 
lmtest::bptest(model6, studentize = F)
# The residuals are then heteroskedastic (p-value < 0.05)

m.VH = sandwich::vcovHC(model6) 
# Then, to perform test procedure with alternative estimates of variance and covariance matrix:

model6.robust <- lmtest::coeftest(model6, vcov = sandwich::vcovHC(model6, type="HC0"), df=3938) # White correction

summary(model6) #No significant difference in the Std. Error among model 6 and the model6.robust. 
# Thus, model6.robust has also improved in significance, with an higher p-value in log(X3) and X4
model6.robust

# Model6.robust is thus the best. 

