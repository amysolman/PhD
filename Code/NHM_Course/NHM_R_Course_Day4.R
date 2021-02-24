# R4ALL Day 4

#visualize your data before doing further
#statistical analysis

#General linear models in R
#lm()

#Special cases of the same model:
#ANOVA ("Analysis of Variance") - all categorical predictors
#Simple regression - only one numeric predictor
#ANCOVA - several numeric predictors
#Multiple Regression - categorical and numeric predictors

#What are they trying to do?
#One response variable (dependent variable)
#One or more predictor variables (independent variables)
#Numeric variables and/or categorical variables

#Assumptions
#What are the key assumptions?
#Independence (one observation doesn't tell you anything about another observation)
#Measurement scale (y variable is a numeric quantity)
#Constant variance (variability left over in data after model 
#fitting doesn't change e.g. get bigger for bigger values)
#Normality
#Is the functional form correct? (Is the shape of the relationship 
#between x and y appropriate?)

#Most assumptions relate to the residuals (what's left over after 
#the model fitting)

#Distributional assumptions (e.g. normality and equal variance) 
#are NOT so important unless your focus is prediction
#Getting the model structure is important
#Is the additivity assumption met?
#Have you captured the interactions?
#If there is a curvy relationship you need a model to capture 
#a curvy relationship (e.g. linear vs logistic model?)
#The most important thing you need to do is make sure the
#model is adequent for the question you have

#Simple Linear Regression

#y=mx+c
#two coefficients define the relationship
#c = intercept and m = slop/gradient

#In R formula speak y ~ x

#New example ----
#Plant growth Data

library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggfortify)

growth<- read_csv("Data/plant.growth.rate.csv")

ggplot(growth, aes(x=soil.moisture.content, y=plant.growth.rate))+
  geom_point(colour="cornflowerblue", size=5)+
  labs(x = expression(paste("Soil Moisture", bar(x) / mu)),y = "Growth Rate (mm/week)")+
  theme_bw()

#Guess the slope and intercept
#Intercept = 17
#Slope = 21

#a = 0.25, 23
#b = 1.75, 49
#(49-23)/(1.75-0.25 )
#y = (17*1) + 17

#Model fitting

mod1_growth <- lm(plant.growth.rate ~ soil.moisture.content, data = growth)

#lots of stuff in here... not what is captured
names(mod1_growth)

#check assumptions using model object
#check the diagnostics

#residuals - what's left over?
#big vs small magnitude, positive vs negative

autoplot(mod1_growth)

#What are these plots for?
#Assess whether the assumed mean-variance relationship is appropriate - A) Scale-Location
#Evaluate whether the shape of the relationship we modelled is appropriate - B) Residuals vs fitted
#Identify data points that may have a large impact on the fitted model - C) Residual vs Leverage
#Evaluate the shape of the distribution we used is appropriate - D) Normal Q-Q

#A) Scale-Location
#B) Residuals vs fitted - yes!
#C) Residual vs Leverage
#D) Normal Q-Q

#Evaluate model conclusions

#The workhorses
#Understanding what your fitted model is telling you
anova(mod1_growth) #<testing overall significance
#global significance of terms - a.k.a. main effects and interactions
#F ratio tests
#sequential sums of squares if unbalanced data

summary(mod1_growth) # about individual coefficients 
#estimates of coefficients and their contrasts (differences)
#significance of these coefficients/contrasts
#t-tests

#is there a significant effect of soil moisture on plant growth?
Yes
#do you understand what is shown in the summary table?
#do you understand what is shown in the ANOVA table?
#F value
#p value  - small number = highly significant effect (calculated from F statistics and degrees of freedom (1,48))
#DF = 1 = 1 slope
#Residuals DF = number of data points minus 2

#summary - you want to look at the coefficients table
#need to know coefficients of your model if you want to predict
#dependent variable with different independent variable values
#t-value calculated from coefficient/standard error
#from t statistic p value is calculated

#regression - plotting with prediction and intervals
#let's add confidence intervals to our plot

#Using predict
x_vals <- seq(from=0.25,to=2,length=50)
x_vals
newX <- expand.grid(soil.moisture.content = x_vals)
newX

#generate predictions
newY <- predict(mod1_growth, newdata = newX, interval='confidence')
newY

#bind the dataframes together
newXY <- cbind(newX, newY)
addThese <- data.frame(newX, newY)

addThese <- rename(addThese, plant.growth.rate = fit)
addThese

#add confidence intervals to graph
ggplot(growth, aes(x=soil.moisture.content, y=plant.growth.rate))+
  geom_point(colour="cornflowerblue", size=5)+
  geom_smooth(data=addThese,
              aes(ymin=lwr, ymax=upr),
              stat='identity')+ #identity means just plot what I give you don't mess about with the data
  labs(x = expression(paste("Soil Moisture", bar(x) / mu)),y = "Growth Rate (mm/week)")+
  theme_bw()

# ANCOVA ----

limpets <- read_csv("Data/limpet.csv")

ggplot(limpets, aes(x=DENSITY, y=EGGS, colour=SEASON))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE)+
  theme_bw()

#ANCOVA model with interaction
limpet_mod <- lm(EGGS ~ DENSITY * SEASON, data = limpets)

#Check the assumptions
autoplot(limpet_mod, smooth.colour = NA)

#examine the results
anova(limpet_mod) #ANOVA tables with sequential sums of squares
#anova says density has an effect, season has an effect, there is no significant interaction between density and season
#so we need to do a different type of model?
limpet_mod2 <- lm(EGGS ~ DENSITY + SEASON, data = limpets)

summary(limpet_mod) #coefficients

