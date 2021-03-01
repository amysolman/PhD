#Three stats in three hours
#Andrew Beckerman
#20.01.21

#Setup ---- <- this is good for segmenting your code!

install.packages("agricolae")
install.packages("ggfortify")

library(agricolae)
library(ggfortify)
library(ggplot2)
library(tidyverse)

#What is your question?
#Why is it interesting?
#What are you expecting?
#Based on previous data, what would we expect to see?

#R4All Philosophy
#Experimental design
#Bad data is bad data
#plot your data < important 
#theory drives expectation of a pattern
#make a model
#model and experimental match
#evaluate assumptions
#check the p-value
#evaluate model conclusions
#now have a look at F and p
#add interpretation to figure
#figures are better than tables

#data I need

ozone <- read_csv("Data/GardenOzone.csv")
daphnia <- read_csv("Data/Daphniagrowth.csv")
ladybirds <- read_csv("Data/ladybirds.csv")

#Test 1: The T-test ----

#data we are using
glimpse(ozone)

#make a figure: two histograms stacked on top of each other.
#Question: Why do we only have aes(x=) and no y=?
ggplot(ozone, aes(x=Ozone, fill=Garden.location))+
  geom_histogram(bins=10)+ #  colour="red", fill="blue"
  facet_wrap(~Garden.location, ncol=1)+ #separates into two graphs
  scale_fill_manual(values=c("deeppink", "navyblue"))+
  theme_bw()

#Null hypothesis: There is no significant difference in the Ozone levels between
#East and West garden locations

#Guestimate
#Mean ozone 80 east gardens
#Mean ozone 65 west gardens
#15 difference <- against the null hypothesis!

#check an assumption
#fill in the blanks with what you think
#we are aiming to calculate the mean and variance of each set of gardens
ozone %>%
  group_by(Garden.location) %>%
  summarise(meanOzone = mean(Ozone),
            varOzone = var(Ozone))

#do the test
#the first blank is the FUNCTION we use to do a t-test
#?t.test <- check notes, because variance isn't equal a Welsh modification of the t-test will be used
#that doesn't rely on the assumption that the variances are equal so it's okay
#to do the t-test
ozone_t_test <- t.test(Ozone ~ Garden.location, data=ozone)

#check out the results
ozone_t_test # p<0.001 to results are significant so we reject the null hypothesis



#what do the means help us to do?
#The means help us form an assumption about what we expect to see in our analysis
#The means were quite different do we made the assumption that the datasets were
#significantly different, with a null hypothesis saying they were NOT
#significantly different

#what does the 95% CI tell us?
#With 95% confidence we can say the true difference between the variables is between 8 and 24 

#what does the t-value tell you?
# t values gives us standard deviations from zero - I think!

#what does the p-value ACTUALLY mean?
# p value tells us how likely it is that we would get a difference of 16 by random under the 
# null hypothesis

#General Linear Model ---- 
#(not generalised linear model)

#include ANOVA, ANCOVA, regression, multiple regression, 
#non-linear via polynomials & transformations

#One response variable
#Covariates and Factors - continuous or discrete independent variable

#Balanced (=orthogonal design) and unbalanced designed
#e.g. two treatments: food and moisture
# 9 replicates in each combination = balanced (order of testing does not matter)
# different numbers of replicated for each combination = unbalanced (order of testing does matter)

#Linear models 

#Distribution assumptions (e.g. normailty & equal variance) are not to important unless your focus is prediction
#Getting the model structure is important:
 # is the additivity assumption met?
#  Have you captured the interactions?
 # The MOST IMPORTANT thing you need to do is to make sure the model addresses the research question!

#1-Way ANOVA (linear model)

#daphnia
#do parasites have an effect on growth?
#Which parasite is the worst?
#  Write down your uess at the NULL hypothesis for these data

#Null hypothesis: there is no significant difference in parasite treatments

glimpse(daphnia)

ggplot(daphnia, aes(y=growth.rate, x=parasite, fill=parasite, colour=parasite))+
  geom_boxplot()+ #  colour="red", fill="blue"
  theme_bw()

#Gestimate
#Control mean = 1.45
#Met. B mean = 0.7 (0.7-1.45)
#Pan. per mean = 1.2 (1.2-1.45)
#Past. ram mean = 0.55 <- most negative on growth (0.55-1.45)

daphnia %>%
  group_by(parasite) %>%
  summarise(meanGrowth = mean(growth.rate),
            varGrowth = var(growth.rate))

#graph seems to suggest I should reject my null hypotheses

#build the model

daph_mod <- lm(growth.rate ~ parasite, data=daphnia)

#view model
daph_mod

#what do the coefficients say?
#Intercept = growth.rate of control
#parasiteMet = -0.41 (reduced growth rate for this treatment)
#parasitePan = -0.14 (reduced growth rate for this treatment)
#parasitePast = -0.73 (reduced growth rate for this treatment)

#check assumptions
autoplot(daph_mod) #residuals, normal Q-Q, scale location, constant leverage

#Residuals (important)
#looking at heterostochasticity (variation around the means versus the fitted values)
#Want this to look like the sky at night. But won't get it because using treatments.
#If there are patterns the model has not been formulated correctly. 

#Normal Q-Q (important)
#Formally an assessment of the distribution of the assumption about the residuals.
#Dashed line represents what you would expect if residuals followed normal distribution exactly

#Scale-Location (important)
#Want to see that there is no relationship between mean and variance in residuals.
#basically line should be roughly straight. Variability should not increase or 
#decrease with mean values

#Constant Leverage: Residuals vs Factor Levels (not so important here)

#Now we examine the model!

#ANOVA(model.name) does not perform an ANOVA
#lm() produces an ANOVA
#ANOVA produces a table of ANOVA results

#Summary(model.name) produces a table of coefficients

#make strong inference
#what does this tell us? Which question does this answer?

anova(daph_mod)

#mean squared column parasite value 1.04, residual value 0.03
#large difference between what's explained by parasite variable (1.04) and what's
#left over
#1.04/0.032 = 32.5 = F value estimate of how much variation is explained by explanatory variable 
#relative to what's left over. We want a high number here! There are big differences between 
#our treatments

#contrasts... the summary table
#let's look at how different each treatment is
summary(daph_mod)

#treatments in summary put automatically in alphabetical order
#so control = intercept

#intercept value is mean growth rate of control
#estimates after are minus intercept to give mean

#These are treatment contrasts! 

#Reject null hypothesis! There is a difference between treatments.
#Two of the parasites were significantly different from the control.

#Test 3: Chi-square contingency table ----

glimpse(ladybirds)

#is there an association between the colourmorph of the ladybirds and the habitat
#they live in

#summarise the data - collect the observations into groups
LBtable <- ladybirds %>%
  group_by(Habitat, colour) %>%
  summarise(sumNum = sum(number))

#make the barplot - use the table!
ggplot(LBtable, aes(x=Habitat, y=sumNum, fill=colour))+
  geom_col(position="dodge")+
  scale_fill_manual(values=c(black="black", red="red"))+
  theme_bw()

#from the barplot it looks like there is a significant association between habitat and colour

#don't need to make assumptions here, we are ready to move on to the test
#must change out data table into a matrix
#chi-test requires a MATRIX
lb_chi_data <- xtabs(sumNum ~ Habitat + colour, data = LBtable)
lb_chi_data

#run the test
lb_mod <- chisq.test(lb_chi_data)
lb_mod
#Chi value 19 is big! Very unlikely to see a pattern like this if there was not an association
#p-vale < 0.001
names(lb_mod)

#build a picture before doing statistics with your data :)

