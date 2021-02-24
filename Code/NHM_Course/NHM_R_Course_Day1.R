# NHM 2021 Getting Started with R course
# With Owen Petchey
# 18 - Jan - 2021

# install dplyr and other packages
install.packages("dplyr")
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

#give me the whole numbers from 1 to 10
1:10
#assign the answer to x
x <- 1:10
#show me x
x

#Import untidy data into R and tidy there 

#Now let's do some work with some data
#and use the compensation.csv data to play with
#using underscore instead of . gives column specification
#of data
#The read_csv function imports data into R as a tibble, 
#while read. csv imports a regular old R data frame instead. 
#Tibbles are better than regular data frames because they: 
#load faster
#fruit changed to comp
comp <- read_csv("data/compensation.csv")

#have a look at the data
glimpse(comp)

#lets make a dataset with only the Fruit and Root variables
new_fruit <- select(fruit, Fruit, Root)
glimpse(new_fruit)

#do it a different way, by removing the Grazing variable
new_fruit <- select(fruit, -Grazing)

#give us particular rows of a dataset
slice(fruit, 1:5)

#lets keep only the rows in which plants were grazed
grazed_fruit <- filter(fruit, Grazing == "Grazed")
glimpse(grazed_fruit)

#make a dataset that contains only the root and fruit measurements
#of the ungrazed plants

ungrazed_fruit <- filter(fruit, Grazing == "Ungrazed")
ungrazed_fruit <- select(ungrazed_fruit, -Grazing)

#use the pipe %>%
# ., is a placeholder for what is being sent
new_data <- filter(fruit, Grazing == "Ungrazed") %>%
  select(., Root, Fruit)

#or you can leave the ., out!
new_data <- filter(fruit, Grazing == "Ungrazed") %>%
  select(Root, Fruit)

#filter on two conditions
filter(fruit, Grazing == "Grazed" & Root >7)

#create a new variable based on one or more existing variables
fruit <- mutate(fruit, log10_Fruit = log10(Fruit),
                l10_Root = log10(Root))
glimpse(fruit)
