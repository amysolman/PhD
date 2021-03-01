# NHM 2021 Getting started with R course
# Owen Petchey
# 18-Jan-2021

#install.packages("dplyr") ## don't do this!!!!

library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

# give me the whole numbers from 1 to 10
1:10
# assign the answer to x
x <- 1:10
# show me x
x

# Now lets do some work with some data
# and use the compensation.csv data to play with

# import the compensation.csv data file
fruit <- read_csv("data/compensation.csv")

## have a look at the data
glimpse(fruit)

# lets make a dataset with only the Fruit and Root variables
new_fruit <- select(fruit, Fruit, Root)
glimpse(new_fruit)

## do it a different way, by removing the Grazing variable
new_fruit <- select(fruit, -Grazing)

## give us particular row of a dataset
slice(fruit, 1:5)

## lets keep only the rows in which plants were grazed
grazed_fruit <- filter(fruit, Grazing == "Grazed")
glimpse(grazed_fruit)


## make a dataset that contains only the Root and Fruit measurements
## of the Ungrazed plants
new_data <- filter(fruit, Grazing == "Ungrazed")
new_data <- select(new_data, Root, Fruit)
new_data

## use the pipe --  %>%
new_data <- filter(fruit, Grazing == "Ungrazed") %>%
  select(., Root, Fruit)
new_data

new_data <- filter(fruit, Grazing == "Ungrazed") %>%
  select(Root, Fruit)
new_data

## filter on two conditions:
filter(fruit, Grazing == "Grazed" & Root > 7)

## create a new variable, based on one or more existing variables
fruit <- mutate(fruit, log10_Fruit = log10(Fruit),
                log10_Root = log10(Root))
fruit

