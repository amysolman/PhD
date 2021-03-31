##WRANGLING MY METADATA TABLE TO INCLUDE CATEGORIAL AND NUMERIC COLUMN NAMES FOR QIIME2

library(tidyverse)

metadata <- read.csv("../Data/Sommers2020/SraRunTable.txt")
x <- colnames(metadata)
df <- data.frame("q2:types", "categorical", "numeric", "numeric", "categorical", "categorical", "categorical",
                 "numeric", "categorical", "numeric", "numeric", "categorical", "categorical", "categorical", "categorical",
                 "numeric", "numeric", "numeric", "numeric", "numeric", "categorical", "categorical", "categorical", "categorical",
                 "categorical", "categorical", "categorical", "categorical", "numeric", "categorical", "categorical", "categorical", "categorical",
                 "numeric", "categorical", "numeric", "categorical", "numeric", "categorical", "categorical", "numeric", "numeric", "categorical",
                 "numeric", "categorical", "categorical", "numeric", "numeric")
colnames(df) <- x
new_df <- rbind(metadata, df)
new_new_df <- new_df %>% slice(160)
new_new_new_df <- new_df %>% slice(1:159)
new_new_new_new_df <- rbind(new_new_df, new_new_new_df)

#export wrangled dataframe
write.csv(new_new_new_new_df, "fastq-sommers2020/SraRunTable-wrangle.tsv")


