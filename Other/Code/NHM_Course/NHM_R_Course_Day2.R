#NHM R Course Day Two 19/01/21

#RStudio Projects are a different way of setting
#working directory.
#It basically links you to the folder where your
#scripts and data are and you can move it
#wherever you like. Means you can send others the 
#folder with the project in it and it will work easily!

#tidyverse contains dyplyr and ggplot
#and lots of other amazing packages

#slice
slice(comp, 2) #slice second row
slice(comp, 2:5) #slice row 2-5
slice(comp, c(2,5)) #slice row 2 and 5

#upiping
comp %>% #name of the dataframe
  select(Fruit) %>% #select column Fruit
  slice(8) #slice row 8

#subsetting
#show me all the root values less than or 
#equal to 10
comp %>%
  filter(Root<=10)

#create a new data frame containing only 
#rows where the Fruit values >80 

Big_Fruit <- comp %>%
  filter(Fruit>80)

#Create a new data frame containing only
#the ungrazed rows
Ungrazed <- comp %>%
  filter(Grazing == "Ungrazed")

#using two dplyr functions, and piping, 
#create a data frame containing the **Root**
#values where Fruit is greater than 80

Two_functions <- comp %>%
  filter(Fruit>80) %>%
  select(Root)

#creating new variables

comp %>%
  mutate(Root_cm = Root/100)

#rearrange data
comp%>%
  arrange(Fruit)

#rename function
comp%>%
  rename(FruityMcFruitFace = Fruit)

#summarising the data
#What functions in R can we use to get the
#following?
#mean, median, SD, variance, sample size, SE

sqrt(var(x))/length(x)

#give me the mean of the Fruit column
summary.data <- comp%>%
  summarise(meanF=mean(Fruit),
sdF = sd(Fruit),
lengthF = length(Fruit),
sumF=sum(Fruit))

#summarise by groups

comp%>%
  group_by(Grazing)%>%
  summarise(meanF=mean(Fruit))

#we can also string together lots of dplyr functions

comp%>%
  mutate(FRratio = Fruit/Root) %>% #give us a fruit to root ratio
  group_by(Grazing)%>% #look at only the grazing group
             summarise(meanFRrat = mean(FRratio)) #what is the mean of those ratios

#calculate the median, variance and sample size of the Fruit/Root Ratio 
#in each Grazing category in the compensation data, and assing this data 
#frame to a word called MrPedantic

MrPedantic <- comp%>%
  mutate(FRratio = Fruit/Root) %>% #give us a fruit to root ratio
  group_by(Grazing)%>%
  summarise(medianFRrat = median(FRratio),
            varFRrat = var(FRratio),
            lenFRrat = length(FRratio)) #can put n instead of length()

MrPedantic

write_csv(MrPedantic, "MrPedantic.csv")

#using ggplot2

library(ggthemes)
install.packages("ggthemes")
#need to tell ggplot about the dataframe,
#the aesthetic components (x and y axis)
#and the geometric objects (points, lines, bars)

ggplot(comp, aes(x=Root, y=Fruit))+ #change colour based on grouping
  geom_point(colour="#6ba292", size=4,
             shape="triangle") #change all colour


ggplot(comp, aes(x=Root, y=Fruit, colour= Grazing))+ #change colour based on grouping
  geom_point(size=4,
             shape="triangle") +
  scale_colour_manual(values=c("springgreen", "deeppink"))+#scale = change, colour = colour, manual = manually
#theme_classic()+
 # theme_bw()
 # theme_excel() #these are from the ggthemes package
  theme_economist()+
  xlab("Hello booger")+
  ylab("Bugger off!")

colors() #will show you all the built-in colours in R
#can also use hex colours #


ggplot(comp, aes(x=Root, y=Fruit, colour=Grazing))+
  geom_point(size=7,
             shape="triangle") +
  scale_colour_manual(values=c("springgreen", "deeppink"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="red", size=5))

ggsave("ugly.png", width=10, height=10)  

ggplot(comp, aes(x=Fruit, fill=Grazing))+
  geom_histogram(bins=11)+ #  colour="red", fill="blue"
  facet_wrap(~Grazing, ncol=1)+ #separates into two graphs
  scale_fill_manual(values=c("deeppink", "navyblue"))+
  theme_bw()
geom_histogram()

#plot boxplot!

ggplot(comp, aes(y=Fruit, x=Grazing, colour=Grazing, fill=Grazing))+
  geom_boxplot()+
  scale_colour_manual(values=c("springgreen", "deeppink"))+
  theme_bw()

install.packages("patchwork") #use to put multiple different plots on same plot
library(patchwork)

plot1 <- ggplot(comp, aes(y=Fruit, x=Grazing, colour=Grazing, fill=Grazing))+
  geom_boxplot()+
  scale_colour_manual(values=c("springgreen", "deeppink"))+
  theme_bw()

plot2 <- ggplot(comp, aes(x=Fruit, fill=Grazing))+
  geom_histogram(bins=11)+ #  colour="red", fill="blue"
  facet_wrap(~Grazing, ncol=1)+ #separates into two graphs
  scale_fill_manual(values=c("deeppink", "navyblue"))+
  theme_bw()
geom_histogram()

#two plots
plot1 + plot2

#three plots
(plot1+plot1)/plot2


  
