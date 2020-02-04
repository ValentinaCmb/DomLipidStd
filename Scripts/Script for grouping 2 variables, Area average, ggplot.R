source("QQQlipidomics_v0.5.R")
df <- import.lipids("StdCurves02.xlsx")
View(df)

library('tidyverse')

names(df)
names(df)[1]<-'Sample_Name'
head(Try)

rename(df, 'Sample_Name' = 'SampleName')

#calls the name of the columns and change the "Sample name" into "Sample_Name'

head(df)
#prints the head(x) 

Try<- df %>%
  select('Sample_Name', 'Lipid_Class', 'Area')
#objects name need the quotation mark 'xxx' otherwise R looks for them in the Environment list.

glimpse(Try)


Try %>% 
  group_by(Lipid_Class, Sample_Name) %>% 
  summarize (avg = mean(Area))
#make the avarege  based on Lipid class and Sample (std concentrations)

df2 <- Try %>% 
  group_by(Lipid_Class, Sample_Name) %>% 
  summarize (avg = mean(Area))


#rename(df2, 'avg' = 'Mean_Area')
#why is it not working? It doesn't change the name!


ggplot(df2, aes(x = Sample_Name, y = avg, colour=Lipid_Class)) +
 #scale_colour_brewer(Lipid_Class, palette="Dark2")+
  geom_point()
