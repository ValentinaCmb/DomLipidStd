library(tidyverse)
library(readr)
getwd()

# LS = Lipid Search Software i.e. table of lipids etc obtained from LipidSearch
# CD =  Compund Discoverer

## ---- 1st function -----------

# to transorm the data obtained in lipid search in a format that is usable/we can work on

readLipidSearch <- function(filename) {
  skiprow <- read_tsv(filename)
  skiprow <- which(grepl("Result", skiprow$`Job Name`)) +3
  df <- read_tsv(filename, skip = skiprow) %>%
    separate(`FA Group Key`, c("TC", "DB"), sep = ":", remove = FALSE)
  }


LSdata <- readLipidSearch("./Data/LS_TableDomBatch01.txt")


ggplot(LSdata, aes(BaseRt, `Calc Mass`, colour = Class)) +geom_point() 

#some compounds look a bit dodgy so lets separate the different lipid classes
#but want to look to see if there are any difference
ggplot(filter(LSdata, Class == "TG"), aes(BaseRt, `Calc Mass`, colour = DB)) +geom_point() +
  facet_wrap(~Class)


# ---------- 2nd function ------------------------------------------------------------------

# following step 1, now we want to pull out  data from the lipidsearch table to make a new table that will run in CD

LS_to_CD_list <- function(df){
  #unassigned <- filter(df, Rej == 0) #need to change to do this only if there is a column called Rej
  calc <- df %>%
    ungroup() %>%
    select(Name = LipidMolec, LSFormula = Formula, RT = BaseRt, MainIon, `Calc Mass`, Class, FA, `FA Group Key`) %>%
    unique() %>%
    mutate(C = as.numeric(str_match(LSFormula, "\\bC(\\d+)")[,2]),
           H = as.numeric(str_match(LSFormula, "\\bH(\\d+)")[,2]),
           N = as.numeric(ifelse(grepl("N", LSFormula), 
                                 ifelse(grepl("N\\d+", LSFormula), 
                                        str_match(LSFormula, "\\bN(\\d+)")[,2], 1), 0)), 
           O = as.numeric(ifelse(grepl("O", LSFormula), 
                                 ifelse(grepl("O\\d+", LSFormula), 
                                        str_match(LSFormula, "\\bO(\\d+)")[,2], 1), 0)),
           P = as.numeric(ifelse(grepl("P", LSFormula), 
                                 ifelse(grepl("P\\d+", LSFormula), 
                                        str_match(LSFormula, "\\bP(\\d+)")[,2], 1), 0)),
           S = ifelse(grepl("S", LSFormula), 
                      ifelse(grepl("S\\d+", LSFormula), 
                             str_match(LSFormula, "\\bS(\\d+)")[,2], 1), ""),
           N = ifelse(MainIon == "+NH4", N + 1, N),
           H = ifelse(MainIon == "+NH4", H + 3, H),
           Name = ifelse(MainIon == "+NH4", paste(Name, "+NH4", sep = " "), Name),
           Formula = paste("C", C, "H", H, ifelse(N > 0, "N", ""), ifelse(N > 1, N, ""), 
                           "O", O, ifelse(P > 0, "P", ""), ifelse(P > 1, P, ""), sep = "")) 
  # %>%
    # mutate(Lipid_Class = ifelse(grepl("NO6", Formula) & !grepl("P", Formula), "TG", 
                               # ifelse(grepl("N2O5P", Formula), "CAEP", NA)),
           # Lipid_Class = ifelse(grepl("NO8P", Formula) & C %% 2 == 0, "PC" ,
               #                 ifelse(grepl("NO8P", Formula) & C %% 2 != 0, "PE", Lipid_Class)),
           # Lipid_Class = ifelse(C>30 & grepl("NO7P", Formula) & C %% 2 == 0, "PC",
                              #  ifelse(grepl("NO7P", Formula) & C %% 2 != 0, "PE", Lipid_Class)),
          # Lipid_Species = ifelse(grepl("NO6", Formula) & !grepl("P", Formula),
                                 # paste(C-3, ((C-3)*2-(H-5))/2, sep = ":"), NA)) 
  
  
  return(calc)
}

LSdata01 <- LS_to_CD_list(LSdata)
write_csv(LSdata01, path = "Data/Table_CD28Nov.csv")

# this below is not needed 
T01 <- LSdata01 %>% 
      select(Name, LSFormula, RT, MainIon, `Calc Mass`, Class, FA, `FA Group Key`) %>% 
      filter(Class %in% c("DG", "FA", "PE", "LPC", "LPE", "LPI", "LPS", "PA", "PC", "PG", "PI", "PS", "TG"))

#############################################################################################
#
#   START FROM HERE FOR STANDARDS ANALYSIS 
#
# ------------------------ 3rd part ------------------------------------------------- ######

# now that CD has run we want to plot this data and see what makes sense


require(tidyverse)
require(RColorBrewer)

CDdata <- read_tsv("Data/CompoundsFromCDUpdate.csv", col_types = cols(Name = col_character()))

key <- read_csv("Data/Standards_UpdatedTable_ComDisco.csv", col_types = cols(Name = col_character())) %>% 
      rename(LSFormula= Formula)

CDdata_gathered <- CDdata %>% 
  gather(sample, area, contains("Area: ")) %>%
  group_by(Name, Formula, `Molecular Weight`, `RT [min]`) %>%
  filter(!is.na(area),
         `RT [min]` < 25) %>% 
  mutate(count = n()) %>% 
  left_join(key) %>% 
  mutate(TC = word(`FA Group Key`, 1, 1, sep = "_"),
         DB = word(`FA Group Key`, -1, -1, sep = "_")) %>% 
  mutate(Batch = str_extract(sample, "[B]\\d+"),
         StandardsName =str_extract(sample, "\\d+"))



PC <- CDdata_gathered %>% 
  filter(Class == "PC")

ggplot( aes(`RT [min]`, `Molecular Weight`, colour = DB)) +
  geom_point()


# ______ THESE BELOW ARE NEEDED ______________ #

p1 <- CDdata_gathered %>% 
    group_by(Batch)

ggplot(p1, aes(`StandardsName`, log(area), colour = Batch)) +
  #geom_boxplot()+
  geom_point() +
  geom_point() +
  facet_wrap(~Class, scales ="free_y")




# --------- Part 4: make calibration table etc out of standards curve ----------------------------------------------------------- #

library(modelr)


# Assign concetrations (0.001; 0.01; 1; 10; 50 mg/L) to the Standards names (3, 4, 5 , 6, 13, 14 etc)
# take p1
p3 <-  p1 %>% 
      mutate(Conc =  as.numeric(ifelse (StandardsName == "3", "0.1",
                     ifelse (StandardsName == "13", "0.1",  
                     ifelse (StandardsName == "4", "1", 
                     ifelse (StandardsName == "14", "1", 
                     ifelse (StandardsName == "5", "10",
                     ifelse (StandardsName == "15", "10",
                     ifelse (StandardsName == "6", "50",
                     ifelse (StandardsName == "16", "50", "NA")))))))))) %>% 
  select(-Checked, -('Annotation Source: Predicted Compositions': 'FISh Coverage'), -(LSFormula : Mass))

    

p3

# Log transform the data 

ggplot(p3, aes(log(Conc), log(area), colour = Batch)) +
  geom_point() +
  geom_point() +
  facet_wrap(~Class, scales ="free_y")


# Add the Linear Model (lm) and remove colour =  Batch to see 1 regression line only 

ggplot(p3, aes(log(Conc), log(area))) +
  geom_point() +
  facet_wrap(~Class, scales ="free_y") +
  geom_smooth( method = "lm")


#  get the linear model thingy for slope and intercept for each Lipid class  
# (when code becomes long copy and paste it again and it will format itself automatically)
# "mod" is part of the liner model output but it isn't important for theis analysis so it gets removed with "select all but mod
stdlm <- p3 %>% 
  group_by(Class) %>% 
  do(mod= lm(log(area) ~ log(Conc), data = .)) %>% 
  mutate(slope = summary(mod)$coeff[2],
         intercept = coef(mod)[1],
                          adj.r.squared = summary(mod)$adj.r.squared) %>% 
  select(-mod)



# recall the summary of the LM for the equation
summary(stdlm)

# create a table with kable function as a figure. cannot be used for calculations! ____________

install.packages("kableExtra")
library(kableExtra)

standards <- stdlm %>% 
  kable() %>% 
  kable_styling()


standards

# ____ Calculate the Concentration of each compound __________

# AREA USED IS "area" NOT "Area max" !

sample.values <-  p3 %>% 
  left_join(stdlm) %>% 
  mutate(conc_mgmL = exp((log(area) - intercept)/ slope)) %>% 
  select( -slope, -intercept, -adj.r.squared)
  

sample.values  

# ________ plot the samples based on the conc and not area anylonger ________

ggplot(sample.values, aes(log(conc_mgmL), log(area))) +
  geom_point() +
  facet_wrap(~Class, scales ="free_y") +
  geom_smooth( method = "lm")

# make a table of all the standards _____

allStandards <- sample.values %>% 
  kable() %>% 
  kable_styling()

allStandards
