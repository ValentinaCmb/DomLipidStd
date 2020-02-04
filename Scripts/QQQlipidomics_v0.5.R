## V0.3 updated for tidyverse packages and more complicated and complete Agilent Mass Hunter output in xlsx format
##V0.4 remove metadata information for single file analysis. 

require(tidyverse)
require(readxl) 
require(stringr) #probably not required


#samples were diluted into a 96 well plate 


import.lipids <- function(StdCurves02) { 
  df1 <- read_excel(StdCurves02, col_names = FALSE, skip = 2)
  heading <- read_excel(StdCurves02, n_max = 2, col_names = FALSE)
  
  tidycolumnnames <-  tibble(x= as.vector(t(heading[1,])), y = as.vector(t(heading[2,]))) %>%
    fill(x) %>%
    filter(!is.na(y)) %>%
    mutate(header=paste(x, y, sep = " "))
  
  df1 <- setNames(df1, tidycolumnnames$header)
  
  df2 <- df1 %>%
    gather(compound, response, contains("Results")) %>%
    separate(compound, c("Lipids", "Measurement"), sep = " Results ") %>%
    mutate(Lipids = str_replace_all(Lipids, "_", ":")) %>%
    spread(Measurement, response) %>%
    filter(!is.na(Area))
  
  PC <- df2%>%
    filter(grepl("PC|FFA|LPA", Lipids),
           !grepl("Chlor", Lipids)) %>%
    separate(Lipids, c("Lipid_Class", "Lipid_Species"), sep =  " ") %>%
    mutate(Acyl_Chain = "na")
  
  lipids <- df2 %>%
    filter(!grepl("PC|FFA|LP|Chlor", Lipids)) %>% #removes lyso phospholipids, PA and FFA and Chlorophyll
    separate(Lipids, c("Lipid_Class", "Lipid_Species", "Acyl_Chain"), sep =  " ") %>%
    mutate(Acyl_Chain = gsub("C", "", Acyl_Chain),
           Lipid_Class = ifelse(Lipid_Class == "TG", "TAG", ifelse(Lipid_Class == "DG", "DAG", Lipid_Class)))
  
  df3 <- rbind(lipids, PC) %>%
    separate(Lipid_Species, c("TotalCarbons", "desaturation"), sep = ":", remove = FALSE)
  #add chlorophyll and other compounds in as well
  Chlorophyll <- df2 %>%
    filter(grepl("Chlor", Lipids)) %>%
    mutate(Lipid_Class = Lipids,
           Lipid_Species = "na",
           Acyl_Chain = "na",
           TotalCarbons = "na",
           desaturation = "na") %>%
    select(-Lipids)
  df4 <- rbind(df3, Chlorophyll)
  return(df4)
  
}

#assign MW

MW.table <- read.table(header = TRUE, text = '
                       Lipid_Class  Lipid_Species  Acyl_Chain  MW  purity
                       TAG  54:3  18:1  884.7833  100
                       DAG  36:2  18:1  620.538  100
                       PC  36:2  na  785.593  100
                       PC  36:2  18:1  785.593  100
                       LPC  18:1  na  521.348  100
                       MGDG  34:6  18:3  746.4969  66.88
                       DGDG  34:3  18:3  914.5967  33.69
                       PE  34:1  16:0  716.6  100 ###check molecular weights, do we need Na adduct for PG and PI?
                       PG  36:1  18:1  775.7  100
                       PI  34:2  16:0  832.5102  50.27
                       MGDG 34:6  16:3  746.4969  66.88
                       PE  34:1  18:1  716.6  100
                       PG  36:1  18:0  775.7  100
                       PI  34:2  18:2  832.5102  50.27
                       DGDG  34:3  16:0  914.5967  33.69
                       PA  36:2  18:1  700.5043  100'
)
MW.table

standards <- read.table(header = TRUE, stringsAsFactors = FALSE, text = '
                        Lipid_Class  Std_Class  Mult_Factor  Lipid_Type
                        PI  PI  1  Phospholipid
                        PG  PG  1  Phospholipid
                        PE  PE  1  Phospholipid
                        DAG  DAG  1  Glycerolipid
                        TAG  TAG  1  Glycerolipid
                        MGDG  MGDG  1  Galactolipid  # galactolipid is not pure. 
                        DGDG  DGDG  1  Galactolipid  #
                        PC  PC  1  Phospholipid
                        LPC  PC  1  Phospholipid
                        PA  PA  1  Phospholipid') 

###function for standards set up correctly
#sample names std 0.01 - std 1000.0  Updated method might be suitable to 0.001 ug/ml

stdfunction <- function(dfinput, AdditionalStdClassGrepl = "") {
  stintard <- as.data.frame(left_join(dfinput, MW.table)) %>%
    filter(!is.na(MW),
           !grepl("Soy", `Sample Name`),
           #!grepl("0.001", `Sample Name`),
           grepl("std|Std|STD|Galacto", `Sample Name`)#,
           #grepl("1", Name)
    ) %>%
    group_by(`Sample Name`, `Sample Data File`, Lipid_Class, Lipid_Species, MW, purity) %>%
    summarise(
      response = sum(Area), 
      n = n()
    ) %>%
    mutate(conc = as.numeric(str_extract(`Sample Name`, "\\d+[.]\\d+"))*purity/100,
           conc = ifelse(is.na(conc), as.numeric(str_extract(`Sample Name`, "\\d+"))*purity/100, conc))
  prep <- function(df) {
    lipidclassSelect <- df %>%
      mutate(concmM = conc/MW) %>%    #µmole/ml
      mutate(conc_nmolml = concmM  * 1000) ## convert result from umol/ml to nmol/ml (need to take into consideration the dilution of dry weight etc later)
  }
  
  glycero <- stintard %>%
    filter(!grepl(AdditionalStdClassGrepl, `Sample Name`),
           !grepl(AdditionalStdClassGrepl , Lipid_Class))%>%
    prep()
  
  galacto <- stintard %>% #probably redundant if the galacto standards are run with other standards
    filter(grepl(AdditionalStdClassGrepl, `Sample Name`),
           grepl(AdditionalStdClassGrepl , Lipid_Class))%>%
    prep()
  
  stdcurve <- rbind(galacto, glycero) 
  
  stdplot <- ggplot(stdcurve, aes(log(conc), log(response))) +
    geom_point() +
    facet_wrap("Lipid_Class", ncol = 3, scales = "free_y") +
    theme_bw() + 
    geom_smooth(aes(group=1,color=Lipid_Class), method="lm", se=FALSE)
  stdlm <- stdcurve %>%
    group_by("Std_Class" = Lipid_Class) %>%
    do(mod = lm(log(response) ~ log(conc_nmolml), data = .)) %>%
    mutate(slope = summary(mod)$coeff[2],
           intercept = coef(mod)[1],
           adj.r.squared = summary(mod)$adj.r.squared) %>%
    select(-mod)
  return(list(stdcurve, stdlm, stdplot))
}


sample.values <- function(dataframe, stds) {
  samples <- dataframe%>%
    filter(!grepl("std|lank|DG|pho", `Sample Name`)) %>%
    left_join(standards) %>%
    left_join(data.frame(stds[2])) %>%
    mutate(conc_nmolmg = exp((log(Area) - intercept)/slope)) %>%
    select(-slope, -intercept, -adj.r.squared)
}

totalclass <- function(samples) {
  lipidclass <- samples %>%
    filter(!is.na(conc_nmolmg)) %>%
    ungroup() %>%
    group_by(`Sample Name`, Lipid_Class, Lipid_Type) %>%
    summarise(
      speciestotal = sum(conc_nmolmg)
    ) 
  return(lipidclass)
}

totalspecies <- function(df) {
  species  <- df%>%
    filter(!grepl("31|37" , TotalCarbons),
           !is.na(conc_nmolmg)) %>% #remove phospho internal standards
    group_by(`Sample Name`, Lipid_Class, Lipid_Species, Lipid_Type, TotalCarbons, desaturation) %>%
    summarise(
      conc.nmol.mg = sum(conc_nmolmg)
    )
  return(species)
}


totalacyl <- function(samples) {
  acyl  <- samples %>%
    group_by(`Sample Name`, Lipid_Class, Acyl_Chain, Lipid_Type) %>%
    summarise(
      conc.nmol.mg = sum(conc_nmolmg)
    ) 
  return(acyl)
}

totallipids <- function(samples) { #this one is now redundant, I think. 
  lipids  <- samples %>%
    group_by(`Sample Name`, Lipid_Class, Lipid_Species, Acyl_Chain, Lipid_Type) %>%
    summarise(
      conc.nmol.mg = sum(conc_nmolmg)
    ) 
  return(lipids)
}


sample.values2 <- function(dataframe, stds) {
  samples <- dataframe%>%
    filter(!grepl("std|lank|DG|pho", `Sample Name`)) %>%
    left_join(standards) %>%
    left_join(data.frame(stds[2])) %>%
    mutate(conc_nmolmg = exp((log(Area) - intercept)/slope)) %>%
    select(-slope, -intercept, -adj.r.squared)
  
  lipidclass <- samples %>%
    filter(!is.na(conc_nmolmg)) %>%
    group_by(`Sample Name`, Lipid_Class, Lipid_Type) %>%
    summarise(
      conc.nmol.mg = sum(conc_nmolmg))
      
  species  <- samples %>%
    filter(!grepl("31|37" , TotalCarbons),
           !is.na(conc_nmolmg)) %>% #remove phospho internal standards
    group_by(`Sample Name`, Lipid_Class, Lipid_Species, Lipid_Type, TotalCarbons, desaturation) %>% 
    summarise(
      conc.nmol.mg = sum(conc_nmolmg)
    )
  
  acyl  <- samples %>%
    group_by(`Sample Name`, Lipid_Class, Acyl_Chain, Lipid_Type) %>%
    summarise(
      conc.nmol.mg = sum(conc_nmolmg)
    ) 
  #this one is now redundant, I think. 
  lipids  <- samples %>%
    group_by(`Sample Name`, Lipid_Class, Lipid_Species, Acyl_Chain, Lipid_Type) %>%
    summarise(
      conc.nmol.mg = sum(conc_nmolmg)
    ) 
  return(list("Original" = samples, "SumLipidClass" = lipidclass, "SumLipidSpecies" = species, "Acyl" = acyl, "SpeciesAcyl" = lipids))
}

classplot <- function(df, grepclass) {
 df <- df %>%
  filter(grepl(grepclass, Lipid_Class))
ggplot(df, aes(Lipid_Class, mean.nmol.mg, fill = `Sample Name`)) +
 geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=mean.nmol.mg, ymax=mean.nmol.mg+se_conc),position=position_dodge())+
facet_wrap(facets = "tissue", scales = "free") +
theme_bw()+
    scale_fill_manual(name = "genotype", values = myColours) +
    theme(axis.text.x = element_text(angle = 90))
}
classplot
