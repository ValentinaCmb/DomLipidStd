source("QQQlipidomics_v0.5.R")
require(RColorBrewer)

df <- import.lipids("Reproducibility Test01.xlsx")
MRMlist <- read_csv("Bactrocera_dMRM.csv")

df1 <- df%>%
  mutate(`Compound Name` = paste(Lipid_Class, Lipid_Species, Acyl_Chain, sep = " "),
         SampleGroup = ifelse(grepl("34|35|36|37", `Sample Name`), "5 mg", 
                              ifelse(grepl("29|30|31|32", `Sample Name`), "2mg", 
                                     ifelse(`Sample Name` == "Blank", "Blank", "Test")))) %>%
  left_join(MRMlist)  %>%
  filter(Area > 1000,
         #!is.na(Polarity),
         !grepl("p|e", Lipid_Species),
         !grepl("PG", Lipid_Class),
         !grepl("PS 34:2|PI 32:1", `Compound Name`),
         grepl("mg", SampleGroup)) 


RTplot <- ggplot(df1, aes(RT, `Precursor Ion`, colour = desaturation)) +
  geom_point() +
  facet_wrap(~Lipid_Class) +
  #theme_bw() +
  scale_colour_brewer(palette = "Paired")

RTplot

## Different graphs showing variation for each different measurement. 
## Stats on different measurements.
df2 <- df1 %>%
  group_by(`Sample Name`, SampleGroup, Lipid_Class, Lipid_Species, TotalCarbons, desaturation, `Precursor Ion`) %>%
  summarise(Area = sum(Area), n = n())

ggplot(df2, aes(Lipid_Species, Area, colour = SampleGroup)) + 
  geom_boxplot() +
  facet_wrap(~Lipid_Class, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


#make into two groups of 5 and 2 mg
#use as colour

# use dplyr to group data and sum lipid species measurements
# group_by(-Acyl_Chain) %>% summarise(n = n(), x = sum(Area))

TAGplot <- ggplot(filter(df2, Lipid_Class == "TAG"), aes(Lipid_Species, Area, colour = SampleGroup)) + 
  geom_boxplot() +
  facet_wrap(~Lipid_Class, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

TAGplot

