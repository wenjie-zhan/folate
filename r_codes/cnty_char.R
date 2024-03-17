################################################################################
#### Set up
################################################################################
# Clean the Environment
rm(list=ls())

# Loading packages
# Packages
library(tidyverse)
library(data.table)
library(dplyr)
library(stargazer)
library(ggplot2)
library(zoo)
library(sandwich)
library(lmtest)
library(usmap)
library(magrittr) # %$%
library(readxl)
library(fredr) # extract data from St. Louis Fed
library(haven) # read_dta
library(fixest)

# Directories
setwd("C:/Files/davis_research/folate_data/")
repo <- "C:/Users/wenji/OneDrive/github_files/GitHub/folate/"

fitstat_register(type = "wt_mean", alias = "Dependent Variable Mean",
                 fun = function(x){
                   wt = weights(x) %>% as.data.frame() %>% filter(!is.na(.)) %>% rowid_to_column()
                   names(wt) = c("rowid", "wt")
                   y = model.matrix(x, type = "lhs") %>% as.data.frame() %>% rowid_to_column() %>%
                     merge(wt, by = "rowid", all.x = T)
                   names(y) = c("rowid", "y", "wt")
                   weighted.mean(y$y, y$wt/sum(y$wt), na.rm = T)} 
)

################################################################################
#### 1. Demographic characteristics (1988) from County Intercensal Estimates
################################################################################
# %black, %age_under_5, \%age_over65, \%female
intcen88<- read_excel("raw_data/cnty_char/county_intercensal_est/pe-02-1988.xls", skip = 5)
names(intcen88) <- c("year", "fips", "race_sex", "under5", "age5_9", "age10_14",
                     "age15_19", "age20_24", "age25_29", "age30_34", "age35_39",
                     "age40_44", "age45_49", "age50_54", "age55_59", "age60_64",
                     "age65_69", "age70_74", "age75_79", "age80_84", "over85")
intcen88 <- intcen88 %>% filter(!is.na(year))

blk88 <- intcen88 %>%
  filter(grepl("Black", race_sex) == TRUE) %>%
  rowwise() %>% 
  mutate(pop88 = sum(c_across("under5":"over85"), na.rm = T)) %>%
  group_by(fips) %>%
  summarise(blk88 = sum(pop88, na.rm = T))

fml88 <- intcen88 %>%
  filter(grepl("female", race_sex) == TRUE) %>%
  rowwise() %>% 
  mutate(pop88 = sum(c_across("under5":"over85"), na.rm = T)) %>%
  group_by(fips) %>%
  summarise(fml88 = sum(pop88, na.rm = T))

age88 <- intcen88 %>%
  rowwise() %>% 
  mutate(over65 = sum(c_across("age65_69":"over85"), na.rm = T)) %>%
  group_by(fips) %>%
  summarise(under5_88 = sum(under5, na.rm = T),
            over65_88 = sum(over65, na.rm = T))

dmgph88_pop <- intcen88 %>%
  rowwise() %>% 
  mutate(pop88 = sum(c_across("under5":"over85"), na.rm = T)) %>%
  group_by(fips) %>%
  summarise(pop88 = sum(pop88, na.rm = T)) %>%
  mutate(logpop88 = log(pop88)) %>%
  merge(blk88, by = "fips", all = T) %>%
  merge(fml88, by = "fips", all = T) %>%
  merge(age88, by = "fips", all = T) %>%
  mutate(f_blk88 = blk88*100/pop88,
         f_fml88 = fml88*100/pop88,
         f_under5_88 = under5_88*100/pop88,
         f_over65_88 = over65_88*100/pop88,
         fips = as.numeric(fips))

# population by county 1985
intcen85 <- read_excel("raw_data/cnty_char/county_intercensal_est/cnty_dmgph85.xls", 
                       skip = 1)
names(intcen85) <- c("year", "fips", "race_sex", "under5", "age5_9", "age10_14",
                     "age15_19", "age20_24", "age25_29", "age30_34", "age35_39",
                     "age40_44", "age45_49", "age50_54", "age55_59", "age60_64",
                     "age65_69", "age70_74", "age75_79", "age80_84", "over85")

pop85 <- intcen85 %>%
  rowwise() %>% 
  mutate(pop85 = sum(c_across("under5":"over85"), na.rm = T)) %>%
  group_by(fips) %>%
  summarise(pop85 = sum(pop85, na.rm = T)) %>%
  mutate(fips = as.numeric(fips))

# population by county 1986
intcen86 <- read_excel("raw_data/cnty_char/county_intercensal_est/cnty_dmgph86.xls", 
                       skip = 1)
names(intcen86) <- c("year", "fips", "race_sex", "under5", "age5_9", "age10_14",
                     "age15_19", "age20_24", "age25_29", "age30_34", "age35_39",
                     "age40_44", "age45_49", "age50_54", "age55_59", "age60_64",
                     "age65_69", "age70_74", "age75_79", "age80_84", "over85")

pop86 <- intcen86 %>%
  rowwise() %>% 
  mutate(pop86 = sum(c_across("under5":"over85"), na.rm = T)) %>%
  group_by(fips) %>%
  summarise(pop86 = sum(pop86, na.rm = T)) %>%
  mutate(fips = as.numeric(fips))

# population by county 1987
intcen87 <- read_excel("raw_data/cnty_char/county_intercensal_est/cnty_dmgph87.xls", 
                       skip = 1)
names(intcen87) <- c("year", "fips", "race_sex", "under5", "age5_9", "age10_14",
                     "age15_19", "age20_24", "age25_29", "age30_34", "age35_39",
                     "age40_44", "age45_49", "age50_54", "age55_59", "age60_64",
                     "age65_69", "age70_74", "age75_79", "age80_84", "over85")

pop87 <- intcen87 %>%
  rowwise() %>% 
  mutate(pop87 = sum(c_across("under5":"over85"), na.rm = T)) %>%
  group_by(fips) %>%
  summarise(pop87 = sum(pop87, na.rm = T)) %>%
  mutate(fips = as.numeric(fips))


################################################################################
#### 2. County Data Book 1988 and 1994
################################################################################
cdb88 <- read_dta("raw_data/cnty_char/county_data_book/DS0079/02896-0079-Data.dta") %>%
  select(fips = fips, 
         stfips = statefip,
         phy_pp85 = var44,
         bed_pp85 = var47,
         ssb_pp85 = var51,
         s_crime85 = var54,
         v_crime85 = var55,
         psch_enrl86n87 = var57,
         income_pp85 = var63,
         unemp_rate86 = var82,
         f_funds_pp86 = var177
  ) %>%
  mutate( fips = as.numeric(fips)
  ) %>%
  filter(fips != 0 & fips %% 1000 != 0) 

cdb94 <- read_dta("raw_data/cnty_char/county_data_book/DS0080/02896-0080-Data.dta") %>%
  select(fips = fips, 
         stfips = statefip,
         birth_rate88 = var048,
         death_rate88 = var051,
         v_farm_pa87 = var161,
         m_size_farm87 = var164,
         local_gov_exp86n87 = var215
  ) %>%
  mutate( fips = as.numeric(fips)
  ) %>%
  filter(fips != 0 & fips %% 1000 != 0) 

################################################################################
#### 3. Merge data sets, state-level baseline characteristics
################################################################################
cnty_char <- dmgph88_pop %>%
  merge(cdb88, by = "fips", all.x = T) %>%
  merge(cdb94, by = "fips", all.x = T) %>%
  merge(pop85, by = "fips", all.x = T) %>%
  merge(pop86, by = "fips", all.x = T) %>%
  merge(pop87, by = "fips", all.x = T) %>%
  mutate(stfips = fips %/% 1000,
         pub_sch_enrl_pp86n87 = psch_enrl86n87/pop87,
         local_gov_exp_pp86n87 = local_gov_exp86n87/pop87)

pop_st <-  cnty_char %>%
  select(stfips, pop85, pop86, pop87, pop88) %>%
  group_by(stfips) %>%
  summarise(pop85_st = sum(pop85, na.rm = T),
            pop86_st = sum(pop86, na.rm = T),
            pop87_st = sum(pop87, na.rm = T),
            pop88_st = sum(pop88, na.rm = T))

cnty_char_wt <- cnty_char %>%
  merge(pop_st, by = "stfips", all.x = T) %>%
  mutate(wt85 = pop85/pop85_st,
         wt86 = pop86/pop86_st,
         wt87 = pop87/pop87_st,
         wt88 = pop88/pop88_st
  )

cnty_char_st <- cnty_char_wt %>%
  group_by(stfips) %>%
  summarise(f_blk88 = weighted.mean(f_blk88, wt88, na.rm = T),
            f_fml88 = weighted.mean(f_fml88, wt88, na.rm = T),
            under5_88 = weighted.mean(under5_88, wt88, na.rm = T),
            over65_88 = weighted.mean(over65_88, wt88, na.rm = T),
            phy_pp85 = weighted.mean(phy_pp85, wt85, na.rm = T),
            bed_pp85 = weighted.mean(bed_pp85, wt85, na.rm = T),
            ssb_pp85 = weighted.mean(ssb_pp85, wt85, na.rm = T),
            s_crime85 = weighted.mean(s_crime85, wt85, na.rm = T),
            v_crime85 = weighted.mean(v_crime85, wt85, na.rm = T),
            psch_enrl86n87 = weighted.mean(psch_enrl86n87, wt87, na.rm = T),
            income_pp85 = weighted.mean(income_pp85, wt85, na.rm = T),
            unemp_rate86 = weighted.mean(f_funds_pp86_rs, wt86, na.rm = T),
            l_exp_pp81n82_rs = weighted.mean(l_exp_pp81n82_rs, wt82, na.rm = T),
            l_exp_health81n82_rs = weighted.mean(l_exp_health81n82_rs, wt82, na.rm = T),
            l_exp_welfare81n82_rs = weighted.mean(l_exp_welfare81n82_rs, wt82, na.rm = T),
            sales_food_store82_rs = weighted.mean(sales_food_store82_rs, wt82, na.rm = T),
            sales_eat_places82_rs = weighted.mean(sales_eat_places82_rs, wt82, na.rm = T),
            v_farmland_pa82 = weighted.mean(v_farmland_pa82, wt82, na.rm = T),
            v_crop_pa82 = weighted.mean(v_crop_pa82, wt82, na.rm = T)
  )

saveRDS(cdb88_wt, paste(dict_wd, "cnty_char.rds", sep = ""))


################################################################################
#### 4. Regression
################################################################################
cnsa_tmp <- read.csv(paste(dict_wd, "cnsa_all_st.csv", sep = "")) %>%
  mutate(stfips = as.numeric(fips))

dta_st <- cdb88_st %>%
  merge(cnsa_tmp[, c("stfips", "cnsar_nc89_93_res", "median")], by = "stfips",all = T) %>%
  merge(pop_st, by = "stfips", all = T ) %>%
  mutate(logpop86 = log(pop86_st))

dta_cnty <- cdb88_wt %>%
  merge(cnsa_tmp[, c("stfips", "cnsar_nc89_93_res", "median")], by = "stfips",all = T) %>%
  mutate(logpop86 = log(pop86))

st_reg_ct <- feols(
  cnsar_nc89_93_res ~ 
    f_blk80_rs
  + f_hsp80_rs
  + br_rate84_rs
  + mar_rate84_rs
  + edu12more80_rs
  + edu16more80_rs
  + ssb_pp85_rs
  + unemp86_rs
  + incpp84_rs
  # + phy_pp85_rs
  # + bed_pp85_rs
  + f_funds_pp86_rs
  + l_exp_pp81n82_rs
  + l_exp_health81n82_rs
  + l_exp_welfare81n82_rs
  + sales_food_store82_rs
  + sales_eat_places82_rs
  + v_farmland_pa82
  + v_crop_pa82
  ,
  data = dta_st,
  vcov = "hetero",
  weights = ~pop80_st
) 
st_reg_ct
summary(st_reg_ct)

st_reg_bi <- feols(
  median ~ 
    f_blk80_rs
  + f_hsp80_rs
  + br_rate84_rs
  + mar_rate84_rs
  + edu12more80_rs
  + edu16more80_rs
  + ssb_pp85_rs 
  + unemp86_rs
  + incpp84_rs
  # + phy_pp85_rs
  # + bed_pp85_rs
  + f_funds_pp86_rs
  + l_exp_pp81n82_rs
  + l_exp_health81n82_rs
  + l_exp_welfare81n82_rs
  + sales_food_store82_rs
  + sales_eat_places82_rs
  + v_farmland_pa82
  + v_crop_pa82,
  data = dta_st,
  vcov = "hetero",
  weights = ~pop80_st
) 
st_reg_bi
summary(st_reg_bi)

cnty_reg_ct <- feols(
  cnsar_nc89_93_res ~ 
    f_blk80_rs
  + f_hsp80_rs
  + br_rate84_rs
  + mar_rate84_rs
  + edu12more80_rs
  + edu16more80_rs
  + ssb_pp85_rs 
  + unemp86_rs
  + incpp84_rs
  # + phy_pp85_rs
  # + bed_pp85_rs
  + f_funds_pp86_rs
  + l_exp_pp81n82_rs
  + l_exp_health81n82_rs
  + l_exp_welfare81n82_rs
  + sales_food_store82_rs
  + sales_eat_places82_rs
  + v_farmland_pa82
  + v_crop_pa82,
  data = dta_cnty,
  cluster = ~stfips,
  weights = ~pop80
) 
cnty_reg_ct
summary(cnty_reg_ct)

cnty_reg_bi <- feols(
  median ~ 
    f_blk80_rs
  + f_hsp80_rs
  + br_rate84_rs
  + mar_rate84_rs
  + edu12more80_rs
  + edu16more80_rs
  + ssb_pp85_rs 
  + unemp86_rs
  + incpp84_rs
  # + phy_pp85_rs
  # + bed_pp85_rs
  + f_funds_pp86_rs
  + l_exp_pp81n82_rs
  + l_exp_health81n82_rs
  + l_exp_welfare81n82_rs
  + sales_food_store82_rs
  + sales_eat_places82_rs
  + v_farmland_pa82
  + v_crop_pa82,
  data = dta_cnty,
  cluster = ~stfips,
  weights = ~pop80,
) 
cnty_reg_bi
summary(cnty_reg_bi)

etable(st_reg_ct,  cnty_reg_ct, st_reg_bi, cnty_reg_bi,
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~  r2 + ar2 + n + wt_mean, tex = TRUE)







