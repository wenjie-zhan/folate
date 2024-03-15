################################################################################
#### Set up
################################################################################
# Clean the Environment
rm(list=ls())

# Load packages
library(tidyverse)
library(tidycensus)
library(dplyr)
library(readxl)
library(stargazer)
library(rlist)
library(ggplot2)
library(haven)
library(zoo)
library(sandwich)
library(lmtest)
library(fixest)
library(fredr)
library(DescTools)

conflicted::conflicts_prefer(dplyr::lag)
conflicted::conflicts_prefer(dplyr::filter)

# Dictionaries
setwd("C:/Files/davis_research/folate_data/")
repo <- "C:/Users/wenji/OneDrive/github_files/GitHub/folate/"

# Manually update the "stargazer" package to fix some problems in the older version
##WZ: R 4.2.0 or later will fail the stargazer if the the model names of the input are too long.
#Solution is from: https://gist.github.com/alexeyknorre/b0780836f4cec04d41a863a683f91b53
#Quick fix for stargazer <= 5.2.3 is.na() issue with long model names in R >= 4.2
# Unload stargazer if loaded
# detach("package:stargazer",unload=T)
# Delete it
# remove.packages("stargazer")
# Download the source
# download.file("https://cran.r-project.org/src/contrib/stargazer_5.2.3.tar.gz", destfile = "stargazer_5.2.3.tar.gz")
# Unpack
# untar("stargazer_5.2.3.tar.gz")
# Read the sourcefile with .inside.bracket fun
# stargazer_src <- readLines("stargazer/R/stargazer-internal.R")
# Move the length check 5 lines up so it precedes is.na(.)
# stargazer_src[1990] <- stargazer_src[1995]
# stargazer_src[1995] <- ""
# Save back
# writeLines(stargazer_src, con="stargazer/R/stargazer-internal.R")
# Compile and install the patched package
# install.packages("stargazer", repos = NULL, type="source")
# library(stargazer)
# The output can be put into the Latex, but you may need to load the package "dcolumn" to your Overleaf if you have not done that yet.
################################################################################
#### Get ACS 2016-2019, 2021-2022 from ipums.org
################################################################################
#Loading and organizing the PUMS, obtaining the state of residency/household
acs16_22 <- read_dta("raw_data/acs/acs16_22.dta")
acs_tmp <- acs16_22 %>%
  filter(birthqtr %in% c(1:4) & bpl <= 56 ) 

################################################################################
#### %Student by age
################################################################################
age <- c()
stu_all <-c()
stu_hs <- c()
stu_clg <- c()
stu_gs <- c()

for (i in 1:19) {
  age <- age %>% 
    append(10 + i) 
  stu_all <- stu_all %>%
    append(nrow(acs_tmp[acs_tmp$school == 2  & acs_tmp$age == 10 + i,])/nrow(acs_tmp[acs_tmp$age == 10 + i,]))
  stu_gs <- stu_gs %>%
    append(nrow(acs_tmp[acs_tmp$gradeatt == 7  & acs_tmp$age == 10 + i,])/nrow(acs_tmp[acs_tmp$age == 10 + i,]))
  stu_clg <- stu_clg %>%
    append(nrow(acs_tmp[acs_tmp$gradeatt == 6  & acs_tmp$age == 10 + i,])/nrow(acs_tmp[acs_tmp$age == 10 + i,]))
  stu_hs <- stu_hs %>% 
    append(nrow(acs_tmp[acs_tmp$gradeatt < 6 & acs_tmp$school == 2  & acs_tmp$age == 10 + i,])/nrow(acs_tmp[acs_tmp$age == 10 + i,]))
}

stu_ts_tmp1 <- data.frame(age, stu_all) %>%
  mutate(group = "All Types of Students") %>%
  rename(stu_rate =stu_all )
stu_ts_tmp2 <- data.frame(age, stu_hs) %>%
  mutate(group = "High School or Lower Grades") %>%
  rename(stu_rate = stu_hs)
stu_ts_tmp3 <- data.frame(age, stu_clg) %>%
  mutate(group = "College Student or Above") %>%
  rename(stu_rate = stu_clg)
stu_ts_tmp4 <- data.frame(age, stu_gs) %>%
  mutate(group = "Graduate Student or Above") %>%
  rename(stu_rate = stu_gs)

stu_ts <- data.table::rbindlist(list(stu_ts_tmp1, stu_ts_tmp2, stu_ts_tmp3, stu_ts_tmp4))
stu_ts$group <- stu_ts$group %>%
  factor(levels = c("All Types of Students", "High School or Lower Grades", 
                    "College Student or Above", "Graduate Student or Above"), ordered = T)

ggplot(data = stu_ts, aes(y = stu_rate, x = age, group = group), ) +
  geom_point(aes(color = group, shape = group), size = 2.5) +
  geom_line(aes(color = group)) +
  scale_color_manual(values = c("#144062", "#1e6091", "#168aad", "#76c893")) + 
  geom_vline(xintercept = 19, linetype = 2)+
  scale_x_continuous(labels = 11:29, breaks = 11:29)+
  xlab("Age") +
  ylab("%Students") +
  theme_classic()+
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), axis.title = element_text(size = 12)) +
  guides(color = guide_legend(nrow = 2))

ggsave("stu_rate_res.png", path = paste(repo, "r_outputs/", sep = ""), width = 6.472, height = 4.5)

################################################################################
#### Load CNSA data
################################################################################
# D1-2 Create variables
# Loading and organizing other data
cnsa_res <- read.csv("raw_data/vital_stat/cnsa_res.csv")
cnsa_fh <- read.csv("raw_data/vital_stat/cnsa_fh.csv")
ntd <- read.csv("raw_data/vital_stat/ntd.csv")
cnsa_brfd <- read.csv("raw_data/vital_stat/cnsa_brfd.csv")
# Merge cnsa_res and cnsa1
cnsa_tmp <- cnsa_res %>%
  select(fips, cnsar_nc89_93_res, cnsar93_res, cnsar_cp89_92_res, br89_93_res, br89_92_res, br93_res) %>%
  merge(cnsa_fh[, c("fips", "cnsar_fh", "br_fh")], by = "fips", all = T) %>%
  merge(ntd[, c("fips", "ntdr")], by = "fips", all = T) %>%
  merge(cnsa_brfd[, c("fips", "cnsar_brfd")], by = "fips", all = T)

cnsa_tmp$fips <- as.numeric(cnsa_tmp$fips) 

# Continuous cnsa to high-low binary cnsa
cnsa_tmp$qt <- case_when(
  cnsa_tmp$cnsar_nc89_93_res > Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.8, na.rm = TRUE) ~ 5,
  cnsa_tmp$cnsar_nc89_93_res > Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.6, na.rm = TRUE) & cnsa_tmp$cnsar_nc89_93_res <= Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.8, na.rm = TRUE) ~ 4,
  cnsa_tmp$cnsar_nc89_93_res > Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.4, na.rm = TRUE) & cnsa_tmp$cnsar_nc89_93_res <= Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.6, na.rm = TRUE) ~ 3,
  cnsa_tmp$cnsar_nc89_93_res > Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.2, na.rm = TRUE) & cnsa_tmp$cnsar_nc89_93_res <= Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.4, na.rm = TRUE) ~ 2,
  cnsa_tmp$cnsar_nc89_93_res <= Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.2, na.rm = TRUE) ~ 1
)

cnsa_tmp$mean <- ifelse(cnsa_tmp$cnsar_nc89_93_res >= weighted.mean(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, na.rm = TRUE), 1, 0) 
cnsa_tmp$median <- ifelse(cnsa_tmp$cnsar_nc89_93_res >= Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.5, na.rm = TRUE), 1, 0)
cnsa_tmp$top40 <- ifelse(cnsa_tmp$cnsar_nc89_93_res >= Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.6, na.rm = TRUE), 1, 0)
cnsa_tmp$top30 <- ifelse(cnsa_tmp$cnsar_nc89_93_res >= Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.7, na.rm = TRUE), 1, 0)
cnsa_tmp$top20 <- ifelse(cnsa_tmp$cnsar_nc89_93_res >= Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.8, na.rm = TRUE), 1, 0)

cnsa_tmp$mean_fh <- ifelse(cnsa_tmp$cnsar_fh >= weighted.mean(cnsa_tmp$cnsar_fh, weights = cnsa_tmp$br_fh, na.rm = TRUE), 1, 0)
cnsa_tmp$mean93_res <- ifelse(cnsa_tmp$cnsar93_res >= weighted.mean(cnsa_tmp$cnsar93_res, weights = cnsa_tmp$br93_res, na.rm = TRUE), 1, 0)
cnsa_tmp$median_fh <- ifelse(cnsa_tmp$cnsar_fh >= Quantile(cnsa_tmp$cnsar_fh, weights = cnsa_tmp$br_fh, probs = 0.5, na.rm = TRUE), 1, 0)
cnsa_tmp$median93_res <- ifelse(cnsa_tmp$cnsar93_res >= Quantile(cnsa_tmp$cnsar93_res, weights = cnsa_tmp$br93_res, probs = 0.5, na.rm = TRUE), 1, 0)

write.csv(cnsa_tmp, paste(repo, "r_workingdata/cnsa_all_st.csv", sep = ""))
# cnsa_tmp <- read.csv(paste(repo, "r_workingdata/cnsa_all_st.csv", sep = ""))

################################################################################
#### Merge CNSA data and ACS data
################################################################################
acs_sub <- acs_tmp %>%
  filter(age >= 19 & (birthyr %in% c(1994:2003)| (birthyr == 1993 & birthqtr %in% c(3:4)) )) %>% #only ones born after the first half of 1993
  rename(fips = bpl) %>%
  merge(cnsa_tmp, by = "fips", all.x = T) %>%
  mutate(qy = birthyr*10 + birthqtr, effective = ifelse(birthqtr == 4, birthyr, birthyr - 1))

################################################################################
#### Create confounders
################################################################################
# Confounder 1: BBA 1997 Medicaid expansion Data 1
# Data is from Hoynes and Luttmer (2011)
acs_sub1 <- read_dta("raw_data/confounds/MedicaidEligPreg.dta") %>%
  rename(fips = stfips) %>%
  merge(acs_sub, by = c("fips", "effective"), all.y = T)

# Confounder 2: PRWORA
# Data is from Schonei and Blank (2000)
# For ADFC, use time of impementation, which is already available
# For TANF, use time of impementation instead of time of authorization
acs_sub2 <- read.csv("raw_data/confounds/welfare_reform.csv") %>%
  mutate(waiver_imp_mo = waiver_imp_mo %>% as.numeric(),
         waiver_imp_yr = waiver_imp_yr %>% as.numeric(),
         tanf_imp_mo = ifelse(is.na(tanf_act_mo), tanf_ofc_mo, tanf_act_mo) %>% as.numeric(),
         tanf_imp_yr = ifelse(is.na(tanf_act_yr), tanf_ofc_yr, tanf_act_yr) %>% as.numeric()) %>%
  merge(acs_sub1, by = "fips", all.y = T) %>% #need birthqtr
  mutate(waiver_imp_my = (waiver_imp_yr + 1900)*100 + waiver_imp_mo,
         tanf_imp_my = (tanf_imp_yr + 1900)*100 + tanf_imp_mo,
         first_exp = case_when(
           birthqtr == 1 ~ (birthyr-1)* 100 + 8,
           birthqtr == 2 ~ (birthyr-1)* 100 + 11,
           birthqtr == 3 ~ birthyr* 100 + 2,
           birthqtr == 4 ~ birthyr* 100 + 5)) %>%
  mutate(waiver_exp = ifelse(first_exp >= waiver_imp_my, 1, 0) %>% replace(is.na(.), 0),
         tanf_exp = ifelse(first_exp >= tanf_imp_my, 1, 0) %>% replace(is.na(.), 0)) 

# Confounder 3: MHPA 1996
# from Buchmuller et al. (2007)
acs_sub3 <- read.csv("raw_data/confounds/mhlaws.csv") %>%
  merge(acs_sub2, by = "fips", all.y = T) %>%
  mutate(mh_my = mh_yr*100 + mh_mo) %>%
  mutate(mhlaws_exp = ifelse(first_exp >= mh_my, 1, 0) %>% replace(is.na(.), 0))

# Confounder 4: Unemployment rate as a proxy for the local economic factors at pregnancy. Bad controls, needs IV.
# Load annual local economic data from St.Louis Fed
fredr_set_key("5b9eaec48eb2797ef0a7238a09163045")
# State-year unemployment rate
unemp_st <- function(stfips) {
  fred_code <- ifelse(nchar(stfips) <2, paste("LAUST0",stfips, "0000000000003A", sep = ""), 
                      paste("LAUST",stfips, "0000000000003A", sep = ""))
  fredr(series_id = fred_code, observation_start = as.Date("1992-01-01")) %>%
    mutate(fips = stfips)
}
unemp <- purrr::map_dfr(cnsa_tmp$fips, unemp_st) %>%
  mutate(year = format(date, format = "%Y")) %>%
  select(fips, year, value) %>%
  rename(unemp = value, effective = year) 

emp_st_tmp <- read.csv("raw_data/confounds/efsy_panel_naics.csv") %>%
  filter(year %in% c(1987:2003) & emp!=0) %>%
  select(fipstate, naics12, emp, year) %>%
  group_by(fipstate, naics12, year) %>%
  summarise(emp = sum(emp))

#1. find the naics ended with "///", separate them out, delete subset
#2. sum the rest and merge the last one
emp_st_3d <- emp_st_tmp %>%
  filter(substr(naics12, 4, 6) == "///") %>%
  mutate(naics12 = substr(naics12, 1, 3))
emp_st_n3d <- emp_st_tmp %>%
  filter(!(substr(naics12, 1, 3) %in% substr(naics12, 1, 3)))
emp_st_3d2 <- emp_st_tmp %>%
  filter(substr(naics12, 4, 6) != "000" & year >= 1998) %>%
  mutate(naics12 = substr(naics12, 1, 3)) %>%
  group_by(fipstate, year, naics12) %>%
  summarise(emp = sum(emp)) %>%
  select(fipstate, naics12, year, emp) %>% 
  rbind(emp_st_3d)

summary(as.factor(emp_st_n3d$naics12))
# 11---- 21---- 22---- 23---- 31---- 42---- 44---- 48---- 51---- 52---- 53---- 54---- 55---- 56---- 61---- 62---- 71---- 72---- 81---- 
#    442    353    368    392    442    360    360    459    443    376    387    443    370    443    351    372    387    372    442 
# 92---- 95---- 99---- 
#    459    459    459 
# All industries with 4 or more digits have a 3-digit superset before 1998
# National employment in each industry
emp_us_tmp <- emp_st_3d2 %>%
  select(naics12, emp, year) %>%
  group_by(naics12, year) %>%
  summarise(emp = sum(emp)) %>%
  mutate(emp_lag1 = emp %>% lag(2)) %>%
  filter(year >= 1992) %>%
  mutate(emp_chg = log(emp) - log(emp_lag1))

emp_st_sum <- emp_st_3d2 %>%
  select(fipstate, naics12, emp, year) %>%
  group_by(fipstate, year) %>%
  summarise(emp_sum = sum(emp)) 

emp_st_3dm <- merge(emp_st_3d2, emp_st_sum, by = c("fipstate", "year")) %>%
  mutate(share = emp/emp_sum) 

emp_st_base <- emp_st_3dm %>%
  filter(year >= 1988 & year <= 1990 ) %>%
  select(fipstate, naics12, share) %>%
  group_by(fipstate, naics12) %>%
  summarise(share89_90 = mean(share))

emp_st_ssiv <- merge(emp_st_3dm, emp_st_base, by = c("fipstate", "naics12")) %>%
  filter(year >= 1992) %>% 
  merge(emp_us_tmp, by = c("year", "naics12")) %>%
  mutate(ss = share89_90*emp_chg) %>%
  select(fipstate, year, naics12, ss) %>%
  group_by(fipstate, year) %>%
  summarise(ssiv = sum(ss)) %>%
  rename(fips = fipstate, effective = year)

unemp_iv <- merge(unemp, emp_st_ssiv, by = c("fips", "effective")) %>%
  mutate(effective = effective %>% as.numeric())

cor(unemp_iv$unemp, unemp_iv$ssiv , method = c("pearson"))
var.test(unemp_iv$unemp, unemp_iv$ssiv, alternative = "two.sided")

write.csv(unemp_iv, paste(repo, "r_workingdata/unemp_iv.csv", sep = ""))
acs_sub4 <- merge(acs_sub3, unemp_iv, by = c("fips", "effective"), all.x = T)

################################################################################
#### Create other variables
################################################################################
acs_mg <- acs_sub4 %>%
  mutate(
    post = ifelse(effective >= 1996, 1, 0),
    race_rc = case_when(
      race == 1 ~ "White",
      race == 2 ~ "Black",
      TRUE ~ "Other"
    ),
    division = case_when(
      fips %in% c(9, 23, 25, 33, 44, 50) ~ 1,
      fips %in% c(34, 35, 42) ~ 2,
      fips %in% c(18, 17, 26, 39) ~ 3,
      fips %in% c(19, 20, 27, 29, 31, 38, 46) ~ 4,
      fips %in% c(10, 11, 12, 13, 24, 37, 45, 51, 54) ~ 5,
      fips %in% c(1, 21, 28, 47) ~ 6,
      fips %in% c(5, 22, 40, 48) ~ 7,
      fips %in% c(4, 8, 16, 35, 30, 49, 32, 56) ~ 8,
      fips %in% c(2, 6, 15, 41, 53) ~ 9
    ),
    hisp = ifelse( hispan !=0, 1, 0), # no missing value
    male = ifelse( sex ==1 , 1, 0), # no missing value
    stu = ifelse(school==2, 1, 0), # school only has values 1 and 2
    clg_stu = ifelse(gradeatt == 6, 1, 0 ),
    gs_stu = ifelse(gradeatt == 7, 1, 0),
    clg_gs = ifelse(gradeatt %in% c(6,7), 1, 0),
    hs_ged = ifelse(educd >= 63, 1, 0), # no missing data
    hs_grad = ifelse(educd >=63 & educd != 64, 1, 0),
    bchl = ifelse(educd >= 101, 1, 0),
    stem1 = ifelse(degfield %in% c(11, 13, 14, 20, 21, 24, 25, 36:40, 50, 51, 57:59, 61), 1, 0) %>%
      replace(degfieldd %in% c(1102,4001, 4007), 0),
    stem2 = ifelse(degfield2 %in% c(11, 13, 14, 20, 21, 24, 25, 36:40, 50, 51, 57:59, 61), 1, 0) %>%
      replace(degfield2d %in% c(1102,4001, 4007), 0),
    lf = ifelse(empstat != 3, 1, 0),
    employed = ifelse(empstat == 1, 1, 0),
    wklw = case_when(
      wrklstwk == 1 ~ 0,
      wrklstwk == 2 ~ 1,
      wrklstwk == 3 ~ NA_real_
    ),
    ft = ifelse(uhrswork >= 40, 1, 0) %>% replace(uhrswork == 0, NA),
    income = inctot/1000,
    wage = incwage/1000,
    earning = incearn/1000
  )

saveRDS(acs_mg, paste(repo, "r_workingdata/acs_mg.rds", sep = ""))
