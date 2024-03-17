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
library(DescTools) # weighted quantile
conflicted::conflicts_prefer(dplyr::lag)
conflicted::conflicts_prefer(dplyr::filter)

# Directories
setwd("C:/Users/wenji/OneDrive/github_files/GitHub/folate/")

fl_theme <- theme(legend.position = 'none', panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), axis.text = element_text(size = 12), 
                  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
lg_theme <- theme(legend.position = 'bottom', panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), axis.text = element_text(size = 12), 
                  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
                  legend.title = element_text(size = 12),
                  legend.text = element_text(size = 12))
################################################################################
#### Load the data
################################################################################
# Load the cleaned data
acs_mg <- readRDS("r_workingdata/acs_mg.rds")
cnsa_tmp <- read.csv("r_workingdata/cnsa_all_st.csv")
mtpler <- Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.75, na.rm = T) - Quantile(cnsa_tmp$cnsar_nc89_93_res, weights = cnsa_tmp$br89_93_res, probs = 0.25, na.rm = T)
mtpler_fh <- Quantile(cnsa_tmp$cnsar_fh, weights = cnsa_tmp$br_fh, probs = 0.75, na.rm = T) - Quantile(cnsa_tmp$cnsar_fh, weights = cnsa_tmp$br_fh, probs = 0.25, na.rm = T)
mtpler93 <- Quantile(cnsa_tmp$cnsar93_res, weights = cnsa_tmp$br93_res, probs = 0.75, na.rm = T) - Quantile(cnsa_tmp$cnsar93_res, weights = cnsa_tmp$br93_res, probs = 0.25, na.rm = T)
# mtpler = 0.498 case per 1000 births

# for etable, return %(Coefficient/Mean)
fitstat_register(type = "wt_mean", alias = "Dependent Variable Mean",
                 fun = function(x){
                   wt = weights(x) %>% as.data.frame() %>% filter(!is.na(.)) %>% rowid_to_column()
                   names(wt) = c("rowid", "wt")
                   y = model.matrix(x, type = "lhs") %>% as.data.frame() %>% rowid_to_column() %>%
                     merge(wt, by = "rowid", all.x = T)
                   names(y) = c("rowid", "y", "wt")
                   weighted.mean(y$y, y$wt/sum(y$wt), na.rm = T)} 
)

fitstat_register(type = "p_eff_ct", alias = "(Coefficient/Mean)",
                 fun = function(x) x$coefficients[["cnsar_nc89_93_res:post"]]*100*mtpler/fitstat(x, type = "wt_mean")[["wt_mean"]])
fitstat_register(type = "coef_rs", alias = "Coefficient Rescaled",
                 fun = function(x) x$coefficients[["cnsar_nc89_93_res:post"]]*mtpler)
fitstat_register(type = "se_rs", alias = "Std.Err Rescaled",
                 fun = function(x) x$se[["cnsar_nc89_93_res:post"]]*mtpler)

fitstat_register(type = "p_eff_mn", alias = "(Coefficient/Mean)",
                 fun = function(x) x$coefficients[["mean:post"]]*100/fitstat(x, type = "wt_mean")[["wt_mean"]])
fitstat_register(type = "p_eff_md", alias = "(Coefficient/Mean)",
                 fun = function(x) x$coefficients[["median:post"]]*100/fitstat(x, type = "wt_mean")[["wt_mean"]])
fitstat_register(type = "p_eff_40", alias = "(Coefficient/Mean)",
                 fun = function(x) x$coefficients[["top40:post"]]*100/fitstat(x, type = "wt_mean")[["wt_mean"]])
fitstat_register(type = "p_eff_30", alias = "(Coefficient/Mean)",
                 fun = function(x) x$coefficients[["top30:post"]]*100/fitstat(x, type = "wt_mean")[["wt_mean"]])

fitstat_register(type = "p_eff_fh", alias = "(Coefficient/Mean)",
                 fun = function(x) x$coefficients[["cnsar_fh:post"]]*100*mtpler_fh/fitstat(x, type = "wt_mean")[["wt_mean"]])
fitstat_register(type = "coef_fh_rs", alias = "Coefficient Rescaled",
                 fun = function(x) x$coefficients[["cnsar_fh:post"]]*mtpler_fh)
fitstat_register(type = "se_fh_rs", alias = "Std.Err Rescaled",
                 fun = function(x) x$se[["cnsar_fh:post"]]*mtpler_fh)
fitstat_register(type = "p_eff_mdfh", alias = "(Coefficient/Mean)",
                 fun = function(x) x$coefficients[["median_fh:post"]]*100/fitstat(x, type = "wt_mean")[["wt_mean"]])

fitstat_register(type = "p_eff_93", alias = "(Coefficient/Mean)",
                 fun = function(x) x$coefficients[["cnsar93_res:post"]]*100*mtpler93/fitstat(x, type = "wt_mean")[["wt_mean"]])
fitstat_register(type = "p_eff_md93", alias = "(Coefficient/Mean)",
                 fun = function(x) x$coefficients[["median93_res:post"]]*100/fitstat(x, type = "wt_mean")[["wt_mean"]])
fitstat_register(type = "coef93_rs", alias = "Coefficient Rescaled",
                 fun = function(x) x$coefficients[["cnsar93_res:post"]]*mtpler93)
fitstat_register(type = "se93_rs", alias = "Std.Err Rescaled",
                 fun = function(x) x$se[["cnsar93_res:post"]]*mtpler93)

acs_mg <- acs_mg %>%
  mutate(
    partial = ifelse(effective == 1996, 1, 0),
    full = ifelse(effective >=  1997, 1, 0)
  )

early_dta <- acs_mg %>%
  filter(age %in% c(19:22))
late_dta <- acs_mg %>%
  filter(age %in% c(23:29) )

dd_clg_stu3 <- feols(gs_stu ~ cnsar_nc89_93_res:partial + cnsar_nc89_93_res:full
                     + factor(statefip):factor(year) + factor(birthyr):factor(year)
                     + threshpreg + mhlaws_exp + waiver_exp + tanf_exp
                     + factor(race_rc) + hisp + male  + ssiv
                     |fips + qy,
                     cluster = ~fips,
                     weights = ~ perwt,
                     data = acs_mg[acs_mg$age >= 23, ])
summary(dd_clg_stu3)

dd_clg_stu3 <- feols(clg_gs ~ cnsar_nc89_93_res:partial + cnsar_nc89_93_res:full
                     + factor(statefip):factor(year) + factor(birthyr):factor(year)
                     + threshpreg + mhlaws_exp + waiver_exp + tanf_exp
                     + factor(race_rc) + hisp + male  + ssiv
                     |fips + qy,
                     vcov = "hetero",
                     weights = ~ perwt,
                     data = acs_mg)
summary(dd_clg_stu3)

dd_clg_stu3 <- feols(clg_gs ~ cnsar_nc89_93_res:partial + cnsar_nc89_93_res:full
                     + factor(statefip):factor(year)
                     + threshpreg + mhlaws_exp + waiver_exp + tanf_exp
                     + factor(race) + hisp + male  + ssiv
                     |fips + qy,
                     cluster = ~fips,
                     weights = ~ perwt,
                     data = acs_mg)
summary(dd_clg_stu3)

dd_clg_stu3 <- feols(clg_gs ~ cnsar_nc89_93_res:post
                     + factor(puma) + factor(year)
                     + threshpreg + mhlaws_exp + waiver_exp + tanf_exp
                     + factor(race) + hisp + male  + ssiv
                     | fips + qy,
                     cluster = ~fips,
                     # vcov = "hetero",
                     weights = ~ perwt,
                     data = acs_mg[acs_mg$year %in% c(2016:2019), ])
summary(dd_clg_stu3)


################################################################################
#### 1. Educational outcomes: DD
################################################################################
# 1.1 College enrollment for 19-24 y.o.
dd_clg_stu1 <- feols(clg_stu ~ cnsar_nc89_93_res:post
                     | fips + effective,
                     cluster = ~fips, 
                     weights = ~ perwt,
                     data = early_dta)
summary(dd_clg_stu1)

dd_clg_stu2 <- feols(clg_stu ~ cnsar_nc89_93_res:post
                     | fips + qy,
                     cluster = ~fips, 
                     weights = ~ perwt,
                     data = early_dta)
summary(dd_clg_stu2)

dd_clg_stu3 <- feols(clg_stu ~ cnsar_nc89_93_res:post
                     + factor(statefip):factor(year) 
                     + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                     + factor(race) + hisp  + ssiv
                     | fips + qy,
                     cluster = ~fips, 
                     weights = ~ perwt,
                     data = early_dta[early_dta$male==0, ])
summary(dd_clg_stu3)

dd_clg_stu4 <- feols(clg_stu ~ cnsar_nc89_93_res:post
                     + factor(statefip):factor(year) 
                     + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                     + factor(race) + hisp  + male + ssiv
                     + factor(fips):qy
                     | qy,
                     cluster = ~fips, 
                     weights = ~ perwt,
                     data = early_dta)
summary(dd_clg_stu4)

# 1.2 Graduate school enrollment for 19-24 y.o.
dd_gs_stu1 <- feols(gs_stu ~ cnsar_nc89_93_res:post
                | fips + effective,
                cluster = ~fips, 
                weights = ~ perwt,
                data = late_dta)
summary(dd_gs_stu1)

dd_gs_stu2 <- feols(gs_stu ~ cnsar_nc89_93_res:post
                | fips + qy,
                cluster = ~fips, 
                weights = ~ perwt,
                data = late_dta)
summary(dd_gs_stu2)

dd_gs_stu3 <- feols(gs_stu ~ cnsar_nc89_93_res:post
                + factor(statefip):factor(year) 
                + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                + factor(race) + hisp  + male + ssiv
                | fips + qy,
                cluster = ~fips, 
                weights = ~ perwt,
                data = late_dta)
summary(dd_gs_stu3)

dd_gs_stu4 <- feols(gs_stu ~ cnsar_nc89_93_res:post
                + factor(statefip):factor(year) 
                + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                + factor(race) + hisp  + male  + ssiv
                + factor(fips):qy
                | qy,
                cluster = ~fips, 
                weights = ~ perwt,
                data = late_dta)
summary(dd_gs_stu4)


etable(dd_clg_stu1, dd_clg_stu2,dd_clg_stu3,dd_clg_stu4, dd_gs_stu1, dd_gs_stu2,dd_gs_stu3,dd_gs_stu4,
       keep = "%cnsar_nc89_93_res:post",
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~  r2 + n + wt_mean + p_eff_ct + coef_rs + se_rs, tex = TRUE)


# 2.1 High school diploma/GED for all >= 19 years old
dd_hs_ged1 <- feols(hs_ged ~ cnsar_nc89_93_res:post
                    | fips + effective,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = acs_mg)
summary(dd_hs_ged1)

dd_hs_ged2 <- feols(hs_ged ~ cnsar_nc89_93_res:post
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = acs_mg)
summary(dd_hs_ged2)

dd_hs_ged3 <- feols(hs_ged ~ cnsar_nc89_93_res:post
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = acs_mg)
summary(dd_hs_ged3)

dd_hs_ged4 <- feols(hs_ged ~ cnsar_nc89_93_res:post
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    + factor(fips):qy
                    | qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = acs_mg)
summary(dd_hs_ged4)

# 2.1 Bachelor's degree for 23-29 years old
dd_bchl1 <- feols(bchl ~ cnsar_nc89_93_res:post
                  | fips + effective,
                  cluster = ~fips, 
                  weights = ~ perwt,
                  data = late_dta)
summary(dd_bchl1)

dd_bchl2 <- feols(bchl ~ cnsar_nc89_93_res:post
                  | fips + qy,
                  cluster = ~fips, 
                  weights = ~ perwt,
                  data = late_dta)
summary(dd_bchl2)

dd_bchl3 <- feols(bchl ~ cnsar_nc89_93_res:post
                  + factor(statefip):factor(year) 
                  + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                  + factor(race) + hisp  + male + ssiv
                  | fips + qy,
                  cluster = ~fips, 
                  weights = ~ perwt,
                  data = late_dta)
summary(dd_bchl3)

dd_bchl4 <- feols(bchl ~ cnsar_nc89_93_res:post
                  + factor(statefip):factor(year) 
                  + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                  + factor(race) + hisp  + male  + ssiv
                  + factor(fips):qy
                  | qy,
                  cluster = ~fips, 
                  weights = ~ perwt,
                  data = late_dta)
summary(dd_bchl4)


etable(dd_hs_ged1, dd_hs_ged2,dd_hs_ged3,dd_hs_ged4,dd_bchl1,dd_bchl2, dd_bchl3,dd_bchl4,
       keep = "%cnsar_nc89_93_res:post",
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~ r2 + n + wt_mean + p_eff_ct + coef_rs + se_rs , tex = TRUE)

# dd_hs_ged4[["coefficients"]][["cnsar_nc89_93_res:post"]]*mtpler
# dd_hs_ged4[["se"]][["cnsar_nc89_93_res:post"]]*mtpler
# 
# dd_bchl4[["coefficients"]][["cnsar_nc89_93_res:post"]]*mtpler
# dd_bchl4[["se"]][["cnsar_nc89_93_res:post"]]*mtpler
# 
# dlt_x <- data.frame(row_dlt = x[["obs_selection"]][["obsRemoved"]]*(-1))
# x_dta <- wic_dta[wic_dta$births >= 25 & wic_dta$group == "LEUM" ,]  %>%
#   rowid_to_column() %>%
#   filter(!(rowid %in% dlt_twfe8$row_dlt))
# mean(nobs_twfe8$wic_par, na.rm = T)
# nrow(nobs_twfe8)

################################################################################
#### 2. Educational outcomes: event study
################################################################################
# 2.1 Event study plot for college enrollment

es_clg_stu <- feols(clg_gs ~ i(effective, cnsar_nc89_93_res, ref = 1995) 
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + male +  ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = acs_mg)

esreg <- data.frame(coef =c(
  es_clg_stu$coefficients[["effective::1992:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$coefficients[["effective::1993:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$coefficients[["effective::1994:cnsar_nc89_93_res"]]*mtpler,
  0,
  es_clg_stu$coefficients[["effective::1996:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$coefficients[["effective::1997:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$coefficients[["effective::1998:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$coefficients[["effective::1999:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$coefficients[["effective::2000:cnsar_nc89_93_res"]]*mtpler
), 
se = c(
  es_clg_stu$se[["effective::1992:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$se[["effective::1993:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$se[["effective::1994:cnsar_nc89_93_res"]]*mtpler,
  0,
  es_clg_stu$se[["effective::1996:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$se[["effective::1997:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$se[["effective::1998:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$se[["effective::1999:cnsar_nc89_93_res"]]*mtpler,
  es_clg_stu$se[["effective::2000:cnsar_nc89_93_res"]]*mtpler
), 
year =c('1992',
        '1993',
        '1994',
        '1995',
        '1996',
        '1997',
        '1998',
        '1999',
        '2000') %>%
  factor(
    levels = c(
      '1992',
      '1993',
      '1994',
      '1995',
      '1996',
      '1997',
      '1998',
      '1999',
      '2000'
    ),
    ordered = T
  ))

esreg$group <- case_when(
  esreg$year < 1996  ~ 'Pre',
  esreg$year >= 1996 ~ 'Post')
write.csv(esreg, paste(dict_wd, "es_clg_stu_coef.csv", sep = ""))

ggplot(data = esreg, aes(y = coef, x = year, color = group)) +
  geom_point(size = 3) +
  geom_line(linewidth =0.5, group = 1) +
  geom_errorbar(aes(ymin = (coef - 1.96 * se), ymax = (coef + 1.96 * se), color = group), width =0.2) +
  scale_color_manual(values = c("#144062", "#168aad")) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = as.numeric(esreg$year)[5] - 0.5, linetype = 2)+
  # scale_y_continuous(limits = c(-.06, .06), breaks = seq(-.06, .06, by = 0.02)) +
  xlab("Effective Year of Exposure") +
  ylab("Coefficients") +
  theme_classic()+
  fl_theme

ggsave("es_clg_stu_res.png", path = "r_outputs/", width = 6.472, height = 4.2)


# 2.2 Event study plot for college enrollment
es_gs_stu <- feols(gs_stu ~ i(effective, cnsar_nc89_93_res, ref = 1995) 
               + factor(statefip):factor(year)
               + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
               + factor(race) + hisp  + male + ssiv
               | fips + qy,
               cluster = ~fips,
               # vcov = "hetero",
               weights = ~ perwt,
               data = acs_mg[acs_mg$age >= 23 & acs_mg$year %in% c(2021:2022) , ])

esreg <- data.frame(coef =c(
  es_gs_stu$coefficients[["effective::1992:cnsar_nc89_93_res"]]*mtpler,
  es_gs_stu$coefficients[["effective::1993:cnsar_nc89_93_res"]]*mtpler,
  es_gs_stu$coefficients[["effective::1994:cnsar_nc89_93_res"]]*mtpler,
  0,
  es_gs_stu$coefficients[["effective::1996:cnsar_nc89_93_res"]]*mtpler,
  es_gs_stu$coefficients[["effective::1997:cnsar_nc89_93_res"]]*mtpler,
  es_gs_stu$coefficients[["effective::1998:cnsar_nc89_93_res"]]*mtpler,
  es_gs_stu$coefficients[["effective::1999:cnsar_nc89_93_res"]]
), 
se = c(
  es_gs_stu$se[["effective::1992:cnsar_nc89_93_res"]]*mtpler,
  es_gs_stu$se[["effective::1993:cnsar_nc89_93_res"]]*mtpler,
  es_gs_stu$se[["effective::1994:cnsar_nc89_93_res"]]*mtpler,
  0,
  es_gs_stu$se[["effective::1996:cnsar_nc89_93_res"]]*mtpler,
  es_gs_stu$se[["effective::1997:cnsar_nc89_93_res"]]*mtpler,
  es_gs_stu$se[["effective::1998:cnsar_nc89_93_res"]]*mtpler,
  es_gs_stu$se[["effective::1999:cnsar_nc89_93_res"]]*mtpler
), 
year =c('1992',
        '1993',
        '1994',
        '1995',
        '1996',
        '1997',
        '1998',
        '1999') %>%
  factor(
    levels = c(
      '1992',
      '1993',
      '1994',
      '1995',
      '1996',
      '1997',
      '1998',
      '1999'
    ),
    ordered = T
  ))

esreg$group <- case_when(
  esreg$year < 1996  ~ 'Pre',
  esreg$year >= 1996 ~ 'Post')
write.csv(esreg, paste(dict_wd, "es_gs_stu_coef.csv", sep = ""))

ggplot(data = esreg, aes(y = coef, x = year)) +
  geom_point(aes(color = group), size = 3) +
  # geom_line(linewidth =0.5, group = 1) +
  geom_ribbon(aes(ymin = (coef - 1.96 * se), ymax = (coef + 1.96 * se), group = 1), alpha =0.1) +
  # geom_errorbar(aes(ymin = (coef - 1.96 * se), ymax = (coef + 1.96 * se), color = group), width =0.2) +
  scale_color_manual(values = c("#144062", "#168aad")) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = as.numeric(esreg$year)[5] - 0.5, linetype = 2)+
  # scale_y_continuous(limits = c(-.02, .02), breaks = seq(-.02, .02, by = .005)) +
  xlab("Effective Year of Exposure") +
  ylab("Coefficients") +
  theme_classic()+
  fl_theme

ggsave("es_gs_stu_ipums.png", path = "r_outputs/", width = 6.472, height = 4.2)

# 2.3 Event study plot for HS Diploma/GED
es_hs_ged <- feols(clg_stu ~ i(effective, cnsar_nc89_93_res, ref = 1995) 
                   + factor(statefip):factor(year) 
                   + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                   + factor(race_rc) + hisp  + male + ssiv
                   | fips + qy,
                   cluster = ~fips,
                   weights = ~ perwt,
                   data = acs_mg)

summary(es_hs_ged)

esreg <- data.frame(coef =c(
  es_hs_ged$coefficients[["effective::1992:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$coefficients[["effective::1993:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$coefficients[["effective::1994:cnsar_nc89_93_res"]]*mtpler,
  0,
  es_hs_ged$coefficients[["effective::1996:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$coefficients[["effective::1997:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$coefficients[["effective::1998:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$coefficients[["effective::1999:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$coefficients[["effective::2000:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$coefficients[["effective::2001:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$coefficients[["effective::2002:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$coefficients[["effective::2003:cnsar_nc89_93_res"]]*mtpler
), 
se = c(
  es_hs_ged$se[["effective::1992:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$se[["effective::1993:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$se[["effective::1994:cnsar_nc89_93_res"]]*mtpler,
  0,
  es_hs_ged$se[["effective::1996:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$se[["effective::1997:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$se[["effective::1998:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$se[["effective::1999:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$se[["effective::2000:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$se[["effective::2001:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$se[["effective::2002:cnsar_nc89_93_res"]]*mtpler,
  es_hs_ged$se[["effective::2003:cnsar_nc89_93_res"]]*mtpler
), 
year =c('1992',
        '1993',
        '1994',
        '1995',
        '1996',
        '1997',
        '1998',
        '1999',
        '2000',
        '2001',
        '2002',
        '2003') %>%
  factor(
    levels = c(
      '1992',
      '1993',
      '1994',
      '1995',
      '1996',
      '1997',
      '1998',
      '1999',
      '2000',
      '2001',
      '2002',
      '2003'
    ),
    ordered = T
  ))

esreg$group <- case_when(
  esreg$year < 1996  ~ 'Pre',
  esreg$year >= 1996 ~ 'Post')
write.csv(esreg, paste(dict_wd, "es_hs_ged_coef.csv", sep = ""))

ggplot(data = esreg, aes(y = coef, x = year)) +
  geom_point(aes(color = group),  size = 3) +
  # geom_line(linewidth =0.5) +
  geom_ribbon(aes(ymin = (coef - 1.96 * se), ymax = (coef + 1.96 * se), group = 1), alpha =0.1) +
  scale_color_manual(values = c("#144062", "#168aad")) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = as.numeric(esreg$year)[5] - 0.5, linetype = 2)+
  # scale_y_continuous(limits = c(-.02, .02), breaks = seq(-.02, .02, by = 0.005)) +
  xlab("Effective Year of Exposure") +
  ylab("Coefficients") +
  theme_classic()+
  fl_theme

ggsave("es_clg_gs_ipums.png", path = "r_outputs/", width = 6.472, height = 4.2)


# 2.4 Event study plot for bachelor's degree
es_bchl <- feols(bchl ~ i(effective, cnsar_nc89_93_res, ref = 1995) 
                 + factor(statefip):factor(year) 
                 + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                 + factor(race) + hisp  + male + ssiv
                 | fips + qy,
                 cluster = ~fips, 
                 weights = ~ perwt,
                 data = late_dta)

esreg <- data.frame(coef =c(
  es_bchl$coefficients[["effective::1992:cnsar_nc89_93_res"]]*mtpler,
  es_bchl$coefficients[["effective::1993:cnsar_nc89_93_res"]]*mtpler,
  es_bchl$coefficients[["effective::1994:cnsar_nc89_93_res"]]*mtpler,
  0,
  es_bchl$coefficients[["effective::1996:cnsar_nc89_93_res"]]*mtpler,
  es_bchl$coefficients[["effective::1997:cnsar_nc89_93_res"]]*mtpler,
  es_bchl$coefficients[["effective::1998:cnsar_nc89_93_res"]]*mtpler,
  es_bchl$coefficients[["effective::1999:cnsar_nc89_93_res"]]
), 
se = c(
  es_bchl$se[["effective::1992:cnsar_nc89_93_res"]]*mtpler,
  es_bchl$se[["effective::1993:cnsar_nc89_93_res"]]*mtpler,
  es_bchl$se[["effective::1994:cnsar_nc89_93_res"]]*mtpler,
  0,
  es_bchl$se[["effective::1996:cnsar_nc89_93_res"]]*mtpler,
  es_bchl$se[["effective::1997:cnsar_nc89_93_res"]]*mtpler,
  es_bchl$se[["effective::1998:cnsar_nc89_93_res"]]*mtpler,
  es_bchl$se[["effective::1999:cnsar_nc89_93_res"]]*mtpler
), 
year =c('1992',
        '1993',
        '1994',
        '1995',
        '1996',
        '1997',
        '1998',
        '1999') %>%
  factor(
    levels = c(
      '1992',
      '1993',
      '1994',
      '1995',
      '1996',
      '1997',
      '1998',
      '1999'
    ),
    ordered = T
  ))

esreg$group <- case_when(
  esreg$year < 1996  ~ 'Pre',
  esreg$year >= 1996 ~ 'Post')
write.csv(esreg, paste(dict_wd, "es_bchl_coef.csv", sep = ""))

ggplot(data = esreg, aes(y = coef, x = year)) +
  geom_point(aes(color = group), size = 3) +
  geom_ribbon(aes(ymin = (coef - 1.96 * se), ymax = (coef + 1.96 * se), group = 1), alpha =0.1) +
  # geom_errorbar(aes(ymin = (coef - 1.96 * se), ymax = (coef + 1.96 * se), color = group), width =0.2) +
  scale_color_manual(values = c("#144062", "#168aad")) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = as.numeric(esreg$year)[5] - 0.5, linetype = 2)+
  scale_y_continuous(limits = c(-.065, .065), breaks = seq(-.06, .06, by = 0.02)) +
  xlab("Effective Year of Exposure") +
  ylab("Coefficients") +
  theme_classic()+
  fl_theme

ggsave("es_bchl_ipums.png", path = "r_outputs/", width = 6.472, height = 4.2)

################################################################################
#### 3. Educational outcomes: high and low (bin), DD
################################################################################
# 3.0 >= Mean
dd_clg_stu_bin0 <- feols(clg_stu ~ mean:post
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)
summary(dd_clg_stu_bin0)

dd_gs_stu_bin0 <- feols(gs_stu ~ mean:post
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)
summary(dd_gs_stu_bin0)

dd_hs_ged_bin0 <- feols(hs_ged ~ mean:post
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)
summary(dd_hs_ged_bin0)

dd_bchl_bin0 <- feols(bchl ~ mean:post
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)
summary(dd_bchl_bin0)

etable(dd_clg_stu_bin0,dd_gs_stu_bin0,dd_hs_ged_bin0, dd_bchl_bin0,
       keep = c("%mean:post"),
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~ r2 + n + wt_mean + p_eff_mn, tex = TRUE)

# 3.1 Top 50th percentile
dd_clg_stu_bin1 <- feols(clg_stu ~ median:post
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)
summary(dd_clg_stu_bin1)

dd_gs_stu_bin1 <- feols(gs_stu ~ median:post
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)
summary(dd_gs_stu_bin1)

dd_hs_ged_bin1 <- feols(hs_ged ~ median:post
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)
summary(dd_hs_ged_bin1)

dd_bchl_bin1 <- feols(bchl ~ median:post
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)
summary(dd_bchl_bin1)

etable(dd_clg_stu_bin1,dd_gs_stu_bin1,dd_hs_ged_bin1, dd_bchl_bin1,
       keep = c("%median:post"),
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~ r2 + n + wt_mean + p_eff_md, tex = TRUE)

# 3.2 Top 40th percentile
dd_clg_stu_bin2 <- feols(clg_stu ~ top40:post
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)
summary(dd_clg_stu_bin2)

dd_gs_stu_bin2 <- feols(gs_stu ~ top40:post
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)
summary(dd_gs_stu_bin2)

dd_hs_ged_bin2 <- feols(hs_ged ~ top40:post
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)
summary(dd_hs_ged_bin2)

dd_bchl_bin2 <- feols(bchl ~ top40:post
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)
summary(dd_bchl_bin2)

etable(dd_clg_stu_bin2,dd_gs_stu_bin2,dd_hs_ged_bin2, dd_bchl_bin2,
       keep = c("%top40:post"),
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~ r2 + n + wt_mean + p_eff_40, tex = TRUE)

#3.3 Top 30th percentile as high CNA anomaly
dd_clg_stu_bin3 <- feols(clg_stu ~ top30:post
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)
summary(dd_clg_stu_bin3)

dd_gs_stu_bin3 <- feols(gs_stu ~ top30:post
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)
summary(dd_gs_stu_bin3)

dd_hs_ged_bin3 <- feols(hs_ged ~ top30:post
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)
summary(dd_hs_ged_bin3)

dd_bchl_bin3 <- feols(bchl ~ top30:post
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)
summary(dd_bchl_bin3)

etable(dd_clg_stu_bin3,dd_gs_stu_bin3,dd_hs_ged_bin3, dd_bchl_bin3,
       keep = c("%top30:post"),
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~ r2 + n + wt_mean + p_eff_30, tex = TRUE)

################################################################################
#### 4. Educational outcomes: high and low (bin), event study
################################################################################
# 4.1 College enrollment
es_clg_stu_bin1 <- feols(clg_stu ~ i(effective, mean, ref = 1995) 
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)

es_clg_stu_bin2 <- feols(clg_stu ~ i(effective, median, ref = 1995)
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)

es_clg_stu_bin3 <- feols(clg_stu ~ i(effective, top40, ref = 1995)
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)

es_clg_stu_bin4 <- feols(clg_stu ~ i(effective, top30, ref = 1995)
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)

rb1_df <- data.frame(coef =c(
  es_clg_stu_bin1$coefficients[["effective::1992:mean"]],
  es_clg_stu_bin1$coefficients[["effective::1993:mean"]],
  es_clg_stu_bin1$coefficients[["effective::1994:mean"]],
  0,
  es_clg_stu_bin1$coefficients[["effective::1996:mean"]],
  es_clg_stu_bin1$coefficients[["effective::1997:mean"]],
  es_clg_stu_bin1$coefficients[["effective::1998:mean"]],
  es_clg_stu_bin1$coefficients[["effective::1999:mean"]],
  es_clg_stu_bin1$coefficients[["effective::2000:mean"]]
), 
se = c(
  es_clg_stu_bin1$se[["effective::1992:mean"]],
  es_clg_stu_bin1$se[["effective::1993:mean"]],
  es_clg_stu_bin1$se[["effective::1994:mean"]],
  0,
  es_clg_stu_bin1$se[["effective::1996:mean"]],
  es_clg_stu_bin1$se[["effective::1997:mean"]],
  es_clg_stu_bin1$se[["effective::1998:mean"]],
  es_clg_stu_bin1$se[["effective::1999:mean"]],
  es_clg_stu_bin1$se[["effective::2000:mean"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000'),
    ordered = TRUE
  ),
group = ">= Mean", position = 1)

rb2_df <- data.frame(coef =c(
  es_clg_stu_bin2$coefficients[["effective::1992:median"]],
  es_clg_stu_bin2$coefficients[["effective::1993:median"]],
  es_clg_stu_bin2$coefficients[["effective::1994:median"]],
  0,
  es_clg_stu_bin2$coefficients[["effective::1996:median"]],
  es_clg_stu_bin2$coefficients[["effective::1997:median"]],
  es_clg_stu_bin2$coefficients[["effective::1998:median"]],
  es_clg_stu_bin2$coefficients[["effective::1999:median"]],
  es_clg_stu_bin2$coefficients[["effective::2000:median"]]
), 
se = c(
  es_clg_stu_bin2$se[["effective::1992:median"]],
  es_clg_stu_bin2$se[["effective::1993:median"]],
  es_clg_stu_bin2$se[["effective::1994:median"]],
  0,
  es_clg_stu_bin2$se[["effective::1996:median"]],
  es_clg_stu_bin2$se[["effective::1997:median"]],
  es_clg_stu_bin2$se[["effective::1998:median"]],
  es_clg_stu_bin2$se[["effective::1999:median"]],
  es_clg_stu_bin2$se[["effective::2000:median"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000'),
    ordered = TRUE
  ),
group = ">= Median", position = 2)

rb3_df <- data.frame(coef =c(
  es_clg_stu_bin3$coefficients[["effective::1992:top40"]],
  es_clg_stu_bin3$coefficients[["effective::1993:top40"]],
  es_clg_stu_bin3$coefficients[["effective::1994:top40"]],
  0,
  es_clg_stu_bin3$coefficients[["effective::1996:top40"]],
  es_clg_stu_bin3$coefficients[["effective::1997:top40"]],
  es_clg_stu_bin3$coefficients[["effective::1998:top40"]],
  es_clg_stu_bin3$coefficients[["effective::1999:top40"]],
  es_clg_stu_bin3$coefficients[["effective::2000:top40"]]
), 
se = c(
  es_clg_stu_bin3$se[["effective::1992:top40"]],
  es_clg_stu_bin3$se[["effective::1993:top40"]],
  es_clg_stu_bin3$se[["effective::1994:top40"]],
  0,
  es_clg_stu_bin3$se[["effective::1996:top40"]],
  es_clg_stu_bin3$se[["effective::1997:top40"]],
  es_clg_stu_bin3$se[["effective::1998:top40"]],
  es_clg_stu_bin3$se[["effective::1999:top40"]],
  es_clg_stu_bin3$se[["effective::2000:top40"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000'),
    ordered = TRUE
  ),
group = "Top 40 Percentile", position = 3)

rb4_df <- data.frame(coef =c(
  es_clg_stu_bin4$coefficients[["effective::1992:top30"]],
  es_clg_stu_bin4$coefficients[["effective::1993:top30"]],
  es_clg_stu_bin4$coefficients[["effective::1994:top30"]],
  0,
  es_clg_stu_bin4$coefficients[["effective::1996:top30"]],
  es_clg_stu_bin4$coefficients[["effective::1997:top30"]],
  es_clg_stu_bin4$coefficients[["effective::1998:top30"]],
  es_clg_stu_bin4$coefficients[["effective::1999:top30"]],
  es_clg_stu_bin4$coefficients[["effective::2000:top30"]]
), 
se = c(
  es_clg_stu_bin4$se[["effective::1992:top30"]],
  es_clg_stu_bin4$se[["effective::1993:top30"]],
  es_clg_stu_bin4$se[["effective::1994:top30"]],
  0,
  es_clg_stu_bin4$se[["effective::1996:top30"]],
  es_clg_stu_bin4$se[["effective::1997:top30"]],
  es_clg_stu_bin4$se[["effective::1998:top30"]],
  es_clg_stu_bin4$se[["effective::1999:top30"]],
  es_clg_stu_bin4$se[["effective::2000:top30"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000'),
    ordered = TRUE
  ),
group = "Top 30 Percentile", position = 4)

rb <- data.table::rbindlist(list(rb1_df, rb2_df, rb3_df, rb4_df))

rb$group <- rb$group %>%
  factor(levels = c(">= Mean", ">= Median", "Top 40 Percentile", "Top 30 Percentile"))
write.csv(rb, paste(dict_wd, "es_clg_stu_bin.csv", sep = ""))

ggplot(data = rb, aes(y = coef, x = yr, group = group, color = group)) +
  geom_point(aes(color = group), size = 1.5, position = position_nudge(x = (rb$position - 2.5)*0.12 )) +
  geom_errorbar(aes(ymin = (coef- 1.96 * se), ymax = (coef + 1.96 * se), color = group), width =0.1, 
                position = position_nudge(x = (rb$position - 2.5)*0.12))+
  scale_color_manual(values = c("#144062", "#1e6091", "#168aad", "#76c893")) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = as.numeric(rb$yr)[5] - 0.5, linetype = 2)+
  scale_y_continuous(limits = c(-.08, .08), breaks = seq(-.08, .08, by = 0.02)) +
  xlab("Effective Year of Exposure") +
  ylab("Coefficients") +
  theme_classic()+
  lg_theme

ggsave("es_clg_stu_bin.png", path = "r_outputs/", width = 6.472, height = 4.5 )

# 4.2 Graduate/Professional school enrollment
es_gs_stu_bin1 <- feols(gs_stu ~ i(effective, mean, ref = 1995) 
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)

es_gs_stu_bin2 <- feols(gs_stu ~ i(effective, median, ref = 1995)
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)

es_gs_stu_bin3 <- feols(gs_stu ~ i(effective, top40, ref = 1995)
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)

es_gs_stu_bin4 <- feols(gs_stu ~ i(effective, top30, ref = 1995)
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)

rb1_df <- data.frame(coef =c(
  es_gs_stu_bin1$coefficients[["effective::1992:mean"]],
  es_gs_stu_bin1$coefficients[["effective::1993:mean"]],
  es_gs_stu_bin1$coefficients[["effective::1994:mean"]],
  0,
  es_gs_stu_bin1$coefficients[["effective::1996:mean"]],
  es_gs_stu_bin1$coefficients[["effective::1997:mean"]],
  es_gs_stu_bin1$coefficients[["effective::1998:mean"]],
  es_gs_stu_bin1$coefficients[["effective::1999:mean"]]
), 
se = c(
  es_gs_stu_bin1$se[["effective::1992:mean"]],
  es_gs_stu_bin1$se[["effective::1993:mean"]],
  es_gs_stu_bin1$se[["effective::1994:mean"]],
  0,
  es_gs_stu_bin1$se[["effective::1996:mean"]],
  es_gs_stu_bin1$se[["effective::1997:mean"]],
  es_gs_stu_bin1$se[["effective::1998:mean"]],
  es_gs_stu_bin1$se[["effective::1999:mean"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999'),
    ordered = TRUE
  ),
group = ">= Mean", position = 1)

rb2_df <- data.frame(coef =c(
  es_gs_stu_bin2$coefficients[["effective::1992:median"]],
  es_gs_stu_bin2$coefficients[["effective::1993:median"]],
  es_gs_stu_bin2$coefficients[["effective::1994:median"]],
  0,
  es_gs_stu_bin2$coefficients[["effective::1996:median"]],
  es_gs_stu_bin2$coefficients[["effective::1997:median"]],
  es_gs_stu_bin2$coefficients[["effective::1998:median"]],
  es_gs_stu_bin2$coefficients[["effective::1999:median"]]
), 
se = c(
  es_gs_stu_bin2$se[["effective::1992:median"]],
  es_gs_stu_bin2$se[["effective::1993:median"]],
  es_gs_stu_bin2$se[["effective::1994:median"]],
  0,
  es_gs_stu_bin2$se[["effective::1996:median"]],
  es_gs_stu_bin2$se[["effective::1997:median"]],
  es_gs_stu_bin2$se[["effective::1998:median"]],
  es_gs_stu_bin2$se[["effective::1999:median"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999'),
    ordered = TRUE
  ),
group = ">= Median", position = 2)

rb3_df <- data.frame(coef =c(
  es_gs_stu_bin3$coefficients[["effective::1992:top40"]],
  es_gs_stu_bin3$coefficients[["effective::1993:top40"]],
  es_gs_stu_bin3$coefficients[["effective::1994:top40"]],
  0,
  es_gs_stu_bin3$coefficients[["effective::1996:top40"]],
  es_gs_stu_bin3$coefficients[["effective::1997:top40"]],
  es_gs_stu_bin3$coefficients[["effective::1998:top40"]],
  es_gs_stu_bin3$coefficients[["effective::1999:top40"]]
), 
se = c(
  es_gs_stu_bin3$se[["effective::1992:top40"]],
  es_gs_stu_bin3$se[["effective::1993:top40"]],
  es_gs_stu_bin3$se[["effective::1994:top40"]],
  0,
  es_gs_stu_bin3$se[["effective::1996:top40"]],
  es_gs_stu_bin3$se[["effective::1997:top40"]],
  es_gs_stu_bin3$se[["effective::1998:top40"]],
  es_gs_stu_bin3$se[["effective::1999:top40"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999'),
    ordered = TRUE
  ),
group = "Top 40 Percentile", position = 3)

rb4_df <- data.frame(coef =c(
  es_gs_stu_bin4$coefficients[["effective::1992:top30"]],
  es_gs_stu_bin4$coefficients[["effective::1993:top30"]],
  es_gs_stu_bin4$coefficients[["effective::1994:top30"]],
  0,
  es_gs_stu_bin4$coefficients[["effective::1996:top30"]],
  es_gs_stu_bin4$coefficients[["effective::1997:top30"]],
  es_gs_stu_bin4$coefficients[["effective::1998:top30"]],
  es_gs_stu_bin4$coefficients[["effective::1999:top30"]]
), 
se = c(
  es_gs_stu_bin4$se[["effective::1992:top30"]],
  es_gs_stu_bin4$se[["effective::1993:top30"]],
  es_gs_stu_bin4$se[["effective::1994:top30"]],
  0,
  es_gs_stu_bin4$se[["effective::1996:top30"]],
  es_gs_stu_bin4$se[["effective::1997:top30"]],
  es_gs_stu_bin4$se[["effective::1998:top30"]],
  es_gs_stu_bin4$se[["effective::1999:top30"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999'),
    ordered = TRUE
  ),
group = "Top 30 Percentile", position = 4)

rb <- data.table::rbindlist(list(rb1_df, rb2_df, rb3_df, rb4_df))

rb$group <- rb$group %>%
  factor(levels = c(">= Mean", ">= Median", "Top 40 Percentile", "Top 30 Percentile"))
write.csv(rb, paste(dict_wd, "es_gs_stu_bin.csv", sep = ""))

ggplot(data = rb, aes(y = coef, x = yr, group = group, color = group)) +
  geom_point(aes(color = group), size = 1.5, position = position_nudge(x = (rb$position - 2.5)*0.12 )) +
  geom_errorbar(aes(ymin = (coef- 1.96 * se), ymax = (coef + 1.96 * se), color = group), width =0.1, 
                position = position_nudge(x = (rb$position - 2.5)*0.12))+
  scale_color_manual(values = c("#144062", "#1e6091", "#168aad", "#76c893")) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = as.numeric(rb$yr)[5] - 0.5, linetype = 2)+
  scale_y_continuous(limits = c(-.04, .04), breaks = seq(-.04, .04, by = 0.01)) +
  xlab("Effective Year of Exposure") +
  ylab("Coefficients") +
  theme_classic()+
  lg_theme

ggsave("es_gs_stu_bin.png", path = "r_outputs/", width = 6.472, height = 4.5 )

# 4.3 HS Diploma/GED
es_hs_ged_bin1 <- feols(hs_ged ~ i(effective, mean, ref = 1995) 
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)

es_hs_ged_bin2 <- feols(hs_ged ~ i(effective, median, ref = 1995)
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)

es_hs_ged_bin3 <- feols(hs_ged ~ i(effective, top40, ref = 1995)
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)

es_hs_ged_bin4 <- feols(hs_ged ~ i(effective, top30, ref = 1995)
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)

rb1_df <- data.frame(coef =c(
  es_hs_ged_bin1$coefficients[["effective::1992:mean"]],
  es_hs_ged_bin1$coefficients[["effective::1993:mean"]],
  es_hs_ged_bin1$coefficients[["effective::1994:mean"]],
  0,
  es_hs_ged_bin1$coefficients[["effective::1996:mean"]],
  es_hs_ged_bin1$coefficients[["effective::1997:mean"]],
  es_hs_ged_bin1$coefficients[["effective::1998:mean"]],
  es_hs_ged_bin1$coefficients[["effective::1999:mean"]],
  es_hs_ged_bin1$coefficients[["effective::2000:mean"]],
  es_hs_ged_bin1$coefficients[["effective::2001:mean"]],
  es_hs_ged_bin1$coefficients[["effective::2002:mean"]]
), 
se = c(
  es_hs_ged_bin1$se[["effective::1992:mean"]],
  es_hs_ged_bin1$se[["effective::1993:mean"]],
  es_hs_ged_bin1$se[["effective::1994:mean"]],
  0,
  es_hs_ged_bin1$se[["effective::1996:mean"]],
  es_hs_ged_bin1$se[["effective::1997:mean"]],
  es_hs_ged_bin1$se[["effective::1998:mean"]],
  es_hs_ged_bin1$se[["effective::1999:mean"]],
  es_hs_ged_bin1$se[["effective::2000:mean"]],
  es_hs_ged_bin1$se[["effective::2001:mean"]],
  es_hs_ged_bin1$se[["effective::2002:mean"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002'),
    ordered = TRUE
  ),
group = ">= Mean", position = 1)

rb2_df <- data.frame(coef =c(
  es_hs_ged_bin2$coefficients[["effective::1992:median"]],
  es_hs_ged_bin2$coefficients[["effective::1993:median"]],
  es_hs_ged_bin2$coefficients[["effective::1994:median"]],
  0,
  es_hs_ged_bin2$coefficients[["effective::1996:median"]],
  es_hs_ged_bin2$coefficients[["effective::1997:median"]],
  es_hs_ged_bin2$coefficients[["effective::1998:median"]],
  es_hs_ged_bin2$coefficients[["effective::1999:median"]],
  es_hs_ged_bin2$coefficients[["effective::2000:median"]],
  es_hs_ged_bin2$coefficients[["effective::2001:median"]],
  es_hs_ged_bin2$coefficients[["effective::2002:median"]]
), 
se = c(
  es_hs_ged_bin2$se[["effective::1992:median"]],
  es_hs_ged_bin2$se[["effective::1993:median"]],
  es_hs_ged_bin2$se[["effective::1994:median"]],
  0,
  es_hs_ged_bin2$se[["effective::1996:median"]],
  es_hs_ged_bin2$se[["effective::1997:median"]],
  es_hs_ged_bin2$se[["effective::1998:median"]],
  es_hs_ged_bin2$se[["effective::1999:median"]],
  es_hs_ged_bin2$se[["effective::2000:median"]],
  es_hs_ged_bin2$se[["effective::2001:median"]],
  es_hs_ged_bin2$se[["effective::2002:median"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002'),
    ordered = TRUE
  ),
group = ">= Median", position = 2)

rb3_df <- data.frame(coef =c(
  es_hs_ged_bin3$coefficients[["effective::1992:top40"]],
  es_hs_ged_bin3$coefficients[["effective::1993:top40"]],
  es_hs_ged_bin3$coefficients[["effective::1994:top40"]],
  0,
  es_hs_ged_bin3$coefficients[["effective::1996:top40"]],
  es_hs_ged_bin3$coefficients[["effective::1997:top40"]],
  es_hs_ged_bin3$coefficients[["effective::1998:top40"]],
  es_hs_ged_bin3$coefficients[["effective::1999:top40"]],
  es_hs_ged_bin3$coefficients[["effective::2000:top40"]],
  es_hs_ged_bin3$coefficients[["effective::2001:top40"]],
  es_hs_ged_bin3$coefficients[["effective::2002:top40"]]
), 
se = c(
  es_hs_ged_bin3$se[["effective::1992:top40"]],
  es_hs_ged_bin3$se[["effective::1993:top40"]],
  es_hs_ged_bin3$se[["effective::1994:top40"]],
  0,
  es_hs_ged_bin3$se[["effective::1996:top40"]],
  es_hs_ged_bin3$se[["effective::1997:top40"]],
  es_hs_ged_bin3$se[["effective::1998:top40"]],
  es_hs_ged_bin3$se[["effective::1999:top40"]],
  es_hs_ged_bin3$se[["effective::2000:top40"]],
  es_hs_ged_bin3$se[["effective::2001:top40"]],
  es_hs_ged_bin3$se[["effective::2002:top40"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002'),
    ordered = TRUE
  ),
group = "Top 40 Percentile", position = 3)

rb4_df <- data.frame(coef =c(
  es_hs_ged_bin4$coefficients[["effective::1992:top30"]],
  es_hs_ged_bin4$coefficients[["effective::1993:top30"]],
  es_hs_ged_bin4$coefficients[["effective::1994:top30"]],
  0,
  es_hs_ged_bin4$coefficients[["effective::1996:top30"]],
  es_hs_ged_bin4$coefficients[["effective::1997:top30"]],
  es_hs_ged_bin4$coefficients[["effective::1998:top30"]],
  es_hs_ged_bin4$coefficients[["effective::1999:top30"]],
  es_hs_ged_bin4$coefficients[["effective::2000:top30"]],
  es_hs_ged_bin4$coefficients[["effective::2001:top30"]],
  es_hs_ged_bin4$coefficients[["effective::2002:top30"]]
), 
se = c(
  es_hs_ged_bin4$se[["effective::1992:top30"]],
  es_hs_ged_bin4$se[["effective::1993:top30"]],
  es_hs_ged_bin4$se[["effective::1994:top30"]],
  0,
  es_hs_ged_bin4$se[["effective::1996:top30"]],
  es_hs_ged_bin4$se[["effective::1997:top30"]],
  es_hs_ged_bin4$se[["effective::1998:top30"]],
  es_hs_ged_bin4$se[["effective::1999:top30"]],
  es_hs_ged_bin4$se[["effective::2000:top30"]],
  es_hs_ged_bin4$se[["effective::2001:top30"]],
  es_hs_ged_bin4$se[["effective::2002:top30"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002'),
    ordered = TRUE
  ),
group = "Top 30 Percentile", position = 4)

rb <- data.table::rbindlist(list(rb1_df, rb2_df, rb3_df, rb4_df))

rb$group <- rb$group %>%
  factor(levels = c(">= Mean", ">= Median", "Top 40 Percentile", "Top 30 Percentile"))
write.csv(rb, paste(dict_wd, "es_hs_ged_bin.csv", sep = ""))

ggplot(data = rb, aes(y = coef, x = yr, group = group, color = group)) +
  geom_point(aes(color = group), size = 1.25, position = position_nudge(x = (rb$position - 2.5)*0.12 )) +
  geom_errorbar(aes(ymin = (coef- 1.96 * se), ymax = (coef + 1.96 * se), color = group), width =0.1, 
                position = position_nudge(x = (rb$position - 2.5)*0.12))+
  scale_color_manual(values = c("#144062", "#1e6091", "#168aad", "#76c893" )) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = as.numeric(rb$yr)[5] - 0.5, linetype = 2)+
  scale_y_continuous(limits = c(-.03, .03), breaks = seq(-.03, .03, by = 0.01)) +
  xlab("Effective Year of Exposure") +
  ylab("Coefficients") +
  theme_classic()+
  lg_theme

ggsave("es_hs_ged_bin.png", path = "r_outputs/", width = 6.472, height = 4.5)

# 4.4 Bachelor's Degree
es_bchl_bin1 <- feols(bchl ~ i(effective, mean, ref = 1995) 
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)

es_bchl_bin2 <- feols(bchl ~ i(effective, median, ref = 1995)
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)

es_bchl_bin3 <- feols(bchl ~ i(effective, top40, ref = 1995)
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)

es_bchl_bin4 <- feols(bchl ~ i(effective, top30, ref = 1995)
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)

rb1_df <- data.frame(coef =c(
  es_bchl_bin1$coefficients[["effective::1992:mean"]],
  es_bchl_bin1$coefficients[["effective::1993:mean"]],
  es_bchl_bin1$coefficients[["effective::1994:mean"]],
  0,
  es_bchl_bin1$coefficients[["effective::1996:mean"]],
  es_bchl_bin1$coefficients[["effective::1997:mean"]],
  es_bchl_bin1$coefficients[["effective::1998:mean"]],
  es_bchl_bin1$coefficients[["effective::1999:mean"]]
), 
se = c(
  es_bchl_bin1$se[["effective::1992:mean"]],
  es_bchl_bin1$se[["effective::1993:mean"]],
  es_bchl_bin1$se[["effective::1994:mean"]],
  0,
  es_bchl_bin1$se[["effective::1996:mean"]],
  es_bchl_bin1$se[["effective::1997:mean"]],
  es_bchl_bin1$se[["effective::1998:mean"]],
  es_bchl_bin1$se[["effective::1999:mean"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999'),
    ordered = TRUE
  ),
group = ">= Mean", position = 1)

rb2_df <- data.frame(coef =c(
  es_bchl_bin2$coefficients[["effective::1992:median"]],
  es_bchl_bin2$coefficients[["effective::1993:median"]],
  es_bchl_bin2$coefficients[["effective::1994:median"]],
  0,
  es_bchl_bin2$coefficients[["effective::1996:median"]],
  es_bchl_bin2$coefficients[["effective::1997:median"]],
  es_bchl_bin2$coefficients[["effective::1998:median"]],
  es_bchl_bin2$coefficients[["effective::1999:median"]]
), 
se = c(
  es_bchl_bin2$se[["effective::1992:median"]],
  es_bchl_bin2$se[["effective::1993:median"]],
  es_bchl_bin2$se[["effective::1994:median"]],
  0,
  es_bchl_bin2$se[["effective::1996:median"]],
  es_bchl_bin2$se[["effective::1997:median"]],
  es_bchl_bin2$se[["effective::1998:median"]],
  es_bchl_bin2$se[["effective::1999:median"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999'),
    ordered = TRUE
  ),
group = ">= Median", position = 2)

rb3_df <- data.frame(coef =c(
  es_bchl_bin3$coefficients[["effective::1992:top40"]],
  es_bchl_bin3$coefficients[["effective::1993:top40"]],
  es_bchl_bin3$coefficients[["effective::1994:top40"]],
  0,
  es_bchl_bin3$coefficients[["effective::1996:top40"]],
  es_bchl_bin3$coefficients[["effective::1997:top40"]],
  es_bchl_bin3$coefficients[["effective::1998:top40"]],
  es_bchl_bin3$coefficients[["effective::1999:top40"]]
), 
se = c(
  es_bchl_bin3$se[["effective::1992:top40"]],
  es_bchl_bin3$se[["effective::1993:top40"]],
  es_bchl_bin3$se[["effective::1994:top40"]],
  0,
  es_bchl_bin3$se[["effective::1996:top40"]],
  es_bchl_bin3$se[["effective::1997:top40"]],
  es_bchl_bin3$se[["effective::1998:top40"]],
  es_bchl_bin3$se[["effective::1999:top40"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999'),
    ordered = TRUE
  ),
group = "Top 40 Percentile", position = 3)

rb4_df <- data.frame(coef =c(
  es_bchl_bin4$coefficients[["effective::1992:top30"]],
  es_bchl_bin4$coefficients[["effective::1993:top30"]],
  es_bchl_bin4$coefficients[["effective::1994:top30"]],
  0,
  es_bchl_bin4$coefficients[["effective::1996:top30"]],
  es_bchl_bin4$coefficients[["effective::1997:top30"]],
  es_bchl_bin4$coefficients[["effective::1998:top30"]],
  es_bchl_bin4$coefficients[["effective::1999:top30"]]
), 
se = c(
  es_bchl_bin4$se[["effective::1992:top30"]],
  es_bchl_bin4$se[["effective::1993:top30"]],
  es_bchl_bin4$se[["effective::1994:top30"]],
  0,
  es_bchl_bin4$se[["effective::1996:top30"]],
  es_bchl_bin4$se[["effective::1997:top30"]],
  es_bchl_bin4$se[["effective::1998:top30"]],
  es_bchl_bin4$se[["effective::1999:top30"]]
), 
yr =c('1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999') %>%
  factor(
    levels = c( '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999'),
    ordered = TRUE
  ),
group = "Top 30 Percentile", position = 4)

rb <- data.table::rbindlist(list(rb1_df, rb2_df, rb3_df, rb4_df))

rb$group <- rb$group %>%
  factor(levels = c(">= Mean", ">= Median", "Top 40 Percentile", "Top 30 Percentile"))
write.csv(rb, paste(dict_wd, "es_bchl_bin.csv", sep = ""))

ggplot(data = rb, aes(y = coef, x = yr, group = group, color = group)) +
  geom_point(aes(color = group), size = 1.5, position = position_nudge(x = (rb$position - 2.5)*0.12 )) +
  geom_errorbar(aes(ymin = (coef- 1.96 * se), ymax = (coef + 1.96 * se), color = group), width =0.1, 
                position = position_nudge(x = (rb$position - 2.5)*0.12))+
  scale_color_manual(values = c("#144062", "#1e6091", "#168aad", "#76c893")) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = as.numeric(rb$yr)[5] - 0.5, linetype = 2)+
  scale_y_continuous(limits = c(-.06, .06), breaks = seq(-.06, .06, by = 0.02)) +
  xlab("Effective Year of Exposure") +
  ylab("Coefficients") +
  theme_classic()+
  lg_theme

ggsave("es_bchl_bin.png", path = "r_outputs/", width = 6.472, height = 4.5 )

################################################################################
#### 5. Educational outcomes: alternative exposure measurement
################################################################################
# 5.1 Jan.-Jun., 1989-1993
dd_clg_stu_exp1 <- feols(clg_stu ~ cnsar_fh:post
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)
summary(dd_clg_stu_exp1)

dd_gs_stu_exp1 <- feols(gs_stu ~ cnsar_fh:post
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)
summary(dd_gs_stu_exp1)

dd_hs_ged_exp1 <- feols(hs_ged ~ cnsar_fh:post
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)
summary(dd_hs_ged_exp1)

dd_bchl_exp1 <- feols(bchl ~ cnsar_fh:post
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)
summary(dd_bchl_exp1)

dd_clg_stu_exp2 <- feols(clg_stu ~ median_fh:post
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)
summary(dd_clg_stu_exp2)

dd_gs_stu_exp2 <- feols(gs_stu ~ median_fh:post
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)
summary(dd_gs_stu_exp2)

dd_hs_ged_exp2 <- feols(hs_ged ~ median_fh:post
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)
summary(dd_hs_ged_exp2)

dd_bchl_exp2 <- feols(bchl ~ median_fh:post
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)
summary(dd_bchl_exp2)

etable(dd_clg_stu_exp1, dd_gs_stu_exp1, dd_hs_ged_exp1, dd_bchl_exp1,
       keep = c("%cnsar_fh:post"),
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~ r2 + n + wt_mean + p_eff_fh + coef_fh_rs + se_fh_rs, tex = TRUE)

etable(dd_clg_stu_exp2, dd_gs_stu_exp2, dd_hs_ged_exp2, dd_bchl_exp2,
       keep = c("%median_fh:post"),
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~ r2 + n + wt_mean + p_eff_mdfh , tex = TRUE)


# 5.2 Jan.-Jun., 93
dd_clg_stu_exp3 <- feols(clg_stu ~ cnsar93_res:post
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)
summary(dd_clg_stu_exp3)

dd_gs_stu_exp3 <- feols(gs_stu ~ cnsar93_res:post
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)
summary(dd_gs_stu_exp3)

dd_hs_ged_exp3 <- feols(hs_ged ~ cnsar93_res:post
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)
summary(dd_hs_ged_exp3)

dd_bchl_exp3 <- feols(bchl ~ cnsar93_res:post
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)
summary(dd_bchl_exp3)

dd_clg_stu_exp4 <- feols(clg_stu ~ median93_res:post
                         + factor(statefip):factor(year) 
                         + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                         + factor(race) + hisp  + male + ssiv
                         | fips + qy,
                         cluster = ~fips, 
                         weights = ~ perwt,
                         data = early_dta)
summary(dd_clg_stu_exp4)

dd_gs_stu_exp4 <- feols(gs_stu ~ median93_res:post
                    + factor(statefip):factor(year) 
                    + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                    + factor(race) + hisp  + male + ssiv
                    | fips + qy,
                    cluster = ~fips, 
                    weights = ~ perwt,
                    data = late_dta)
summary(dd_gs_stu_exp4)

dd_hs_ged_exp4 <- feols(hs_ged ~ median93_res:post
                        + factor(statefip):factor(year) 
                        + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                        + factor(race) + hisp  + male + ssiv
                        | fips + qy,
                        cluster = ~fips, 
                        weights = ~ perwt,
                        data = acs_mg)
summary(dd_hs_ged_exp4)

dd_bchl_exp4 <- feols(bchl ~ median93_res:post
                      + factor(statefip):factor(year) 
                      + threshpreg + mhlaws_exp + waiver_exp + tanf_exp 
                      + factor(race) + hisp  + male + ssiv
                      | fips + qy,
                      cluster = ~fips, 
                      weights = ~ perwt,
                      data = late_dta)
summary(dd_bchl_exp4)

etable(dd_clg_stu_exp3, dd_gs_stu_exp3, dd_hs_ged_exp3, dd_bchl_exp3,
       keep = c("%cnsar93_res:post"),
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~ r2 + n + wt_mean + p_eff_93 + coef93_rs + se93_rs, tex = TRUE)

etable(dd_clg_stu_exp4, dd_gs_stu_exp4, dd_hs_ged_exp4, dd_bchl_exp4,
       keep = c("%median93_res:post"),
       style.tex = style.tex("aer"), 
       signif.code = c("***"=0.01, "**"=0.05, "*"=0.10, "+"=0.15 ),
       fitstat = ~ r2 + n + wt_mean + p_eff_md93 , tex = TRUE)