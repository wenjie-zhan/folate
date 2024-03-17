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

# Directories
setwd("C:/Files/davis_research/folate_data/")
repo <- "C:/Users/wenji/OneDrive/github_files/GitHub/folate/"

################################################################################
#### Get PSID 2017 and 2019
################################################################################
# Load the data
psid_tmp <- read_xlsx("raw_data/psid/J318308.xlsx") %>%
  select(
    int_id = ER30001,
    per_id = ER30002,
    sex = ER66018,
    br_state = ER70874,
    hsp = ER70881,
    race1 = ER70882,
    race2 = ER70883,
    race3 = ER70884,
    race4 = ER70885,
    fam_id17 = ER34501,
    seq_id17 = ER34502,
    br_mo = ER34505,
    br_yr = ER34506,
    take_test = TA170785,
    sat_reading = TA170787,
    sat_math = TA170788,
    act_composite = TA170789,
    sex = 
  )



