## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(RSTr)

## -----------------------------------------------------------------------------
head(maexample)

## -----------------------------------------------------------------------------
ma_mort <- maexample[which(!is.na(maexample$Year)), ]

## -----------------------------------------------------------------------------
head(ma_mort)

## -----------------------------------------------------------------------------
ma_data <- long_to_list_matrix(ma_mort, Deaths, Population, County.Code, Sex.Code, Year.Code)

## -----------------------------------------------------------------------------
Y <- xtabs(Deaths ~ County.Code + Sex.Code + Year.Code, data = ma_mort)
n <- xtabs(Population ~ County.Code + Sex.Code + Year.Code, data = ma_mort)
ma_data <- list(Y = Y, n = n)

## -----------------------------------------------------------------------------
ma_mort_mcar <- ma_mort[ma_mort$Year == 1979, ] # filter dataset to only show 1979 data
ma_data_mcar <- long_to_list_matrix(ma_mort_mcar, Deaths, Population, County.Code, Sex.Code)

## -----------------------------------------------------------------------------
ma_mort_car <- ma_mort[ma_mort$Year == 1979 & ma_mort$Sex == "Male", ] # filter dataset to only show 1979 data for men
ma_data_car <- long_to_list_matrix(ma_mort_car, Deaths, Population, County.Code)

