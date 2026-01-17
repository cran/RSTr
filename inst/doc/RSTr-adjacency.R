## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ggplot2)
library(RSTr)
is_cran <- identical(Sys.getenv("NOT_CRAN"), "false")
if (Sys.getenv("NOT_CRAN") == "") is_cran <- TRUE

## ----results = "hide"---------------------------------------------------------
head(mamap)

## -----------------------------------------------------------------------------
ma_shp <- sf::st_as_sf(mamap[order(mamap$GEOID), ])
ma_adj <- spdep::poly2nb(ma_shp)

## -----------------------------------------------------------------------------
ma_adj

## -----------------------------------------------------------------------------
no_neigh <- spdep::card(ma_adj) == 0

## -----------------------------------------------------------------------------
ma_shp[no_neigh, ]

## ----warning = FALSE, eval = !is_cran-----------------------------------------
# ggplot(mamap) +
#   geom_sf(aes(fill = NAME)) +
#   geom_sf_label(aes(label = NAME))

## ----eval = is_cran, echo = FALSE---------------------------------------------
# This map is shown on the pkgdown website, but skipped for CRAN checks
NULL

## -----------------------------------------------------------------------------
county_key <- seq_along(ma_shp$NAME)
names(county_key) <- ma_shp$NAME
county_key

## -----------------------------------------------------------------------------
ma_adj <- add_neighbors(ma_adj, c(1, 4, 10))

