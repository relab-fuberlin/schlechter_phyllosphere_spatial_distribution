#!/usr/bin/env Rscript


##  LIBRARIES USED FOR SPATIAL DISTRIBUTION PAPER

# library
# tidy data
library(here)
library(tidyverse)
library(reshape2)
library(broom)
library(magrittr)

# plots
library(patchwork)
library(gridExtra)
library(RColorBrewer)
library(ggh4x)
library(ggpubr)
library(ggtern)
library(ggbeeswarm)
library(wesanderson)
library(ggdist)

# stats
library(car)
library(emmeans)
library(lsr)
library(multcompView)
library(spatstat)
library(MESS)
library(MASS)
library(betareg)

# dependencies
source("code/theme_rs_spatial.R")
