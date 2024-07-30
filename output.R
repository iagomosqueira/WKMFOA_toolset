# output.R - Performance evaluation and output tables
# WKREBUILD_toolset/output.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(icesTAF)
mkdir("output")

library(mse)

source("utilities.R")

# LOAD model.R outputs

load("model/model.rda")


# --- TABLES

tables <- list()
perf <- performance(runs)
unique(perf$statistic)

# IS P(SB < Blim) > 0.95 in any year?
perf[statistic == 'PBlim',]
perf[statistic == 'PBlim', .(PBlim=max(data)), by=.(mp)]

# WHAT are the yearly average catch levels?

tables$catch_mp <- dcast(perf[statistic == 'Cy', .(data=mean(data)),
  by=.(year, mp, name, desc)], mp ~ year, value.var='data')

# WHAT are the long-term catch levels?

# WHAT is the catch variability?

  dcast(perf[statistic == 'AAVC', .(data=mean(data)),
  by=.(year, mp, name, desc)], mp ~ year, value.var='data')

# WHAT is the catch variability?

  dcast(perf[statistic == 'IACC', .(data=mean(data)),
  by=.(year, mp, name, desc)], mp ~ year, value.var='data')


# CREATE table of all statistics by mp and year (statistic + mp ~ year)

tables$stats_mp <- dcast(perf[, .(data=mean(data)),
  by=.(year, mp, name, desc)], name + mp ~ year, value.var='data')

# SAVE
save(tables, file="output/output.rda", compress="xz")





# --- TRACK decisions (EXAMPLES)

# TRACK decision for a single iter

decisions(runs[[1]], year=2024:2026, iter=1)
plot(iter(om, 1), iter(advice,1))

# TRACK decisions for multiple years and all iters
decisions(runs[[1]])

