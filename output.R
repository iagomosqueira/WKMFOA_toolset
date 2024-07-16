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

# COMPUTE yearly performance statistics, ADD F0 as reference

performance(runs) <- performance(runs, statistics=annualstats,
  years=list(short=2024:2034, all=2024:2042))


# --- TABLES

tables <- list()

# WHEN does stock recover (P(SB>Blim) >= 95%) by mp?

tables$recovery <- perf_year[statistic == "PBlim" & data > 0.95, .SD[1], 
  by=mp][order(year),]

perf[statistic=='firstyear' & year == 'all', data, by=.(mp)][order(data),]

# WHEN is P(B>Btrigger) > 50% by mp?

tables$status <- perf_year[statistic == "PBtrigger" & data > 0.50, .SD[1],
  by=mp]

# CREATE table of catch by mp and year (mp ~ year | C)

tables$catch_mp <- dcast(perf_year[statistic == 'C', .(data=mean(data)),
  by=.(year, mp, name, desc)], mp ~ year, value.var='data')

# CREATE table of all statistics by mp and year (statistic + mp ~ year)

tables$stats_mp <- dcast(perf_year[, .(data=mean(data)),
  by=.(year, mp, name, desc)], name + mp ~ year, value.var='data')


# --- TRACK decisions (EXAMPLES)

# TRACK decision for a single iter

decisions(advice, year=2024:2026, iter=1)
plot(iter(om, 1), iter(advice,1))

decisions(advice, year=2024:2026, iter=2)
plot(iter(om, 2), iter(advice, 2))

# TRACK decisions for multiple years and all iters
decisions(advice)

# SAVE
save(perf_year, perf, tables, file="output/output.rda", compress="xz")
