# model.R - Running rebuilding MPs
# WKREBUILD_toolset/model.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(icesTAF)
mkdir("model")

library(mse)

# CHOOSE number of cores for doFuture
plan(multisession, workers=4)

# LOAD oem and oem
load('data/data.rda', verbose = TRUE)

sr(om)

# - SET UP MP runs

# SET intermediate year, other args as default
mseargs <- list(iy=2024)

# SETUP standard ICES advice rule
arule <- mpCtrl(list(

  # (est)imation method: shortcut.sa + SSB deviances
  est = mseCtrl(method=shortcut.sa,
    args=list(SSBdevs=sdevs$SSB)),

  # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0, trigger=refpts(om)$Btrigger, target=refpts(om)$Fmsy,
    min=0, metric="ssb", output="fbar")),

  # (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances
  isys = mseCtrl(method=tac.is,
    args=list(recyrs=-2, fmin=0, Fdevs=sdevs$F))
  ))

# plot HCR
plot_hockeystick.hcr(arule$hcr, labels=c(lim="Blim", trigger="MSYBtrigger",
  min="", target="Ftarget")) +
  xlab(expression(hat(SSB))) + ylab(expression(bar(F)))+theme_bw()

# - RUN applying ICES advice rule
system.time(
  advice <- mp(om, ctrl=arule, args=mseargs)
)

# PLOT
plot(om, advice)


# --- RUN over alternative advice frequencies (1, 2, 3, 5)

library(future.apply)

# TODO: MOVE to mps()
system.time(
runs <- FLmses(future_lapply(setNames(nm=c(1, 2, 3, 5)), function(x)
  mp(om, ctrl=arule, args=list(iy=2024, fy=2050, frq=x), parallel=FALSE)))
)

plot(window(om, start=2000), runs)

# COMPUTE performance statistics

performance(runs) <- rbind(
  # annual
  performance(runs, statistics=annualstats, years=2024:2042),
  # by period
  performance(runs, statistics=fullstats, years=list(all=2024:2042)))

# SAVE
save(runs, file="model/model.rda", compress="xz")

# CLOSE cluster
plan(sequential)
