# utilities.R - Extra functions
# WKREBUILD_toolset/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# PARALLEL setup via doFuture

if(exists("cores")) {
  plan(multisession, workers=cores)
  options(doFuture.rng.onMisuse="ignore")
}

# icesmetrics {{{

# NAME = function ~ refpt, e.g. FMSY = fbar(om) / refpts(om)$Fmsy

icesmetrics <- list(FMSY=fbar~Fmsy, SBMSY=ssb~Btrigger,
  SBPA=ssb~Bpa, SBlim=ssb~Blim)

# }}}

# WKREBUILD2 performance statistics {{{

annualstats <- list(

  # P(SB>SBlim)
  PBlim=list(~iterMeans((SB/Blim) > 1), name="P(SB>SB[lim])",
    desc="Probability that spawner biomass is above Blim"),

  # P(SB>SBtrigger)
  PBtrigger=list(~iterMeans((SB/Btrigger) > 1), name="P(SB>B[trigger])",
    desc="Probability that spawner biomass is above Btrigger"),

  # mean(C)
  Cy=list(~iterMeans(C), name="mean(C)",
    desc="Mean catch per year"),

  # cv(C)
  cvCy=list(~sqrt(iterVars(C)) / iterMeans(C), name="cv(C)",
    desc="CV of catch per year")
)

fullstats <- list(

  # mean(C)
  C=list(~yearMeans(C), name="mean(C)",
    desc="Mean catch over years"),

  # AVVC
  # AVVC
  AAVC=list(~yearMeans(abs(C[, -1] - C[, -dim(C)[2]])/C[, -1]),
    name="AAV(C)", desc="Average annual variability in catch"),
  
  # IACC
  IACC=list(~100 * yearSums(abs(C[, -1] - C[, -dim(C)[2]]))/yearSums(C),
    name="IAC(C)",
    desc="Percentage inter-annual change in catch"),

  # P(SB < SBlim) at least once
  risk2=list(~yearMeans(iterMeans(((SB/Blim) < 1) > 0)),
    name="once(P(SB<B[limit]))",
    desc="ICES Risk 2, probability that spawner biomass is above Blim once")
)

# }}}

# firstyear {{{

# firstYear(iterMeans(SB/Blim > 1) >= 0.95)

firstYear <- function(x) {
  year <- as.numeric(dimnames(x)$year[match(TRUE, x)])
  return(FLQuant(year))
}
# }}}

# decisions {{{

decisions <- function(x, year=NULL, iter=NULL) {

  # EXTRACT tracking and args
  trac <- tracking(x)
  args <- args(x)

  # SET years if null
  if(is.null(year))
    year <- head(args$vy, -args$management_lag)

  # SET iters if not given
  if(is.null(iter))
    iter <- seq(dims(x)$iter)

  # FUNCTION to compute table along years
  .table <- function(d) {

    its <- dims(d)$iter
    dmns <- dimnames(d)

    if(its == 1) {
      data.frame(metric=dmns$metric, year=dmns$year, value=prettyNum(d))
    } else {
      data.frame(metric=dmns$metric, year=dmns$year,
        value=sprintf("%s (%s)", 
          prettyNum(apply(d, 1:5, median, na.rm=TRUE)),
          prettyNum(apply(d, 1:5, mad, na.rm=TRUE))))
    }
  }

  # COMPUTE tables
  res <- lapply(year, function(y) {
  
    # GET advice, data and management years
    ay  <-  an(y)
    dy <- ay - args$data_lag
    my  <- ay + args$management_lag

    # EXTRACT data year metrics
    dmet <- c("SB.om", "SB.obs", "SB.est", "met.hcr")
    dmet <- c("SB.om", "SB.obs", "SB.est")

    dout <- trac[dmet, ac(dy),,,, iter]

    # EXTRACT advice year metrics
    amet <- c("decision.hcr", "fbar.hcr", "hcr", "fbar.isys", "isys",
      "fwd", "C.om")

    aout <- trac[amet, ac(ay),,,, iter]

    # EXTRACT management year metrics
    mmet <- "SB.om"
   
    mout <- trac[mmet, ac(my),,,, iter]

    # COMPUTE management year metrics effect, my / ay
    eout <- trac[mmet, ac(my),,,,iter] / trac[mmet, ac(ay),,,,iter]
    
    dimnames(eout)$metric <- paste0("diff(", mmet, ")")

    # BIND into single table
    rbind(.table(dout), .table(aout), .table(mout), .table(eout))
  })

  if(length(res) > 1)
    res <- cbind(res[[1]], do.call(cbind,
      lapply(res[-1], function(i) i[, -1])))
  else
    res <- res[[1]]

  return(res)
}
# }}}

# tac.is {{{

#' TAC implementation system module
#'
#' Performs a short term forecast (STF) for the target fishing mortality to
#' obtain the corresponding catch.
#'
#' A `fwdControl` object obtained from the 'hcr' step is applied in the
#' management year (`ay + mlag`) or years (`seq(ay + mlag, ay + mlag + freq`).
#' An assumption is made on the mortality in the assessment year (`ay`), which
#' becomes the intermediate year in this projection. By default this is set
#' to Fbar = Fsq, that is, the same fishing mortality estimated in the 
#' last data year (`ay - data_lag`).
#'
#' The projection applies a constant recruitment, equal to the geometric mean
#' over an specified number of years. By default all years minus the last two
#' are included in the calculation. An specific set of years can be employed,
#' by specifying a character vector of year names, or two values can be given
#' for the number of years to be inlcuded, counting from the last, and how many
#' years to exclude at the end. For example, `c(30, 2)` will use the last 30
#' years but excluding the last two, usually worst estimated.
#'
#' @param stk The perceived FLStock.
#' @param ctrl The fwdControl output by the *hcr* step, target must be 'fbar'.
#' @param args The MSE run arguments.
#' @param recyrs Years to use for geometric mean recruitment if projection. Defaults to all years minus the last two.
#' @param Fdevs Deviances on the fbar input to incorporate error and bias when MP is run using the pseudo-estimators 'perfect.sa' or 'shortcut.sa'.
#' @param dtaclow Limit to decreases in output catch, as a proportional change (0.85 for 15%). Applied only when metric > lim, as set by 'hcr' step.
#' @param dtacupp Limit to increases in output catch, as a proportional change (1.15 for 15%). Applied only when metric > lim, as set by 'hcr' step.
#' @param fmin Minimum fbar to apply when catch change limits are use.
#' @param initac Initial catch from which to compute catch change limits. Defaults to previous observed catch.
#' @param tracking The tracking object.
#' @examples
#' data(sol274)
#' # Setup control with tac.is
#' control <- mpCtrl(list(est=mseCtrl(method=perfect.sa),
#'   hcr=mseCtrl(method=hockeystick.hcr,
#'     args=list(lim=0, trigger=4.3e5, target=0.21)),
#'   isys=mseCtrl(method=tac.is, args=list(recyrs=-3, output='landings'))))
#' # Run MP until 2025
#' run <- mp(om, oem, ctrl=control, args=list(iy=2021, fy=2024))
#' # Plot run time series
#' plot(om, TAC.IS=run)

tac.is <- function(stk, ctrl, args, output="catch", recyrs=-2,
  Fdevs=fbar(stk) %=% 1, dtaclow=NA, dtacupp=NA, fmin=0, reuse=TRUE,
  initac=metrics(stk, output)[, ac(iy - 1)], tracking) {

  # EXTRACT args
  spread(args)

  # SET control years
  cys <- seq(ay + management_lag, ay + management_lag + frq - 1)

  # PREPARE stk for cys, biology as in last nsqy years
  fut <- fwdWindow(stk, end=cys[length(cys)], nsq=nsqy)

  # PARSE recyrs if numeric
  id <- dimnames(stk)$year

  # COERCE to list
  if(!is.list(recyrs)) {
    recyrs <- list(recyrs)
  }
  
  # PARSE list
  for(i in recyrs) {
    if(is(i, 'character')) {
      id <- id[!id %in% i]
    } else if(all(i < 0)) {
      if(length(i) == 1)
        id <- rev(rev(id)[-seq(abs(i))])
      else
        id <- rev(rev(id)[i])
    } else if(all(i > 0)) {
      id <- rev(rev(id)[seq(abs(i))])
    }
  }

  # SET years to use
  recyrs <- id

  # CHECK recyrs
  if(!all(recyrs %in% dimnames(stk)$year)) {
    stop("'recyrs' cannot be found in input stk")
  }

  # TODO: OTHER rec options
  
  # SET GM recruitment from past

  gmnrec <- exp(yearMeans(log(rec(stk)[, recyrs])))

  srr <- predictModel(model=rec~a, params=FLPar(a=gmnrec))

  # STORE geomeanrec value
  track(tracking, "gmrec.isys", ay + management_lag) <- gmnrec

  # ADD F deviances for 1 year
  
  # reuse = TRUE
  if(isTRUE(reuse) | toupper(reuse) == 'F') {
    ftar <- rep(c(ctrl[1,]$value * Fdevs[, ac(cys[1])]), length(cys))
  # reuse = FALSE
  } else {
    ftar <- c(ctrl$value * Fdevs[, ac(cys)])
  }

  # TRACK Ftarget
  track(tracking, "fbar.isys", cys) <- ftar

  # FORECAST for iyrs and my IF mlag > 0,
  if(management_lag > 0) {
 
    # SET F for intermediate year
    fsq <- fbar(stk)[, ac(dy)]

    # TODO: ADD TAC option

    # CONSTRUCT fwd control
    fctrl <- fwdControl(
      # ay as intermediate with Fsq TODO: Other options
      list(year=seq(ay - data_lag + 1, length=management_lag), quant="fbar",
        value=rep(c(fsq), management_lag)),
      # target
                        # TODO: 2025 and 2026 same
      list(year=cys, quant="fbar", value=c(ftar))
    )

  # else only for my
  } else {
    fctrl <- fwdControl(
      list(year=ay + management_lag, quant="fbar", value=ftar))
  }

  # RUN STF ffwd
  fut <- ffwd(fut, sr=srr, control=fctrl)

  # ID iters where hcr set met trigger and F > fmin
  id <- c(tracking[[1]]["decision.hcr", ac(ay)] > 2) &
    c(fbar(fut)[, ac(ay + management_lag)] > fmin)

  # EXTRACT catches
  if(isTRUE(reuse) | toupper(reuse) == "C") {
    TAC <- expand(catch(fut)[, ac(cys)[1]], year=seq(length(cys)))
  } else {
    TAC <- catch(fut)[, ac(cys)]
  }

  # GET TAC dy / ay - 1
  if(ay == iy)
    prev_tac <- rep(c(initac), length=args$it)
  else
    prev_tac <- c(tracking[[1]]["isys", ac(ay)])

  # APPLY upper and lower TAC limit, if not NA and only for id iters
  if(!is.na(dtacupp)) {
    iter(TAC, id) <- pmin(c(iter(TAC, id)), prev_tac[id] * dtacupp)
  }
  if(!is.na(dtaclow)) {
    iter(TAC, id) <- pmax(c(iter(TAC, id)), prev_tac[id] * dtaclow)
  }

  # CONSTRUCT fwdControl
  ctrl <- fwdControl(lapply(seq(length(cys)), function(x)
    list(year=cys[x], quant=output, value=TAC[,x])))
    
  return(list(ctrl=ctrl, tracking=tracking))
}

# }}}
