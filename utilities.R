# utilities.R - Extra functions
# WKREBUILD_toolset/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# SET doFuture to ignore warning on RNG
options(doFuture.rng.onMisuse="ignore")

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

# fwd function with process error added {{{

pefwd <- function(object, catch=NULL, fbar = NULL, deviances = NULL, sr = NULL  ){

  if(!class(om)== "FLom") {print("make sure that class of object is om"); return()}
  if(is.null(deviances)){print("For process error calculation, deviances should be provided"); return()}
  if(is.null(catch)&is.null(fbar)){print("'catch' or 'fbar' missing"); return()}
  if(is.null(sr)){print("missing 'sr'"); return()} 

  stock <- object@stock
  om <- object

  nages <- dims(stock)$age; print(paste0("nages = ", nages))
  # DATA year
  dy <- dims(stock)$maxyear; print(paste0("dy = ", dy))
  #HINDCAST YEAR
  hiny <- dims(catch)$minyear; print(paste0("hiny = ", hiny))

  # DETERMINISTIC SEQUENTIAL hindcast with process error (perr) added ----
      # COMPUTE process error, e = y+t/(y exp(-z_in_y)) across all ages except recruitment.
      expected_n_at_age <- (stock.n(stock)[-nages, ac(hiny:dy-1)] * exp(-z(stock)[-nages, ac(hiny:dy-1)]))
      observed_n_at_age <- stock.n(stock)[-1, ac(hiny:dy)]

      perr <- observed_n_at_age / expected_n_at_age
      
      #correct plusgroup
      perr[ac(dims(stock)$max), ] <- stock.n(stock)[nages, ac(hiny:dy)]/(quantSums(stock.n(stock)[(nages-1):nages, ac(hiny:dy -1)] *
                                       exp(-z(stock)[(nages-1):nages, ac(hiny:dy -1)])))
      
      dhind_perr <- stock
      
      # The process error is the right one for the first year, but now needs to be 
      # recalculated as the numbers at ages in years are corrected, so the population in year
      # y included process error.
      perry <- perr
      #print("Process error recalculation ---------")

      for(Y in hiny:dy) {

        # FWD(catch[y])
        
        #print(paste0("Fbar in ", Y, "= ",(as.numeric(fbar(stock)[, ac(Y)])[1])))
        dhind_perr <- fwd(dhind_perr, sr=rec(stock), f=fbar(stock)[, ac(Y)]) # no deviances in this case added by Dorleta, so the process error will not rely on deviations in recruitment

        stock.n(dhind_perr)[-1, ac(Y)] <- stock.n(dhind_perr)[-1, ac(Y)]*perry[, ac(Y)] 
        
        #print(stock.n(dhind_perr)[,ac(hiny:dy)] / stock.n(stock)[,ac(hiny:dy)]) #should be the same since no deviances added here
        
        # perr for next year
        if(Y < dy){
          # all ages minus recruitment
          expected_n_at_age_nextY <- (stock.n(dhind_perr)[-nages, ac(Y)] * exp(-z(dhind_perr)[-nages, ac(Y)]))
          observed_n_at_age_nextY <- stock.n(stock)[-1, ac(Y+1)]

          perry[,ac(Y+1)] <- observed_n_at_age_nextY/expected_n_at_age_nextY
          # correct plusgroup
          perry[ac(dims(stock)$max),ac(Y+1)] <- stock.n(stock)[nages, ac(Y+1)]/(quantSums(stock.n(dhind_perr)[(nages-1):nages, ac(Y)] *
                                                                        exp(-z(dhind_perr)[(nages-1):nages, ac(Y)])))
        }
         

      }

  # STOCHASTIC SEQUENTIAL hindcast with process error (perr) added ----

      # Stochastic recruitment + Process error with lognormal deviations
      hind_perr_log     <- stock

      for(YY in seq(hiny,dy)) {
      
        # FWD(catch[y])
        hind_perr_log     <- fwd(hind_perr_log,     sr=sr[, ac(hiny:dy)], f=fbar(stock)[, ac(YY)], deviances = deviances[, ac(YY)])
        hind_perr_log@stock.n[-1, ac(YY)]     <- stock.n(hind_perr_log)[-1, ac(YY)]*perry[, ac(YY)] 

        print(YY)

      }

      om@stock <- hind_perr_log

      return(om)

}

# }}}
