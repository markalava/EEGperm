################################################################################
###
### DATE:               2007
###
### AUTHOR:             Mark C. Wheldon (biostatmark@gmail.com)
###
### DESC:               F tests used in Section 3.2 of the paper:
###                     Wheldon M. C., Anderson M. J. & Johnson
###                     B. W. "Identifying Treatment Effects In
###                     Multi-Channel Measurements In
###                     Electroencephalographic Studies: Multivariate
###                     Permutation Tests And Multiple Comparisons."
###                     Aust. N.Z. J. Stat 49(4), 2007, 397--413.
###
### LICENCE:            Copyright (c) 2020 Mark C Wheldon,
###                     Released under the MIT Licence. See
### https://github.com/markalava/multi-channel-eeg-perm-tests/blob/master/LICENSE
###
################################################################################

### * SYNOPSIS
################################################################################

## 1. FUNCTIONS
## ___ F Ratios Functions: Calculates resampled F statistics
## ___ Quantile Calculation: Calculates quantiles of the resampling
## distribution.
## ___ Make Plots: Produces plots like Figure 2 in the paper.
##
## 2. EXAMPLE
## ___ Re-creates Figure 2 using a dummy data set with the same shape as that
## used in Section 3.2 (so it won't look exactly the same!).
##
## WARNING: The code is NOT optimized for speed! It takes a LONG time to run ...


### * 1. FUNCTIONS
################################################################################

### ** F Ratios Functions
################################################################################

##-------------------------------------------------##
#... Sums of Squares Function (needs be within FSS
#    function because uses variables defined there)

#... OUTPUT is a (LenTLevs X 3) matrix with F stats
#    for factor A, B and interaction as columns, one
#    row for each time point

#.......................BEGIN.......................#

SumSqs = function(resp, FA, FB,
        Q, sample.size, LenALevs, LenBLevs, PersonT,
        dfT, dfA, dfB, dfI, dfRep, dfRes) {

  resp.df = as.data.frame(resp)

    C = (colSums(resp)^2)/Q
    SS = colSums(resp^2)
    SST = SS - C

    AList = split(resp.df, FA)
    csumA = sapply(AList, FUN = function(el.byT) {
      colSums(el.byT)^2 } )
    rsumA = rowSums(csumA)
    SSA = rsumA/(sample.size*LenBLevs) - C
    MSA = SSA/dfA

    BList = split(resp.df, FB)
    csumB = sapply(BList, FUN = function(z) {
      colSums(z)^2 } )
    rsumB = rowSums(csumB)
    SSB = rsumB/(sample.size*LenALevs) - C
    MSB = SSB/dfB

    TreatGrpList = split(resp.df,
      list(FA, FB))
    csumtg = sapply(TreatGrpList, FUN = function(z) {
      colSums(z)^2} )
    rsumtg = rowSums(csumtg)
    SSI = rsumtg/sample.size - SSA - SSB - C
    MSI = SSI/dfI

    NameList = split(resp.df, PersonT)
    csumN = sapply(NameList, FUN = function(z) {
      colSums(z)^2 } )
    rsumN = rowSums(csumN)
    SSRep = rsumN/(LenALevs*LenBLevs) - C
    MSRep = SSRep/dfRep

    SSRes = SST - SSA - SSB - SSI - SSRep
    MSRes = SSRes/dfRes

    FA = MSA/MSRes
    FB = MSB/MSRes
    FI = MSI/MSRes

  Fmat = cbind(FA, FB, FI)

  return(Fmat)
  }

#........................END........................#
##-------------------------------------------------##



##-------------------------------------------------##
#... Fstat resample function forINTERACTIONS - calls
#    the SumSqs function

#... OUTPUT is a LIST with elements "obs.FI", "res.FI",
#    "TimeLevs" giving observed and resampled FI values
#    and the time points vector.
#
#    "obs.FI": a named vector with interaction FI stats
#    for each time point
#    "res.FI": a LenTLevs X res matrix with resampled
#    FI stats

#...................BEGIN FSSIntv2..................#

FSSInt = function(eeg.data, e.cols=129, res=5) {

  #... Variables needed

  Time = eeg.data$Time
  TimeLevs = unique(eeg.data$Time)
  LenTLevs = length(TimeLevs)

  Person = factor(eeg.data$Name)
  PLevs = levels(Person)
  sample.size = length(PLevs)

  FactorA = factor(eeg.data$Treatment) # make generic
  ALevs = levels(FactorA)
  LenALevs = length(ALevs)

  FactorB = factor(eeg.data$Listening) # make generic
  BLevs = levels(FactorB)
  LenBLevs = length(BLevs)

  Q = sample.size*LenALevs*LenBLevs

  dfT = nrow(eeg.data)/LenTLevs - 1
  dfA = LenALevs - 1
  dfB = LenBLevs - 1
  dfI = dfA*dfB
  dfRep = sample.size - 1
  dfRes = dfT - dfA - dfB - dfI - dfRep


  #... Create within-Time factor vectors

  edT = eeg.data[Time==TimeLevs[1],]
  PersonT = factor(edT$Name)
  FactorAT = factor(edT$Treatment) # TODO: make generic
  FactorBT = factor(edT$Listening)

  # Will this make things faster ??
  FacAT.n = as.numeric(FactorAT)
  FacBT.n = as.numeric(FactorBT)

  repA = rep(1:LenALevs, LenBLevs)
  FacA.mat = matrix(repA, nrow=LenALevs*LenBLevs,
    ncol=sample.size)
  repB = rep(1:LenBLevs, rep(LenALevs, LenBLevs))
  FacB.mat = matrix(repB, nrow=LenALevs*LenBLevs,
    ncol=sample.size)


  #... Convert data frame to list with one component
  #    for each Time value

  eeg.list = list(0)
  for (i in 1:LenTLevs) {
    eeg.list[[i]] = eeg.data[Time==TimeLevs[i],]
  }


  #... Calculate OBSERVED sum of FI values

  cat(paste("Observed stats..."))

  obs.FI.vals = lapply(eeg.list, FUN = function(el.byT) {
    resp = as.matrix(el.byT[,1:e.cols])
    SumSqs(resp, FA=FactorAT, FB=FactorBT,
         Q, sample.size, LenALevs, LenBLevs, PersonT,
        dfT, dfA, dfB, dfI, dfRep, dfRes)
  })

  obs.sumFI.vals = lapply(obs.FI.vals, FUN = function(z) {
    sum(z[,3])
  })
  # could do this as part of first lapply() but this
  # allows checking of actual values for testing
  # purposes

  names(obs.sumFI.vals) = paste("T", TimeLevs)

  cat(paste(" done"))


  #... Fit reduced model

  resids.list = lapply(eeg.list, FUN = function(el.byT) {
    resp = as.matrix(el.byT[,1:e.cols])
    red.m = lm(resp ~ PersonT + FactorAT + FactorBT)
    residuals(red.m)
  })


  #... Set up matrix for RESAMPLED Interaction F values

  res.FI.mat = matrix(0, nrow = LenTLevs, ncol=res)


  #... BEGIN RESAMPLING ............................#

  cat(paste("\n", "Beginning resampling... ", "\n"))

  for (r in 1:res) {

    res.FacA.mat = apply(FacA.mat, 2, sample)
    res.FacA = as.vector(res.FacA.mat)

    res.FacB.mat = apply(FacB.mat, 2, sample)
    res.FacB = as.vector(res.FacB.mat)

    Fvals.list = lapply(resids.list,
      FUN = function(resids.byT) {
        Fmat = SumSqs(resids.byT, FA=res.FacA,
          FB=res.FacB, Q, sample.size,
          LenALevs, LenBLevs, PersonT,
          dfT, dfA, dfB, dfI, dfRep, dfRes)
        FsumI = sum(Fmat[,3]) # third col is FI
      })

    Fvals.vec = unlist(Fvals.list)
    # make resampled sumF values into a vector

    if (sum(Fvals.vec<0) != 0) {
      warning ("Negative F values calculated")
    }

    res.FI.mat[,r] = Fvals.vec

    cat(paste(r, " "))

  }

  #... END RESAMPLING ..............................#


  #... Put observed and re-sampled F.I values into a
  #    list for output

  res.FI.list = list(obs.FI = unlist(obs.sumFI.vals),
    res.FI = res.FI.mat, TimeLevs = TimeLevs)

  return(res.FI.list)

}

#..................END FSS.I...................#
##-------------------------------------------------##


##-------------------------------------------------##
#... Fstat resample function for MAIN EFFECT A - calls
#    the SumSqs function

#... OUTPUT is a List with one element for the
#    OBSERVED sumF mainA values, one for resampled
#    sumF mainA values. The resampled values are
#    in a LenTLevs X res matrix, with one row for
#    each time-point

#.......................BEGIN.......................#

FSSMainA = function(eeg.data, e.cols=129, res=5) {

  #... Variables needed

  Time = eeg.data$Time
  TimeLevs = unique(eeg.data$Time)
  LenTLevs = length(TimeLevs)

  Person = factor(eeg.data$Name)
  PLevs = levels(Person)
  sample.size = length(PLevs)

  FactorA = factor(eeg.data$Treatment) # make generic
  ALevs = levels(FactorA)
  LenALevs = length(ALevs)

  FactorB = factor(eeg.data$Listening) # make generic
  BLevs = levels(FactorB)
  LenBLevs = length(BLevs)

  Q = sample.size*LenALevs*LenBLevs

  dfT = nrow(eeg.data)/LenTLevs - 1
  dfA = LenALevs - 1
  dfB = LenBLevs - 1
  dfI = dfA*dfB
  dfRep = sample.size - 1
  dfRes = dfT - dfA - dfB - dfI - dfRep


  #... Create within-Time factor vectors

  edT = eeg.data[Time==TimeLevs[1],]
  PersonT = factor(edT$Name)
  FactorAT = factor(edT$Treatment) # TODO: make generic
  FactorBT = factor(edT$Listening)

  # Will this make things faster ??
  FacAT.n = as.numeric(FactorAT)
  FacBT.n = as.numeric(FactorBT)

  Afactors.mat = matrix(1:LenALevs, nrow=LenALevs,
        ncol=Q/LenALevs, byrow=F)
  Bfactors.hlf.mat = matrix(1:LenBLevs, nrow=LenBLevs,
        ncol=Q/(2*LenBLevs))


  #... Convert data frame to list with one component
  #    for each Time value

  eeg.list = list(0)
  for (i in 1:LenTLevs) {
    eeg.list[[i]] = eeg.data[Time==TimeLevs[i],]
  }


  #... Calculate OBSERVED sum of F.Main.A values

  cat(paste("Observed stats..."))

  obs.FA.vals = lapply(eeg.list, FUN = function(el.byT) {
    resp = as.matrix(el.byT[,1:e.cols])
    SumSqs(resp, FA=FactorAT, FB=FactorBT,
         Q, sample.size, LenALevs, LenBLevs, PersonT,
        dfT, dfA, dfB, dfI, dfRep, dfRes)
  })

  obs.sumFA.vals = lapply(obs.FA.vals, FUN = function(z) {
    sum(z[,1]) # first col is FA
  })

  cat(paste(" done"))


  #... Set up matrix for resampled Main Effect F.A values

  res.FA.mat = matrix(0, nrow=LenTLevs, ncol=res)
  rownames(res.FA.mat) = paste("T", TimeLevs)


  #... BEGIN RESAMPLING ............................#

  cat(paste("\n", "Beginning resampling... ", "\n"))

  for (r in 1:res) {

    res.Afactors.mat = apply(Afactors.mat, 2,
        FUN=function(z) sample(z))
    res.Afactors.vect = as.vector(res.Afactors.mat)

    Fvals.list = lapply(eeg.list, FUN = function(el.byT) {
      resp = as.matrix(el.byT[,1:e.cols])
      Fmat = SumSqs(resp, FA=factor(res.Afactors.vect),
        FB=factor(FacBT.n), Q, sample.size, LenALevs,
        LenBLevs, PersonT,
        dfT, dfA, dfB, dfI, dfRep, dfRes)
      FsumA = sum(Fmat[,1]) # first col is FA
    })

    Fvals.vec = unlist(Fvals.list)
    # make resampled sumF values into a vector

    if (sum(Fvals.vec<0) != 0) {
      warning ("Negative F values calculated")
    }

    res.FA.mat[,r] = Fvals.vec

    cat(paste(r, " "))

  }

  #... END RESAMPLING ..............................#


  #... Put observed and re-sampled F.A values into a
  #    list for output

  res.FA.list = list(obs.FA = unlist(obs.sumFA.vals),
    res.FA = res.FA.mat, TimeLevs = TimeLevs)

  return(res.FA.list)

}

#..................END FSS.Main.A...................#
##-------------------------------------------------##



##-------------------------------------------------##
#... Fstat resample function for MAIN EFFECT B - calls
#    the SumSqs function

#... OUTPUT is a List with one element for the
#    OBSERVED sumF mainA values, one for resampled
#    sumF mainA values. The resampled values are
#    in a LenTLevs X res matrix, with one row for
#    each time-point

#................BEGIN FSSMainB...................#

FSSMainB = function(eeg.data, e.cols=129, res=5) {

  #... Variables needed

  Time = eeg.data$Time
  TimeLevs = unique(eeg.data$Time)
  LenTLevs = length(TimeLevs)

  Person = factor(eeg.data$Name)
  PLevs = levels(Person)
  sample.size = length(PLevs)

  FactorA = factor(eeg.data$Treatment) # make generic
  ALevs = levels(FactorA)
  LenALevs = length(ALevs)

  FactorB = factor(eeg.data$Listening) # make generic
  BLevs = levels(FactorB)
  LenBLevs = length(BLevs)

  Q = sample.size*LenALevs*LenBLevs

  dfT = nrow(eeg.data)/LenTLevs - 1
  dfA = LenALevs - 1
  dfB = LenBLevs - 1
  dfI = dfA*dfB
  dfRep = sample.size - 1
  dfRes = dfT - dfA - dfB - dfI - dfRep


  #... Create within-Time factor vectors

  edT = eeg.data[Time==TimeLevs[1],]
  PersonT = factor(edT$Name)
  FactorAT = factor(edT$Treatment) # TODO: make generic
  FactorBT = factor(edT$Listening)

  # Will this make things faster ??
  FacAT.n = as.numeric(FactorAT)
  FacBT.n = as.numeric(FactorBT)

  Afactors.mat = matrix(1:LenALevs, nrow=LenALevs,
        ncol=Q/LenALevs, byrow=F)
  Bfactors.hlf.mat = matrix(1:LenBLevs, nrow=LenBLevs,
        ncol=Q/(2*LenBLevs))


  #... Convert data frame to list with one component
  #    for each Time value

  eeg.list = list(0)
  for (i in 1:LenTLevs) {
    eeg.list[[i]] = eeg.data[Time==TimeLevs[i],]
  }


  #... Calculate OBSERVED sum of F.Main.B values

  cat(paste("Observed stats..."))

  obs.FB.vals = lapply(eeg.list, FUN = function(el.byT) {
    resp = as.matrix(el.byT[,1:e.cols])
    SumSqs(resp, FA=FactorAT, FB=FactorBT,
         Q, sample.size, LenALevs, LenBLevs, PersonT,
        dfT, dfA, dfB, dfI, dfRep, dfRes)
  })

  obs.sumFB.vals = lapply(obs.FB.vals, FUN = function(z) {
    sum(z[,2]) # second col is FA
  })

  cat(paste(" done"))


  #... Set up matrix for resampled Main Effect F.B values

  res.FB.mat = matrix(0, nrow=LenTLevs, ncol=res)
  rownames(res.FB.mat) = paste("T", TimeLevs)


  #... BEGIN RESAMPLING ............................#

  cat(paste("\n", "Beginning resampling... ", "\n"))

  for (r in 1:res) {

    res.Bfactors.hlf = apply(Bfactors.hlf.mat, 2,
      FUN=function(z) sample(z))

    res.Bfactors.mat = rbind(as.vector(res.Bfactors.hlf),
      as.vector(res.Bfactors.hlf))

    res.Bfactors.vect = as.vector(res.Bfactors.mat)

     Fvals.list = lapply(eeg.list, FUN = function(el.byT) {
      resp = as.matrix(el.byT[,1:e.cols])
      Fmat = SumSqs(resp, FA=factor(FacAT.n),
        FB=factor(res.Bfactors.mat), Q, sample.size,
        LenALevs, LenBLevs, PersonT,
        dfT, dfA, dfB, dfI, dfRep, dfRes)
      FsumB = sum(Fmat[,2]) # second col is FB
    })

    Fvals.vec = unlist(Fvals.list)
    # make resampled sumF values into a vector

    if (sum(Fvals.vec<0) != 0) {
      warning ("Negative F values calculated")
    }

    res.FB.mat[,r] = Fvals.vec

    cat(paste(r, " "))

  }

  #... END RESAMPLING ..............................#


  #... Put observed and re-sampled F.B values into a
  #    list for output

  res.FB.list = list(obs.FB = unlist(obs.sumFB.vals),
    res.FB = res.FB.mat, TimeLevs = TimeLevs)

  return(res.FB.list)

}

#..................END FSS.Main.B...................#
##-------------------------------------------------##


### ** Quantile Calculation
################################################################################

##-------------------------------------------------##
## Quantile calculation for INTERACTION
##-------------------------------------------------##

##-------------------------------------------------##
#... INPUT Flist.I is a list with one element per
#    time-point. Within elements, a vector with the
#    resampled interaction F values.
#    FI.obs should be a vector of length LenTLevs,
#    each element the observed F Int for each time-
#    point.

#... OUTPUT

#.......................BEGIN.......................#

EEGFquant2.I = function(FIlist, FI.obs=NULL,
  alpha=0.05) {

  FI.res = FIlist$res.FI
  if (is.null(FI.obs)) FI.obs = FIlist$obs.FI

  LenTLevs = length(FIlist$TimeLevs)
  res = ncol(FI.res)

  pvals = seq(from=1/res, to=1, length=res)
  # this is for calculation of EW CIs
  # NOTE limit to precision - machine zero, etc


  ##... Pointwise CIs ...##

  PWp = rep(0, LenTLevs)
  for (t in 1:LenTLevs) {
    PWp[t] = sum(FI.res[t,] >= FI.obs[t]) / res
  }

  PWQuants = apply(FI.res, 1, FUN=function(z) {
    quantile(z, probs=1-alpha)
  }
    )


  ##... Exptwise CIs ...##

  Time.Ps = apply(FI.res, 1, FUN=function(z) {
    ranks = rank(z, ties.method="min")
    rank = (res-ranks)+1
    pvals[rank]
  }
    )

  # Time.Ps is a matrix with one column for each time-
  # point, and one row for each permutation

  MinPVect = apply(Time.Ps, 1, min)

  # MinPVect is a sample of size NRESAMPLE from the
  # distribution of minimum p-values

  MinP = quantile(MinPVect, probs=alpha)


  ##... EWp ...##

  # No times the MINIMUM p-value in each re-sample is
  # LESS than the the observed (point-wise) p-value
  # for EACH PWp

  EWp = rep(0, LenTLevs)
  for (i in 1:LenTLevs) {
    EWp[i] = min(sum(MinPVect < PWp[i]) / res, 1)
  }

  # not sure if this is valid under permutation...
  seEWp = sqrt(EWp*(1-EWp)/res)


  ##... EWQuants ...##

  EWQuants = apply(FI.res, 1, FUN=function(z) {
    quantile(z, probs=1-MinP)
  })


  ##... Bonferroni CIs ...##

  BOQuants = apply(FI.res, 1, FUN=function(z) {
    quantile(z, probs=1-(alpha/LenTLevs))
    #Sidak: quantile(z, probs=1- (1 - (1 - alpha)^(1/LenTLevs)))
  })



  ## OUTPUT ##

  ResQuants = list(PWp = PWp, PWQuants = PWQuants,
    EWp = EWp, seEWp = seEWp, EWQuants = EWQuants,
    BOQuants = BOQuants, BOp=alpha/LenTLevs,
    alpha = alpha, FI.obs = FI.obs,
    TimeLevs = FIlist$TimeLevs)

  return(ResQuants)

}



##-------------------------------------------------##
## Quantile calculation for MAIN EFFECTS
##-------------------------------------------------##

EEGFquant.M = function(FMlist, FM.res=NULL, FM.obs=NULL,
  alpha=0.05) {

  # List elements a little different from interaction
  if (is.null(FM.res)) FM.res = FMlist$res.FM
  if (is.null(FM.obs)) FM.obs = FMlist$obs.FM

  LenTLevs = length(FMlist$TimeLevs)
  res = ncol(FM.res)

  pvals = seq(from=1/res, to=1, length=res)
  # this is for calculation of EW CIs
  # NOTE limit to precision - machine zero, etc


  ##... Pointwise CIs ...##

  PWp = rep(0, LenTLevs)
  for (t in 1:LenTLevs) {
    PWp[t] = sum(FM.res[t,] >= FM.obs[t]) / res
  }

  PWQuants = apply(FM.res, 1, FUN=function(z) {
    quantile(z, probs=1-alpha)
  }
    )


  ##... Experimentwise CIs ...##

  Time.Ps = apply(FM.res, 1, FUN=function(z) {
    ranks = rank(z, ties.method="min")
    rank = (res-ranks) + 1
    pvals[rank]
  }
    )

  # Time.Ps is a matrix with one column for each time-
  # point, and one row for each permutation

  MinPVect = apply(Time.Ps, 1, min)

  MinP = quantile(MinPVect, probs=alpha)


  ##... EWp ...##

  # No times the MINIMUM p-value in each re-sample is
  # LESS than the the observed (point-wise) p-value
  # for EACH PWp

  EWp = rep(0, LenTLevs)
  for (i in 1:LenTLevs) {
    EWp[i] = min(sum(MinPVect < PWp[i]) / res, 1)
  }

  # not sure if this is valid under permutation...
  seEWp = sqrt(EWp*(1-EWp)/res)


  ##... EW Quants ...##

  EWQuants = apply(FM.res, 1, FUN=function(z) {
    quantile(z, probs=1-MinP)
  })


  ## Bonferroni CIs ##

  BOQuants = apply(FM.res, 1, FUN=function(z) {
    quantile(z, probs=1-(alpha/LenTLevs))
    #Sidak: quantile(z, probs=1- (1 - (1 - alpha)^(1/LenTLevs)))
  })


  ## OUTPUT ##

  ResQuants = list(PWp = PWp, PWQuants = PWQuants,
    EWp = EWp, seEWp = seEWp, EWQuants = EWQuants,
    BOQuants = BOQuants, BOp=alpha/LenTLevs,
    alpha = alpha, FM.obs = FM.obs,
    TimeLevs = FMlist$TimeLevs)

  return(ResQuants)

}


### ** Make Plots (as in Figure 2 of the paper)
################################################################################

##-------------------------------------------------##
## FUNCTION FOR INTERACTION EFFECTS
##-------------------------------------------------##

EEGFplot.I = function(fquant, data.name="", sub="",
  TimeLevs=NULL, alpha=NULL) {

  if (is.null(TimeLevs)) TimeLevs = fquant$TimeLevs
  if (is.null(alpha)) alpha = fquant$alpha

  obs.I = fquant$FI.obs

  BOQuants = fquant$BOQuants
  PWQuants = fquant$PWQuants
  EWQuants = fquant$EWQuants

  par(cex=0.8)

  ymax=max(BOQuants)+50
  ymin=min(obs.I)-50

  plot(x=TimeLevs, obs.I, type="l", lwd=2,
       ylim=c(ymin, ymax),
       xlab="time (ms)", ylab=paste("sum of F Int"))
  abline(v=0, lty=2, col="blue")

  lines(x=TimeLevs, PWQuants, col="grey", lty=2)

  lines(x=TimeLevs, EWQuants, col="red", lty=5)

  lines(x=TimeLevs, BOQuants, col="green3", lty=4)

  diffs = EWQuants-obs.I
  sig.diff = diffs

  sig.times = TimeLevs[sig.diff<=0]
  sig.obs = obs.I[sig.diff<=0]
  if (length(sig.times > 0)) {
    for (i in 1:length(sig.times)) {
      lines(x=c(sig.times[i], sig.times[i]),
            y=c(-10,sig.obs[i]))
    }
  }

  title(main=paste("Sum of F Int", " with ", alpha*100,
          "% Critical Values", "\n", data.name, sep=""),
        sub=sub)

  legend(10, ymax, ncol=4, cex=0.75,
         lwd=c(2, 1, 1, 1),
         lty=c(1, 2, 5, 4),
         col=c("black", "grey", "red", "green3"),
         legend=c(paste("sum of F Int"),
           "Point-wise CVs",
           "Expt-wise CVs", "Bonferroni CVs"
           ))
}




##-------------------------------------------------##
## FUNCTION FOR MAIN EFFECTS
##-------------------------------------------------##

EEGFplot.M = function(fquant, data.name="", sub="",
  TimeLevs=NULL, alpha=NULL) {

  if (is.null(TimeLevs)) TimeLevs = fquant$TimeLevs
  if (is.null(alpha)) alpha = fquant$alpha

  obs.main = fquant$FM.obs

  BOQuants = unlist(fquant$BOQuants)
  PWQuants = unlist(fquant$PWQuants)
  EWQuants = unlist(fquant$EWQuants)

  par(cex=0.8)

  ymax=max(BOQuants)+50
  ymin=min(obs.main)-50

  plot(x=TimeLevs, obs.main, type="l", lwd=2,
       ylim=c(ymin, ymax),
       xlab="time (ms)", ylab=paste("sum of F Main"))
  abline(v=0, lty=2, col="blue")

  lines(x=TimeLevs, PWQuants, col="grey", lty=2)

  lines(x=TimeLevs, EWQuants, col="red", lty=5)

  lines(x=TimeLevs, BOQuants, col="green3", lty=4)

  diffs = EWQuants-obs.main
  sig.diff = diffs

  sig.times = TimeLevs[sig.diff<=0]
  sig.obs = obs.main[sig.diff<=0]
  if (length(sig.times > 0)) {
    for (i in 1:length(sig.times)) {
      lines(x=c(sig.times[i], sig.times[i]),
            y=c(-10,sig.obs[i]))
    }
  }

  title(main=paste("Sum of F Main", " with ", alpha*100,
          "% Critical Values", "\n", data.name, sep=""),
        sub=sub)

  legend(10, ymax, ncol=4, cex=0.75,
         lwd=c(2, 1, 1, 1),
         lty=c(1, 2, 5, 4),
         col=c("black", "grey", "red", "green3"),
         legend=c(paste("sum of F"),
           "Point-wise CVs",
           "Expt-wise CVs", "Bonferroni CVs"
           ))
}


### * 2. EXAMPLE
################################################################################

### ** Dummy data
################################################################################

## This isn't supposed to represent anything in particular. I made this data set
## simply to demonstrate that the code runs without error.

## Labels for times, treatments, listening condition and subjects ('Name').
dummy.data.labels <-
    expand.grid(Time = seq(from = 0, to = 500, by = 4)
                ,Treatment = c("N", "T")
                ,Listening = c("A", "B", "C")
                ,Name = letters[1:11]
                )

## AR(1) dummy data
set.seed(21)
dummy.data.electrodes <-
    replicate(n = 129
              ,expr = arima.sim(model = list(ar=0.99), n = nrow(dummy.data.labels)) + 175
              )

## Add something to timepoints 400--450, for Treatment=="N".
index1 <-
    dummy.data.labels$Time > 399 & dummy.data.labels$Time < 451 &
      dummy.data.labels$Treatment == "N"
dummy.data.electrodes[index1,] <-
    dummy.data.electrodes[index1,] + runif(n = sum(index1), min = 0.5, max = 1.5)

index2 <-
    dummy.data.labels$Time > 399 & dummy.data.labels$Time < 451 &
      dummy.data.labels$Treatment == "N" & dummy.data.labels$Listening == "A"
dummy.data.electrodes[index2,] <-
    dummy.data.electrodes[index2,] + runif(n = sum(index2), min = 0.5, max = 1)

## Check
plot(dummy.data.electrodes[order(dummy.data.labels$Time), 4]
     ,type = "l")

## Combine to create complete dataset
dummy.data <- cbind(dummy.data.electrodes, dummy.data.labels)


### ** Do the tests and plot
################################################################################

set.seed(22)

## 'res' sets number of resamples. It is set to '5' because this takes
## only minutes. In the paper we did 15000. This takes a /long/ time.

res <- 5


### F-test of interaction
###

FInt <- FSSInt(dummy.data, e.cols=129, res=res)


### F-tests for main effects
###

FmainA <- FSSMainA(dummy.data, e.cols=129, res=res)
FmainB <- FSSMainB(dummy.data, e.cols=129, res=res)


### Plots
###

#... Interaction ...#

Iquant = EEGFquant2.I(FInt)

dev.new()
EEGFplot.I(Iquant)


#... Main Effects ...#


# Time Levs not put in FmainA, so do now

FmainA$TimeLevs = FInt$TimeLevs

Mquant.A = EEGFquant.M(FmainA, FM.res=FmainA$res.FA,
  FM.obs=FmainA$obs.FA)

dev.new()
EEGFplot.M(Mquant.A)


FmainB$TimeLevs = FInt$TimeLevs

Mquant.B = EEGFquant.M(FmainB, FM.res=FmainB$res.FB,
  FM.obs=FmainB$obs.FB)

dev.new()
EEGFplot.M(Mquant.B)

