################################################################################
###
### DATE:               2007
###
### AUTHOR:             Mark C. Wheldon (biostatmark@gmail.com)
###
### DESC:               t tests used in Section 3.1 of the paper:
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
## ___ Permutation Functions: Calculates resampled t, t|sum| statistics
## ___ Quantile Calculation: Calculates quantiles of the resampling
## distribution.
## ___ Make Plots: Produces plots like Figure 1 in the paper.
##
## 2. EXAMPLE
## ___ Re-creates Figure 1 using a dummy data set with the same shape as that
## used in Section 3.1 (so it won't look exactly the same!).
##
## WARNING: The code is NOT optimized for speed! It takes a LONG time to run ...


### * 1. FUNCTIONS
################################################################################

##-------------------------------------------------##
## Permutation Function
##-------------------------------------------------##

# Venables & Ripley's (2002) binary function:

    binary.v <- function(x, no.pairs) {
       if(missing(no.pairs)) {
           mx <- max(x)
           no.pairs <- if(mx > 0) 1 +
             floor(log(mx, base = 2)) else 1
       }
       ans <- 0:(no.pairs - 1)
       lx <- length(x)
       x <- rep(x, rep(no.pairs, lx))
       x <- (x %/% 2^ans) %% 2
       # x%%y is x mod y, %/% is integer division
       dim(x) <- c(no.pairs, lx)
       x
       # new x is a matrix with each original x element
       # in its own column, expressed in binary to
       # no.pairs number of digits.
    }


# V & R's perm.t.test function with modification to
# stop re-calcultion of binary matrix each time:
# (Venables, W.N. & Ripley, B.D. (2002). Modern Applied Statistics with S, 4th edn. New York: Springer.)


perm.t.test3 <- function(d, bm0.5) {
# ttest is function(x) mean(x)/sqrt(var(x)/length(x))
# d is vector of observed differences, bm0.5 is described
# below.
    no.pairs <- length(d) # number of pairs or differences
    n <- 2^no.pairs # number of permutations
    x <- d * bm0.5
      # bm0.5 is a matrix with (n = no. of permutations)
      # columns (no.pairs = no. of pairs or differences per
      # permutation) rows, each column containing
      # binary repsntn of numbers 1:2^no.pairs to length(d)
      # digits.
      # NB that since have only length(d) digits, and
      # first binary posn is 2^0, the number 2^no.pairs
      # (the last column of matrix) will be a column of
      # zeros (since cannot represent 2^no.pairs in
      # no.pairs-1
      # binary digits).

      # When the above matrix is multiplied d, it gives
      # all combinations of
      # positive and negative d values as columns
    mx <- matrix(1/no.pairs, 1, no.pairs) %*% x
      # means of all permutions of differences
    s <- matrix(1/(no.pairs - 1), 1, no.pairs)
      # 1/n-1 divisor where n is number of pairs
    vx <- s %*% (x - matrix(mx, no.pairs, n, byrow=T))^2
      # 1/(n-1) * ( x-mean(x) )^2 ie: variance
    as.vector(mx/sqrt(vx/no.pairs))
      # t-value: mean / std error
}


#.....Begin function to build permutation distribution

EEGttest2 = function(data, stat="abs.t",
  NVar=dim(data)[2]-2,
  obs.t="values", trace=T) {

    # data must be a data frame of DIFFERENCES
    # NVar is used to specify number of electrodes
    # obs.t controls how observed t values are returned
    #   if at all

  if (length(NVar) > dim(data)[2]) {
    stop ("Number of variables must be less than number of columns in data frame")
  }

    if (!is.element(el=stat,
                  set=c("abs.t", "t"))) {
    stop ("stat must be one of 'abs.t' or 't'")
  }

  if (!is.element(el=obs.t,
                  set=c("none", "index", "values"))) {
    stop ("obs.t must be one of 'none', 'index', 'values'")
  }

  # TO DO:
  # Set up additional checks to make sure columns
  # "Time" and "Name" are present

  DiffTime = as.numeric(as.vector(data$Time))
  TimeLevs = unique(DiffTime)
  LenTLevs = length(TimeLevs)

  DiffName = factor(data$Name)
  NameLevs = levels(DiffName)
  LenNLevs = length(NameLevs)

  nperms = 2^LenNLevs


  # Binary matrix:

  Name.bm1 = 2*(binary.v(1:nperms, LenNLevs) - 0.5)

  # Set up matrix to contain permutation sum-of-t-values

  tperm.vect = rep(0, nperms * LenTLevs)


  # Set up time index so the permutations
  # can be dumped in the correct rows

  time.vector = rep(TimeLevs, rep(nperms, LenTLevs))

  for (TLevIND in TimeLevs) {

    if (trace==T) {
      cat(paste("t", TLevIND, " ", "e1", sep=""))
    }

    s.sum = rep(0, nperms)

    for (var in 1:NVar) {

      # pull out electrodes one at a time
      e = data[DiffTime==TLevIND,var]

      # use the V&R staterm.t.test function, putting
      # results in correct rows and adding to prev
      # electrode's value
      s = perm.t.test3(e, bm0.5=Name.bm1)
      s[is.na(s)] = 0

      if (stat == "abs.t") s.sum = s.sum + abs(s)
      if (stat == "t") s.sum = s.sum + s

    }

    tperm.vect[time.vector==TLevIND] = s.sum

    if (trace==T) {
      cat(paste(" ... ", "e", NVar, "\n", sep=""))
    }
  }

  # Bind time column
  tperm.mat = cbind(tperm.vect, time=time.vector)

  # Make colnames a bit prettier
  colnames(tperm.mat) = c("Sum t", "time")


  if (obs.t=="none") {
    tperm.list = list(perms=tperm.vect, Times=DiffTime,
      nperms=nperms, stat=stat)
    return(tperm.list)
  }

  if (obs.t=="index") {

  # Pull out the observed t-stats: due to nature of
  #   the "binary matrix" this will be the penultimate
  #   value in each permuation set, or (-1)*the last
  #   value

    index.vect = rep(c(rep(F, nperms - 2), T, F),
      LenTLevs)
    tperm.list = list(perms=tperm.vect,
      obs.index=index.vect, Times=DiffTime,
      nperms=nperms, stat=stat)
    return(tperm.list)
  }

  if (obs.t=="values") {
    index.vect = rep(c(rep(F, nperms - 2), T, F),
      LenTLevs)
      sample.t = tperm.vect[index.vect]
    tperm.list = list(perms.sumt=tperm.vect,
      obs.sumt=sample.t, Times=DiffTime,
      nperms=nperms, stat=stat)
    return(tperm.list)
  }

}
#.....END


##---------------------///////---------------------##
##
## Take result of EEGttest2() and calculate PW, EW
## p values and confidence intervals
##
##---------------------///////---------------------##

##-------------------------------------------------##
## Quantile calculation
##-------------------------------------------------##

EEGtquant2 = function(PermT, alpha=0.05) {
  SumPerms = PermT$perms.sumt
  SumObs = PermT$obs.sumt

  DiffTime = PermT$Times
  TimeLevs = unique(DiffTime)
  LenTLevs = length(TimeLevs)

  nperms = PermT$nperms

  PermTLevs = rep(TimeLevs, rep(nperms, LenTLevs))
  # this indexes time-points in the SumPerms vector:
  # the EEGttest_2 program arranges time-points in this
  # way

  pvals = seq(from=1/nperms, to=1, length=nperms)
  # this is for calculation of EW CIs
  # NOTE limit to precision - machine zero, etc


  ##... Absolute t ...###

  if (PermT$stat == "abs.t") {

    ##... Pointwise CIs ...##

    cat(paste(" Pointwise CIs ..."))

    PWp = rep(0, LenTLevs)
    for (t in 1:LenTLevs) {
      PWp[t] = sum(SumPerms[PermTLevs==TimeLevs[t]] >=
           SumObs[TimeLevs==TimeLevs[t]]) / nperms
    }

    PWQuants = as.vector(
      tapply(SumPerms, INDEX=PermTLevs,
             FUN=function(z){
               quantile(z, probs=1-alpha) # one-tailed test
             }
             )
      )


    ##... Exptwise CIs ...##

    cat(paste(" done", "\n", "Experimentwise CIs ..."))



    Time.Ps = tapply(SumPerms, INDEX=PermTLevs, simplify=F,
      FUN=function(z) {
        ranks = rank(z, ties.method="min")
        rank = (nperms-ranks)+1
        # method "min" means that nperms - ranks + 1 gives
        # number of values greater than or equal to
        # each permutation sum of t within time-points
        pvals[rank]
      }
      )

    PermsMatrix = matrix(0, nrow=nperms, ncol=LenTLevs)

    for (i in 1:LenTLevs) PermsMatrix[,i] = Time.Ps[[i]]
    # make Time.Ps list into a matrix with rows for
    # permutations and columns for time-points

    MinPVect = apply(PermsMatrix, 1, min)

    MinP = quantile(MinPVect, probs=alpha)


    # EWp = ????


    EWQuants = as.vector(
      tapply(SumPerms, INDEX=PermTLevs,
             FUN=function(z) {
               quantile(z, probs=1-MinP)
             }
             )
      )


    ## Bonferroni CIs ##

    cat(paste(" done", "\n", "Bonferroni CIs ..."))

    BOQuants = as.vector(
      tapply(SumPerms, INDEX=PermTLevs,
             FUN=function(z) {
               quantile(z, probs=1-(alpha/LenTLevs))
             }
             )
      )

    cat(paste(" done", "\n"))


  } else if (PermT$stat == "t") {


  ##... Raw t ...##

    SumPerms = PermT$perms.sumt
    SumObs = PermT$obs.sumt


    ##... Pointwise CIs ...##

    cat(paste(" Pointwise CIs ..."))

    PWp = rep(0, LenTLevs)
    for (t in 1:LenTLevs) {
      PWp[t] =
        sum(abs(SumPerms[PermTLevs==TimeLevs[t]]) >=
             abs(SumObs[TimeLevs==TimeLevs[t]])) / nperms
    }

    PWQuants = tapply(SumPerms, INDEX=PermTLevs,
      FUN=function(z){
        quantile(z, probs=alpha/2) # two-tailed test
      }
      )


    ##... Exptwise CIs ...##

    cat(paste(" done", "\n", "Experimentwise CIs ..."))

    Time.Ps = tapply(SumPerms, INDEX=PermTLevs, simplify=F,
      FUN=function(z) {
        ranks = rank(z, ties.method="min")
        rank = (nperms-ranks)+1
        # method "min" means that nperms - ranks + 1 gives
        # number of values greater than or equal to
        # each permutation sum of t within time-points
        pvals[rank]
      }
      )

    PermsMatrix = matrix(0, nrow=nperms, ncol=LenTLevs)
    for (i in 1:LenTLevs) PermsMatrix[,i] = Time.Ps[[i]]

    MinPVect = apply(PermsMatrix, 1, min)

    MinP = quantile(MinPVect, probs=alpha)


    # EWp = ????


    EWQuants = as.vector(
      tapply(SumPerms, INDEX=PermTLevs,
             FUN=function(z) {
               quantile(z, probs=MinP)
             }
             )
      )


    ## Bonferroni CIs ##

    cat(paste(" done", "\n", "Bonferroni CIs ..."))

    BOQuants = as.vector(
      tapply(SumPerms, INDEX=PermTLevs,
             FUN=function(z) {
               quantile(z, probs=alpha/LenTLevs)
             }
             )
      )

    cat(paste(" done", "\n"))

  }

  PermQuants = list(
    PWQuants=PWQuants, PWp=PWp,
    EWQuants=EWQuants, EWp=MinP,
    BOQuants = BOQuants, BOp=alpha/LenTLevs,
    TimeLevs=TimeLevs, nperms=nperms, alpha=alpha)

    return(PermQuants)
}


##---------------------///////---------------------##
##
## Plot observed and critical sum of t
## Add bands where significant
##
## Must have a EEGPermT and EEGtquant objects
##
## 5 October
##
##---------------------///////---------------------##


EEGtqplot2 = function(PermT, QuantsT,
  data.name="", sub="") {

    SumObs = PermT$obs.sumt
    TimeLevs = unique(PermT$Times)

    stat = PermT$stat
    alpha = QuantsT$alpha

    BOQuants = QuantsT$BOQuants
    PWQuants = QuantsT$PWQuants
    EWQuants = QuantsT$EWQuants


  ##... Absolute t ...###

  if (stat == "abs.t") {

    par(cex=0.8)

    ymax=max(BOQuants)+50
    ymin=min(SumObs)-50

    plot(x=TimeLevs, SumObs, type="l", lwd=2,
         ylim=c(ymin, ymax),
         xlab="time (ms)", ylab=paste("sum of", stat))
    abline(v=0, lty=2, col="blue")

    lines(x=TimeLevs, PWQuants, col="grey", lty=2)

    lines(x=TimeLevs, EWQuants, col="red", lty=5)

    lines(x=TimeLevs, BOQuants, col="green3", lty=4)

    diffs = EWQuants-SumObs
    sig.diff = diffs

    if (sum(sig.diff<=0) > 0 ) {

      sig.times = TimeLevs[sig.diff<=0]
      sig.obs = SumObs[sig.diff<=0]
      for(i in 1:length(sig.times)) {
        lines(x=c(sig.times[i], sig.times[i]),
              y=c(-10,sig.obs[i]))
      }
    }


    title(main=paste("Sum of ", stat, " with ", alpha*100,
            "% Critical Values", "\n", data.name, sep=""),
          sub=sub)

    legend(10, ymax, ncol=4,
           lwd=c(2, 1, 1, 1),
           lty=c(1, 2, 5, 4),
           col=c("black", "grey", "red", "green3"),
           legend=c(paste("sum of", stat), "Point-wise CVs",
             "Family-wise CVs", "Bonferroni CVs"
             ))
  }


  else if (stat == "t") {


  ##... Raw t ...##

    par(cex=0.8)

    ymax=max(-BOQuants)+20
    ymin=min(BOQuants)-20

    plot(x=TimeLevs, SumObs, type="l", lwd=2,
         ylim=c(ymin, ymax),
         xlab="time (ms)", ylab="sum of t")
    abline(v=0, lty=2, col="blue")
    abline(h=0, lty=2, col="blue")

    lines(x=TimeLevs, PWQuants, col="grey", lty=2)
    lines(x=TimeLevs, -PWQuants, col="grey", lty=2)

    lines(x=TimeLevs, EWQuants, col="red", lty=5)
    lines(x=TimeLevs, -EWQuants, col="red", lty=5)

    lines(x=TimeLevs, BOQuants, col="green3", lty=4)
    lines(x=TimeLevs, -BOQuants, col="green3", lty=4)

    title(main=paste("Sum of ", stat, " with ", alpha*100,
            "% Critical Values", "\n", data.name, sep=""))

    legend(10, ymax, ncol=4,
           lty=c(1, 2, 5, 4), lwd=c(2, 1, 1, 1),
           col=c("black", "grey", "red", "green3"),
           legend=c("sum t", "Point-wise CIs",
             "Expt-wise CIs", "Bonferroni CIs"
             ))
  }

}



### * 2. EXAMPLE
################################################################################

### ** Dummy data
################################################################################

## This isn't supposed to represent anything in particular. I made this data set
## simply to demonstrate that the code runs without error.

## Labels for times, treatments, listening condition and subjects ('Name').
dummy.data.labels <-
    expand.grid(Time = seq(from = 0, to = 700, by = 4)
                ,Treatment = c("N", "T")
                ,Name = letters[1:16]
                )

## AR(1) dummy data
set.seed(21)
dummy.data.electrodes <-
    replicate(n = 129
              ,expr = arima.sim(model = list(ar=0.99), n = nrow(dummy.data.labels)) + 175
              )

## Add something to timepoints 40--60, for Treatment=="N".
index1 <-
    dummy.data.labels$Time > 39 & dummy.data.labels$Time < 61 &
      dummy.data.labels$Treatment == "N"
dummy.data.electrodes[index1,] <-
    dummy.data.electrodes[index1,] + runif(n = sum(index1), min = 1, max = 3.5)

## Check
plot(dummy.data.electrodes[order(dummy.data.labels$Time), 4]
     ,type = "l")

## Combine to create complete dataset
dummy.data <- cbind(dummy.data.electrodes, dummy.data.labels)

## Actually need /differences/ within participants
txN <- dummy.data$Treatment == "N"
txT <- dummy.data$Treatment == "T"
dummy.data <- data.frame(dummy.data[txN, 1:129] - dummy.data[txT, 1:129],
                         Time = dummy.data$Time[txN],
                         Name = dummy.data$Name[txN],
                         check.names = FALSE)


### ** Do the tests and plot
################################################################################

set.seed(22)

## This takes a /long/ time so for example use only a subset of the time points.
time_subset <- seq(from = 0, to = 100, by = 4)
dummy.data <- dummy.data[dummy.data$Time %in% time_subset, ]


### t-test
###

Ttest <- EEGttest2(dummy.data)


### Quantiles and Adjusted p-values
###

## Quantiles are in elements $PWQuants, $EWQuants, $BOQuants.
## Adjusted p-values are in elements $PWp, $EWp, $BOp.

tquants <- EEGtquant2(Ttest)


### Plots
###

EEGtqplot2(Ttest, tquants)
