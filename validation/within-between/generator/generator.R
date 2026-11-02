# generator.R -- data generator for a methods test of clustered audit data.
# Base R and stats only.
#
# simulate_cell(cell, seed, fraction) makes one cohort from one row of grid.csv.
# Generation order: design -> covariate -> truth -> outcome -> recording errors
# -> audit draw. The audit is drawn last, so a given seed yields the same truth,
# outcome and recorded values at both audit fractions.
#
# The target depends only on the truth and outcome steps. Recording errors and
# the audit never change it.

ROW_CAP <- 1500L
SIZE_MIN <- 2L
SIZE_MAX <- 30L

# Standardized noise: mean 0, variance 1, whatever the shape.
draw_std <- function(n, tail) {
  switch(tail,
    normal    = rnorm(n),
    t3        = rt(n, 3) / sqrt(3),            # no finite kurtosis: rare huge values
    t4        = rt(n, 4) / sqrt(2),
    skew      = rgamma(n, 1) - 1,              # exponential shape, one-sided outliers
    lognormal = (exp(rnorm(n)) - exp(0.5)) / sqrt((exp(1) - 1) * exp(1)),
    stop("unknown tail: ", tail))
}

cluster_mean <- function(v, cl, sizes) as.vector(rowsum(v, cl, reorder = FALSE)) / sizes

# Clusters per stratum are fixed by design. A registry's list of facilities by
# type does not change between cohorts, only which facility sits in which type.
stratum_counts <- function(K, H, skew) {
  w <- if (skew > 0) seq_len(H) else rep(1, H)
  raw <- K * w / sum(w)
  n <- floor(raw)
  left <- K - sum(n)
  if (left > 0) {
    top <- order(raw - n, decreasing = TRUE)[seq_len(left)]
    n[top] <- n[top] + 1
  }
  n
}

# Cluster sizes.
#  balanced  : a fixed form per site, a fixed number of charts pulled.
#  gamma     : ordinary spread in clinic volume.
#  lognormal : volume across facilities is typically right-skewed, many small.
#  two_tier  : many health posts and a few district hospitals.
draw_sizes <- function(p, h) {
  K <- p$K
  centre <- h - mean(seq_len(max(h)))
  mu <- p$size_mean * exp(p$size_shift * centre)
  s <- switch(p$size_type,
    balanced  = rep(p$size_mean, K),
    gamma     = rgamma(K, shape = 1 / p$size_cv^2, scale = mu * p$size_cv^2),
    lognormal = mu * exp(rnorm(K, -p$size_cv^2 / 2, p$size_cv)),
    two_tier  = ifelse(runif(K) < 0.15, runif(K, 20, 31), runif(K, 2, 7)),
    stop("unknown size_type: ", p$size_type))
  s <- as.integer(pmin(pmax(round(s), SIZE_MIN), SIZE_MAX))
  excess <- sum(s) - ROW_CAP
  while (excess > 0L) {
    i <- which.max(s)
    s[i] <- s[i] - 1L
    excess <- excess - 1L
  }
  s
}

# Covariate. z is the internal scale; the reported covariate is x_mean + x_sd * z
# (binary stays 0/1).
#  x_icc 1           : a site attribute, such as staffing level or ownership.
#  x_icc 0, centered : varies only within site, such as a within-site allocation.
#  binary            : sex, treated yes/no, with site-specific mix.
#  skew              : age or comorbidity-like, one long side.
draw_covariate <- function(p, cl, sizes, h, zs) {
  K <- p$K
  N <- length(cl)
  if (p$x_type == "binary") {
    if (p$x_icc >= 1) {
      arm <- numeric(K)
      for (s in unique(h)) {
        idx <- which(h == s)
        arm[idx] <- sample(rep_len(0:1, length(idx)))
      }
      return(arm[cl])
    }
    if (p$x_center > 0) {
      ord <- order(cl, runif(N))
      z <- numeric(N)
      z[ord] <- as.numeric(sequence(sizes) <= floor(sizes[cl] / 2))
      return(z)
    }
    if (p$x_icc <= 0) return(as.numeric(runif(N) < 0.5))
    conc <- 1 / p$x_icc - 1
    pk <- rbeta(K, conc * 0.5, conc * 0.5)
    return(as.numeric(runif(N) < pk[cl]))
  }
  tail <- if (p$x_type == "skew") "skew" else "normal"
  z <- sqrt(p$x_icc) * draw_std(K, tail)[cl] + sqrt(1 - p$x_icc) * draw_std(N, tail)
  if (p$x_center > 0) z <- z - cluster_mean(z, cl, sizes)[cl]
  z + (p$x_size * zs + p$x_strat * (h - mean(seq_len(max(h)))))[cl]
}

# True value: stratum level + slope * z + contextual + curvature + site + row.
draw_truth <- function(p, z, cl, sizes, h, zs) {
  K <- p$K
  N <- length(cl)
  zbar <- cluster_mean(z, cl, sizes)
  u <- p$sd_u * draw_std(K, p$tail_u) + p$u_size * zs
  slope <- p$b_x + p$hslope * (h - mean(seq_len(max(h)))) + p$s_slope * rnorm(K)
  nominal <- p$b0 + (p$hfx * (h - 1) + p$b_ctx * zbar + u)[cl] + slope[cl] * z +
    p$b_x2 * z^2 + p$nom_sd * rnorm(N)
  e <- p$sd_e * draw_std(N, p$tail_e) * (1 + p$e_het * abs(z))
  truth <- switch(p$truth_type,
    linear    = nominal + e,
    lognormal = exp(nominal + e),              # costs, biomarkers: multiplicative
    count     = as.numeric(rpois(N, exp(nominal + e))),  # monthly register totals
    stop("unknown truth_type: ", p$truth_type))
  list(truth = truth, nominal = nominal)
}

# Outcome for the predictor role, a function of the TRUE value. v is returned
# because a site's quality can also bend its recording.
draw_outcome <- function(p, truth, z, cl, sizes, h, zs) {
  K <- p$K
  N <- length(cl)
  g <- switch(p$y_link,
    identity = truth,
    log      = log(truth),                     # dose response on the log scale
    log1p    = log1p(truth),
    step     = as.numeric(truth > p$y_thr),    # a clinical cut point
    stop("unknown y_link: ", p$y_link))
  tbar <- cluster_mean(truth, cl, sizes)
  v <- p$sd_v * draw_std(K, p$tail_v) + p$v_size * zs
  gk <- p$gamma + p$g_slope * rnorm(K)
  y <- p$c0 + (p$h_y * (h - 1) + p$ctx_truth * tbar + v)[cl] + gk[cl] * g +
    p$delta * z + p$sd_eps * draw_std(N, p$tail_eps)
  list(y = y, v = v)
}

# Recording process. Each block is a mechanism; a cell turns on one or several.
draw_recorded <- function(p, tr, y, v, z, cl, sizes, h, zs) {
  K <- p$K
  N <- length(cl)
  truth <- tr$truth
  mt <- mean(truth)
  st <- sd(truth)

  # Berkson: the register holds the assigned or nominal value (prescribed dose,
  # kit lot value), and the true value scatters around it.
  lab <- if (p$berkson > 0) tr$nominal else truth

  # Value-dependent shrinkage: recall and self-report pull high values down and
  # low values up. Negative err_value compresses toward the mean.
  if (p$err_value != 0) lab <- lab + p$err_value * (truth - mt)

  # Site calibration slope: each analyser or lab reads a fixed percent off, so
  # the error grows with the value and is shared within site.
  if (p$prop_sd != 0 || p$prop_mean != 0)
    lab <- lab * (1 + p$prop_mean + p$prop_sd * rnorm(K))[cl]

  # Multiplicative row noise: billed amounts and concentrations err in percent.
  if (p$mult_sd != 0) lab <- lab * exp(p$mult_sd * rnorm(N))

  # Site offset: uncorrected instrument bias, a local reference range, a local
  # coding habit. err_size lets busier sites differ. err_v ties it to the site
  # quality that also moves the outcome.
  if (p$clu_sd != 0 || p$clu_mean != 0 || p$err_size != 0 || p$err_v != 0) {
    offset <- p$clu_mean + p$clu_sd * draw_std(K, p$tail_clu) + p$err_size * zs +
      if (length(v) == K) p$err_v * v else 0
    lab <- lab + offset[cl]
  }

  # Bias by covariate: older patients, or sicker or busier shifts, get recorded
  # differently. err_x is common; err_x_clu differs by site (a clinic that
  # misreads one subgroup).
  if (p$err_x != 0) lab <- lab + p$err_x * z
  if (p$err_x_clu != 0) lab <- lab + p$err_x_clu * rnorm(K)[cl] * z

  # Differential error: the person recording knows how the patient did, and
  # severity colours the entry. Only meaningful with an outcome.
  if (p$err_diff != 0) lab <- lab + p$err_diff * (y - mean(y)) / sd(y)

  # Within-run drift: reagent or calibration ages over a run, at a rate specific
  # to the site, reset between runs. Row order in the cluster stands for time.
  if (p$drift_sd != 0) {
    pos <- (sequence(sizes) - 1) / pmax(sizes[cl] - 1, 1)
    lab <- lab + p$drift_sd * rnorm(K)[cl] * pos
  }

  # Row noise: assay imprecision and transcription jitter. noise_het_value makes
  # it grow with the value (meters are worse at high readings). noise_clu_sd
  # makes some sites sloppier than others.
  if (p$noise_sd != 0) {
    w <- 1 + p$noise_het_value * abs(truth - mt) / st
    if (p$noise_clu_sd != 0) w <- w * exp(p$noise_clu_sd * rnorm(K))[cl]
    lab <- lab + p$noise_sd * draw_std(N, p$tail_noise) * w
  }

  # One bad stratum: for example private facilities outside the reporting
  # support system have both bias and noise.
  if (p$err_s1_mean != 0 || p$err_s1_sd != 0) {
    is1 <- h[cl] == 1
    lab[is1] <- lab[is1] + p$err_s1_mean + p$err_s1_sd * rnorm(sum(is1))
  }

  # Under-reporting: each event reaches the summary form with a site-specific
  # completeness.
  if (p$compl_mean > 0) {
    pc <- plogis(qlogis(p$compl_mean) + p$compl_sd * rnorm(K))[cl]
    lab <- as.numeric(rbinom(N, round(lab), pc))
  }

  # Over-reporting: where pay follows volume, some sites inflate, more so where
  # the covariate is high (inflate_x).
  if (p$inflate_share > 0) {
    bad <- runif(K) < p$inflate_share
    lab <- lab * (1 + p$inflate * bad[cl] * exp(p$inflate_x * z))
    if (p$truth_type == "count") lab <- round(lab)
  }

  # Threshold gaming: values just short of a payment or eligibility cut get
  # nudged over it at some sites.
  if (p$gaming_share > 0) {
    bad <- (runif(K) < p$gaming_share)[cl]
    near <- lab < p$gaming_thr & lab >= p$gaming_thr - p$gaming_band
    hit <- bad & near & runif(N) < p$gaming_prob
    lab[hit] <- p$gaming_thr + runif(sum(hit), 0, p$gaming_band / 4)
  }

  # Assay limits: below the limit of detection, the lab reports the limit.
  if (p$lod > 0) lab <- pmax(lab, p$lod)
  if (p$ceil > 0) lab <- pmin(lab, p$ceil)

  # Fabrication: some sites fill forms at the desk with believable values that
  # bear no relation to the patient.
  if (p$fab_share > 0) {
    bad <- which((runif(K) < p$fab_share)[cl])
    lab[bad] <- mt + 0.6 * st * rnorm(length(bad))
  }

  # Unit mix-up: a whole site reports in a different unit (mmol/L against mg/dL).
  if (p$unit_share > 0) {
    bad <- (runif(K) < p$unit_share)[cl]
    lab[bad] <- lab[bad] * p$unit_factor
  }

  # Sparse gross errors: decimal slips, or blanks keyed as zero. The rate is
  # clumped by site, since data clerks differ.
  if (p$gross_rate > 0) {
    rk <- if (p$gross_disp > 0) {
      conc <- 1 / p$gross_disp - 1
      rbeta(K, p$gross_rate * conc, (1 - p$gross_rate) * conc)
    } else rep(p$gross_rate, K)
    hit <- which(runif(N) < rk[cl])
    if (length(hit)) {
      if (p$gross_type == "zero") {
        lab[hit] <- 0
      } else {
        up <- runif(length(hit)) < 0.5
        lab[hit] <- lab[hit] * ifelse(up, p$gross_factor, 1 / p$gross_factor)
      }
    }
  }

  # Digit preference: staff at some sites round to 5s or 10s by hand.
  if (p$round_to > 0) {
    bad <- (runif(K) < p$round_share)[cl]
    lab[bad] <- round(lab[bad] / p$round_to) * p$round_to
  }

  # Copy-forward: a tired clerk repeats the entry on the line above. First row
  # of each site is never copied.
  if (p$copy_rate > 0) {
    conc <- 1 / 0.2 - 1
    rk <- rbeta(K, p$copy_rate * conc, (1 - p$copy_rate) * conc)
    copied <- runif(N) < rk[cl]
    copied[cumsum(sizes) - sizes + 1L] <- FALSE
    lab <- lab[cummax(seq_len(N) * (!copied))]
  }

  lab
}

# Whole clusters audited at random within strata. Allocation is proportional to
# clusters per stratum, tilted toward stratum 1 by audit_tilt (auditors often
# oversample the riskier stratum), at least 1 per stratum, at least 5 overall.
draw_audit <- function(p, h, fraction) {
  K <- p$K
  m <- max(5L, as.integer(floor(fraction * K + 0.5)))
  aud <- logical(K)
  H <- max(h)
  nh <- tabulate(h, H)
  w <- nh * c(p$audit_tilt, rep(1, H - 1L))
  raw <- m * w / sum(w)
  a <- pmin(pmax(1, floor(raw)), nh)
  while (sum(a) < m) {
    room <- which(a < nh)
    i <- room[which.max((raw - a)[room])]
    a[i] <- a[i] + 1
  }
  while (sum(a) > m) {
    can <- which(a > 1)
    i <- can[which.max((a - raw)[can])]
    a[i] <- a[i] - 1
  }
  for (s in seq_len(H)) {
    idx <- which(h == s)
    aud[idx[sample.int(length(idx), a[s])]] <- TRUE
  }
  aud
}

simulate_cell <- function(cell, seed, fraction) {
  if (!fraction %in% c(0.25, 0.5)) stop("fraction must be 0.25 or 0.5")
  set.seed(seed)
  p <- as.list(cell)
  K <- p$K
  has_strata <- p$n_strata > 1

  h <- if (has_strata) {
    sample(rep.int(seq_len(p$n_strata), stratum_counts(K, p$n_strata, p$strata_skew)))
  } else rep(1L, K)

  sizes <- draw_sizes(p, h)
  cl <- rep.int(seq_len(K), sizes)
  ls <- log(sizes)
  zs <- if (sd(ls) > 0) (ls - mean(ls)) / sd(ls) else numeric(K)

  z <- draw_covariate(p, cl, sizes, h, zs)
  tr <- draw_truth(p, z, cl, sizes, h, zs)
  is_outcome <- p$role == "outcome"
  y <- NULL
  v <- 0
  if (!is_outcome) {
    o <- draw_outcome(p, tr$truth, z, cl, sizes, h, zs)
    y <- o$y
    v <- o$v
  }
  lab <- draw_recorded(p, tr, y, v, z, cl, sizes, h, zs)

  # Audit last, so the data above do not depend on the fraction.
  audited <- draw_audit(p, h, fraction)[cl]

  covariate <- if (p$x_type == "binary") z else p$x_mean + p$x_sd * z
  input <- data.frame(covariate = covariate, cluster_id = cl, lab_value = lab,
                      gold_value = ifelse(audited, tr$truth, NA_real_))
  predictors <- "covariate"
  rhs <- if (is_outcome) "covariate" else "lab_value + covariate"
  lhs <- if (is_outcome) "lab_value" else "outcome"
  if (has_strata) {
    input$stratum <- factor(paste0("S", h[cl]), levels = paste0("S", seq_len(p$n_strata)))
    predictors <- c(predictors, "stratum")
    rhs <- paste(rhs, "+ stratum")
  }
  if (!is_outcome) {
    input$outcome <- y
    predictors <- c(predictors, "outcome")
  }
  form <- as.formula(paste(lhs, "~", rhs), env = globalenv())

  list(input = input,
       truth = tr$truth,
       target = if (is.null(p$target)) NA_real_ else p$target,
       audited = audited,
       formula = form,
       predictors = predictors,
       term = if (is_outcome) "covariate" else "lab_value",
       imputed_outcome = is_outcome)
}

# Convenience: read the grid so that strings stay strings.
read_grid <- function(path = "grid.csv") read.csv(path, stringsAsFactors = FALSE)
