# Synthetic audit cohorts. simulate_cell gives the audited value as the outcome; simulate_covariate_cell gives it as a predictor.
config <- list(simulation_parameters = list(b0 = 2L, b1 = 0.6, size_log_mean = 3L, size_log_sd = 0.5,      x_noise_sd = 3L, cluster_variance = 20L, row_variance = 80L,      stratum_sd = 2L))
DESIGNS <- list(rwanda = c(14L, 9L, 8L, 13L, 5L, 9L, 11L), pairs = rep(8L, 
    10L), none = 70L, none_210 = 210L)
simulate_cell <- function(cell, seed, fraction) {
    set.seed(seed)
    a <- config$simulation_parameters
    sizes <- DESIGNS[[cell$design]]
    rows <- cell$rows
    g <- sum(sizes)
    n <- g * rows
    stratified <- length(sizes) > 1L
    id <- rep(seq_len(g), each = rows)
    stratum_of <- rep(seq_along(sizes), sizes)
    row_stratum <- stratum_of[id]
    size <- rep(stats::rlnorm(g, a$size_log_mean, a$size_log_sd), 
        each = rows)
    x <- size + stats::rnorm(n, 0, a$x_noise_sd)
    if (cell$errors == "rare") {
        level <- if (stratified) 
            stats::runif(length(sizes), 0, 0.2)[row_stratum]
        else 0
        frailty <- rep(stats::rgamma(g, 0.5, 0.5), each = rows)
        truth <- as.numeric(stats::rpois(n, pmax(0.01, 0.1 + 
            0.02 * x + level) * frailty))
        slip <- stats::runif(n) < 0.2
        recorded <- pmax(0, truth + slip * sample(c(-1, 1), n, 
            replace = TRUE))
        target <- 0.02
    }
    else {
        level <- if (stratified) 
            stats::rnorm(length(sizes), 0, a$stratum_sd)[row_stratum]
        else 0
        residual <- if (cell$errors == "row_only") 
            stats::rnorm(n, 0, sqrt(a$cluster_variance + a$row_variance))
        else rep(stats::rnorm(g, 0, sqrt(a$cluster_variance)), 
            each = rows) + stats::rnorm(n, 0, sqrt(a$row_variance))
        truth <- a$b0 + a$b1 * x + level + residual
        recorded <- switch(cell$errors, shift_low = truth + rep(stats::rnorm(g, 
            0, 3), each = rows) + stats::rnorm(n, 0, 6), shift = truth + 
            rep(stats::rnorm(g, 0, 5), each = rows) + stats::rnorm(n, 
            0, 3), under = (1 - rep(stats::rbeta(g, 2, 6), each = rows)) * 
            truth + stats::rnorm(n, 0, 2), outlier = truth + 
            rep((stats::runif(g) < 0.1) * stats::rnorm(g, 0, 
                20), each = rows) + stats::rnorm(n, 0, 3), row_only = truth + 
            stats::rnorm(n, 0, 6), stop("Unknown errors"))
        target <- a$b1
    }
    chosen <- logical(g)
    start <- 0L
    for (nh in sizes) {
        chosen[start + sample.int(nh, max(1L, as.integer(round(fraction * 
            nh))))] <- TRUE
        start <- start + nh
    }
    audited <- chosen[id]
    gold <- truth
    gold[!audited] <- NA_real_
    input <- data.frame(covariate = x, cluster_id = id, lab_value = recorded, 
        gold_value = gold)
    if (stratified) 
        input$stratum <- factor(paste0("s", row_stratum))
    list(input = input, truth = truth, target = target, audited = audited, 
        strata = if (stratified) "stratum" else NULL, formula = if (stratified) lab_value ~ 
            covariate + stratum else lab_value ~ covariate, predictors = if (stratified) c("covariate", 
            "stratum") else "covariate")
}
simulate_covariate_cell <- function(cell, seed, fraction) {
    p <- simulate_cell(cell, seed, fraction)
    d <- p$input
    id <- d$cluster_id
    h <- as.integer(d$stratum)
    d$outcome <- 1 + 0.5 * p$truth + 0.3 * d$covariate + stats::rnorm(nlevels(d$stratum))[h] + 
        rep(stats::rnorm(max(id), 0, sqrt(5)), each = cell$rows) + 
        stats::rnorm(nrow(d), 0, sqrt(20))
    list(input = d, truth = p$truth, target = 0.5, audited = p$audited, 
        formula = outcome ~ lab_value + covariate + stratum, 
        predictors = c("covariate", "stratum", "outcome"), term = "lab_value", 
        imputed_outcome = FALSE)
}
