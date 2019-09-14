### Export Mike's colors
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")


#' Remove the trailing lp__ from a vector or list
#'
#' @param x A Vector.
#'
#' @return The same vector without the trailing "lp__", if present
#' @export
#'
#' @examples
remove_lp__ <- function(x) {
  if ("lp__" %in% names(x))
    x[names(x) != "lp__"]
}

##' Take a vector, bin it and compute the frequencies
##'
##' This function takes a vector of numbers and bins them in bins defined by
##' the \code{br} argument.
##' @title Bin a vector and compute the frequencies
##' @param x Numeric. A vector of values.
##' @param br Numeric. A vector of bin breaks.
##' @return A vector of bin counts.
##' @author Giovanni d'Ario
bin_and_count <- function(x, B, br) {
  hist(x, breaks = br, plot = FALSE)$counts
}

##' Plot the quantiles of histogram bins.
##'
##' Given a numeric matrix, this function computes the bin frequencies for
##' each row, and plots the resulting distribution.
##' @title Plot the quantiles of histogram bins.
##' @param simu_ys Numeric matrix with individual observations in the rows.
##' @param Numeric. The largest value for for hte histogram range.
##' @param br Numeric. Vector of breaks.
##' @return No return value. A plot is produced as a side effect.
##' @author Giovanni d'Ario
plot_bin_quantiles <- function(simu_ys,
                               B=max(simu_ys),
                               br=0:(B + 1) - 0.5) {
  ## Vector of quantiles
  pbs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  ## Compute the bin frequency for each row of simu_ys
  bin_counts <- apply(simu_ys, 1, bin_and_count, B = B, br = br)
  ## Compute the quantiles of each bin across the columns
  quants <- apply(bin_counts, 1, quantile, probs = pbs)
  df_quants <- as_tibble(t(quants))
  names(df_quants) <- paste("q", as.character(10 * pbs), sep = "_")
  df_quants$x <- 0:B
  ## Plot the bin quantiles
  ggplot(df_quants, aes(x = x)) +
    geom_ribbon(aes(ymin = q_1, ymax = q_9), fill = c_light) +
    geom_ribbon(aes(ymin = q_2, ymax = q_8), fill = c_light_highlight) +
    geom_ribbon(aes(ymin = q_3, ymax = q_7), fill = c_mid) +
    geom_ribbon(aes(ymin = q_4, ymax = q_6), fill = c_mid_highlight) +
    geom_line(aes(y = q_5)) +
    xlab("y") + ylab("") +
    ggtitle("Prior Predictive Distribution")
}

create_simulation_matrix <- function(params) {
  t(data.matrix(data.frame(params[-length(params)])))
}

get_sbc_rank <- function(sim_param, param, thin) {
  sum(sim_param < param[seq(1, length(param) - thin, thin)])
}

get_z_score <- function(post_mean, simul, post_sd) {
  (post_mean - simul) / post_sd
}

get_shrinkage <- function(post_sd, prior_sd) {
  1 - (post_sd / prior_sd) ** 2
}

get_sbc_parameters <- function(sim_vector,
                               model,
                               param_names=NULL,
                               prior_sd=NULL,
                               thin=8,
                               N,
                               seed=4938483) {

  ## Consistency checks
  stopifnot(!is.null(param_names))
  stopifnot(!is.null(prior_sd))
  n_params <- length(param_names)
  stopifnot(length(prior_sd) == n_params)
  # ## Make sure that the parameters appear in the correct order in sim_vector
  # if (!identical(param_names, names(sim_vector)[seq_len(n_params)]))
  #   stop("The parameter names are not matching with those in `sim_vector`")

  ## The simulated ys come after the n_params parameters
  sim_params <- sim_vector[seq_len(n_params)]
  sim_y <- sim_vector[-seq_len(n_params)]
  input_data <- list("N" = N, "y" = sim_y)

  ## Sample from the model. These values will be compared with those stored in
  ## sim_params
  capture.output(
    fit <- sampling(model, data = input_data, seed = seed)
  )

  ## Check if there are problematic cases
  util <- new.env()
  source("stan_utility.R", local = util)
  warning_code <- util$check_all_diagnostics(fit, quiet = TRUE)

  ## Simulation Based Calibration ranking
  params <- rstan::extract(fit)
  params <- remove_lp__(params)

  ## List of SBC ranks for each parameter
  sbc_ranks <- purrr::map2(.x = as.list(sim_params),
                           .y = params,
                           .f = get_sbc_rank,
                           thin = thin)
  names(sbc_ranks) <- paste("rank", names(sbc_ranks), sep = "_")

  ## Compute the posterior mean and std.dev for each parameter
  summaries <- summary(fit, probs = c(), pars = param_names)$summary
  post_mean <- summaries[, 1]
  post_sd <- summaries[, 3]

  ## Compute posterior z score and posterior shrinkage
  z_score <- get_z_score(post_mean,
                         sim_vector[seq_len(n_params)],
                         post_sd)
  names(z_score) <- paste("z_score", names(z_score), sep = "_")

  shrinkage <- get_shrinkage(post_sd, prior_sd)
  names(shrinkage) <- paste("shrinkage", names(shrinkage), sep = "_")

  c("warning_code" = warning_code, unlist(sbc_ranks), z_score, shrinkage)
}
