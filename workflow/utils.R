bin_and_count <- function(x, B, br) {
  hist(x, breaks = br, plot = FALSE)$counts
}

plot_bin_quantiles <- function(simu_ys, 
                               B=40, 
                               br=0:(B + 1) - 0.5) {
  pbs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  ## Compute the bin frequency for each row of simu_ys
  bin_counts <- apply(simu_ys, 1, bin_and_count, B = B, br = br)
  ## Compute the quantiles of each bin across the columns 
  quants <- apply(bin_counts, 1, quantile, probs = pbs)
  df_quants <- as_tibble(t(quants))
  names(df_quants) <- paste("q", as.character(10 * pbs), sep = "_")
  df_quants$x <- 0:40
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

get_sbc_rank <- function(param, thin) {
  sum(param < param[seq(1, length(param) - thin, thin)])
}

get_z_score <- function(post_mean, simul, post_sd) {
  (post_mean - simul) / post_sd
}

get_shrinkage <- function(post_sd, prior_sd) {
  1 - (post_sd / prior_sd) ** 2
}

get_sbc_parameters <- function(simu_list,
                               model,
                               params=NULL,
                               thin=8,
                               N,
                               prior_sd_lambda=5.82337,
                               seed=4938483) {

  stopifnot(!is.null(params))
  n_params <- length(params)
  ## The simulated ys come after the n_params parameters
  simu_y <- simu_list[-seq_len(n_params)]
  input_data <- list("N" = N, "y" = simu_y)
  
  capture.output(
    fit <- sampling(model, data = input_data, seed = seed)
  )

  ## Check if there are problematic cases
  util <- new.env()
  source("stan_utility.R", local = util)
  warning_code <- util$check_all_diagnostics(fit, quiet = TRUE)

  ## Simulation Based Calibration ranking
  params <- extract(fit)
  sbc_ranks <- purrr::map(params, get_sbc_rank, thin = thin)
  summaries <- purrr::map(params, ~summary(fit, probs = c(), pars = .)$summary)
  post_mean <- purrr::map(summaries, extract(.[, 1]))
  post_sd <- purrr::map(summaries, extract(.[, 3]))

  ## Compute posterior z score and posterior shrinkage
  z_score <- purrr::pmap(list(post_mean, simu_list, post_sd))
  # shrinkage <- 1 - (post_sd_lambda / prior_sd_lambda)**2

  # c(warning_code, sbc_rank, z_score, shrinkage)
}
