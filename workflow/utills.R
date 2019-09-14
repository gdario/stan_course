get_sbc_rank <- function(param, thin) {
  sbc_rank <- sum(param < lambda[seq(1, length(lambda) - thin, thin)])
}

get_z_score <- function(post_mean, simul, post_sd) {
  (post_mean - simul) / post_sd
}

get_shrinkage <- function(post_sd, prior_sd) {
  1 - (post_sd / prior_sd)**2
}

get_sbc_parameters <- function(simu_list,
                               model,
                               thin=8,
                               prior_sd_lambda=5.82337,
                               seed=4938483) {

  capture.output(
    fit <- sampling(model, data = input_data, seed = seed)
  )

  ## Check if there are problematic cases
  util <- new.env()
  warning_code <- util$check_all_diagnostics(fit, quiet = TRUE)

  ## Simulation Based Calibration ranking
  params <- extract(fit)
  sbc_ranks <- purrr::map(params, get_sbc_rank, thin = thin)
  summaries <- purrr::map(params, ~summary(fit, probs = c(), pars = .)$summary)
  post_mean <- purrr::map(summaries, extract(.[, 1]))
  post_sd <- purrr::map(summaries, extract(.[, 3]))

  z_score <- purrr::pmap(list(post_mean, simu_list, post_sd))
  shrinkage <- 1 - (post_sd_lambda / prior_sd_lambda)**2

  c(warning_code, sbc_rank, z_score, shrinkage)
}
