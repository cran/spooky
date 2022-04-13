#' spooky
#'
#' @param df A data frame with time features on columns
#' @param seq_len Positive integer. Time-step number of the forecasting sequence. Default: NULL (automatic selection between 1 and the square root of full length).
#' @param lno Positive integer. Number of data points to leave out for resampling (using jack-knife approach). Default: NULL (automatic selection between 1 and the square root of full length).
#' @param n_samp Positive integer. Number of samples for random search. Default: 30.
#' @param smoother Logical. Flag to TRUE for loess smoothing. Default: FALSE.
#' @param n_windows Positive integer. Number of validation windows to test prediction error. Default: 10.
#' @param ci Confidence interval for prediction. Default: 0.8
#' @param dates Date. Vector with dates for time features.
#' @param error_scale String. Scale for the scaled error metrics. Two options: "naive" (average of naive one-step absolute error for the historical series) or "deviation" (standard error of the historical series). Default: "naive".
#' @param error_benchmark String. Benchmark for the relative error metrics. Two options: "naive" (sequential extension of last value) or "average" (mean value of true sequence). Default: "naive".
#' @param seed Positive integer. Random seed. Default: 42.
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @return This function returns a list including:
#' \itemize{
#' \item exploration: list of all not-null models, complete with predictions, test metrics, prediction stats and plot
#' \item history: a table with the sampled models, hyper-parameters, validation errors, weighted average rank
#' \item best_model: results for the best selected model according to the weighted average rank, including:
#' \itemize{
#' \item testing_errors: testing errors for each time feature for the best selected model (me, mae, mse, rmsse, mpe, mape, rmae, rrmse, rame, mase, smse, sce, gmrae)
#' \item preds: min, max, q25, q50, q75, quantiles at selected ci, mean, sd, mode, skewness, kurtosis, IQR to range, risk ratio, upside probability and divergence for each point fo predicted sequences
#' \item plots: standard plot with confidence interval for each time feature
#' }
#' \item time_log
#' }
#'
#' @export
#'
#' @import purrr
#' @import tictoc
#' @importFrom readr parse_number
#' @importFrom lubridate seconds_to_period is.Date as.duration
#' @importFrom stats weighted.mean
#' @importFrom imputeTS na_kalman
#' @importFrom fANCOVA loess.as

#'@examples
#'spooky(time_features, seq_len = c(10, 30), lno = c(1, 30), n_samp = 1)
#'

spooky <- function(df, seq_len = NULL, lno = NULL, n_samp = 30, n_windows = 3, ci = 0.8, smoother = FALSE, dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42)
{

  tic.clearlog()
  tic("time")

  set.seed(seed)

  n_length <- nrow(df)
  n_feats <- ncol(df)
  feat_names <- colnames(df)

  if(anyNA(df)){df <- as.data.frame(purrr::map(df, ~ na_kalman(.x))); message("kalman imputation on target and/or regressors\n")}
  if(smoother==TRUE){df <- as.data.frame(purrr::map(df, ~ suppressWarnings(loess.as(x=1:n_length, y=.x)$fitted))); message("performing optimal smoothing\n")}

  if(length(seq_len)==1 & length(lno)==1){n_samp <- 1}
  if(is.null(seq_len)){seq_len <- c(1, round(sqrt(n_length)))}
  if(is.null(lno)){lno <- c(1, round(sqrt(n_length)))}
  sqln_set <- round(runif(n_samp, min(seq_len), max(seq_len)))
  lno_set <- round(runif(n_samp, min(lno), max(lno)))

  exploration <- purrr::map2(sqln_set, lno_set, ~ windower(df, seq_len = .x, lno = .y, n_windows, ci, dates, error_scale, error_benchmark))
  errors <- round(Reduce(rbind, purrr::map(exploration, ~ .x$max_errors)), 2)
  if(n_samp == 1){errors <- t(as.data.frame(errors))}
  colnames(errors) <- paste0("max_", colnames(errors))
  history <- as.data.frame(cbind(seq_len = sqln_set, lno = lno_set, errors))
  rownames(history) <- NULL

  if(n_samp > 1){
    weights <- apply(errors, 2, function(x) {abs(sd(x[is.finite(x)], na.rm = TRUE)/mean(x[is.finite(x)], na.rm = TRUE))})
    finite_w <- is.finite(weights)
    history$wgt_avg_rank <- round(apply(apply(abs(errors[, finite_w, drop = FALSE]), 2, rank), 1, weighted.mean, w = weights[finite_w]), 2)
    history <- history[order(history$wgt_avg_rank),]
    best_index <- as.numeric(rownames(history[1,]))
    best_model <- exploration[[best_index]]
  }

  if(n_samp == 1){history$wgt_avg_rank <- 1; best_model <- exploration[[1]]}

  if(n_feats > 1){testing_errors <- Reduce(rbind, best_model$errors)}
  if(n_feats == 1){testing_errors <- t(as.data.frame(best_model$errors))}
  rownames(testing_errors) <- feat_names
  preds <- purrr::map(best_model$model, ~ .x$quantile_preds)
  names(preds) <- feat_names
  plots <- purrr::map(best_model$model, ~ .x$plot)
  names(plots) <- feat_names

  best_model <- list(testing_errors = testing_errors, preds = preds, plots = plots)

  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best_model = best_model, time_log = time_log)
  return(outcome)
}
