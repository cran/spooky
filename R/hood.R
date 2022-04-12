#' support functions for spooky
#'
#' @param df A data frame with time features on columns
#' @param seq_len Positive integer. Time-step number of the forecasting sequence. Default: NULL (automatic selection between 1 and the square root of full length).
#' @param lno Positive integer. Number of data points to leave out for resampling (using jack-knife approach). Default: NULL (automatic selection between 1 and the square root of full length).
#' @param n_windows Positive integer. Number of validation windows to test prediction error. Default: 10.
#' @param ci Confidence interval for prediction. Default: 0.8
#' @param dates Date. Vector with dates for time features.
#' @param error_scale String. Scale for the scaled error metrics. Two options: "naive" (average of naive one-step absolute error for the historical series) or "deviation" (standard error of the historical series). Default: "naive".
#' @param error_benchmark String. Benchmark for the relative error metrics. Two options: "naive" (sequential extension of last value) or "average" (mean value of true sequence). Default: "naive".
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @import purrr
#' @importFrom modeest mlv1
#' @import ggplot2
#' @importFrom moments skewness kurtosis
#' @importFrom stats quantile sd lm rnorm pnorm fft runif
#' @importFrom scales number
#' @importFrom utils tail head
#' @import greybox

#########
windower <- function(df, seq_len, lno = 1, n_windows = 5, ci = 0.8, dates = NULL, error_scale, error_benchmark)
{
  n_length <- nrow(df)
  n_feats <- ncol(df)
  feat_names <- colnames(df)
  too_short <- floor(n_length/(n_windows + 1)) <= seq_len
  if(too_short){stop("not enough data for the validation windows")}
  idx <- c(rep(1, n_length%%(n_windows + 1)), rep(1:(n_windows + 1), each = n_length/(n_windows + 1)))
  models <- sapply(1:n_feats, function(f) purrr::map(1:n_windows, ~ engine(df[idx <= .x, f, drop = T], seq_len, lno, testing = head(df[idx == .x + 1, f, drop = T], seq_len), ci, dates, feat_names[f], error_scale, error_benchmark)), simplify = F)
  errors <- purrr::map(purrr::map_depth(models, 2, ~ .x$errors), ~ colMeans(Reduce(rbind, .x)))
  errors <- purrr::map(errors, ~ round(.x, 3))
  if(n_feats > 1){max_errors <- apply(Reduce(rbind, errors), 2, max)}
  if(n_feats == 1){max_errors <- unlist(errors)}
  model <- sapply(1:n_feats, function(f) engine(df[, f, drop = T], seq_len, lno, testing = NULL, ci, dates, feat_names[f], error_scale, error_benchmark), simplify = F)
  outcome <- list(errors = errors, max_errors = max_errors, model = model)
  return(outcome)
}

#########
engine <- function(ts, seq_len, lno = 1, testing = NULL, ci = 0.8, dates = NULL, feat_name, error_scale, error_benchmark)
{
  orig <- ts
  diff_model <- recursive_diff(ts, best_deriv(ts))
  ts <- diff_model$vector
  all_positive_check <- all(orig >= 0)
  all_negative_check <- all(orig < 0)
  n_length <- length(ts)
  transformed <- spectral_transformation(ts)
  jacknife_index <- function(v, n) {ifelse(n <= 0, n <- 1, n); n <- n - 1; purrr::map(1:(length(v) - n), ~ v[- (.x:(.x+n))])}
  idx_samples <- jacknife_index(1:length(ts), lno)
  n_trials <- length(idx_samples)
  cmplx <- transformed$complex
  preds <- Reduce(rbind, purrr::map(idx_samples, ~ tail(spectral_inversion(c(cmplx[.x], vector("complex", seq_len))), seq_len)))
  preds <- t(apply(preds, 1, function(x) tail(invdiff(x, diff_model$tail_value), seq_len)))
  if(all_positive_check){preds[preds < 0] <- 0}
  if(all_negative_check){preds[preds > 0] <- 0}

  errors <- NULL
  if(!is.null(testing)){errors <- my_metrics(testing, colMeans(preds), ts, error_scale, error_benchmark)}

  quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))
  p_stats <- function(x){stats <- c(min = suppressWarnings(min(x, na.rm = T)), quantile(x, probs = quants, na.rm = T), max = suppressWarnings(max(x, na.rm = T)), mean = mean(x, na.rm = T), sd = sd(x, na.rm = T), mode = tryCatch(suppressWarnings(modeest::mlv1(x[is.finite(x)], method = "shorth")), error = function(e) NA), kurtosis = tryCatch(suppressWarnings(moments::kurtosis(x[is.finite(x)], na.rm = T)), error = function(e) NA), skewness = tryCatch(suppressWarnings(moments::skewness(x[is.finite(x)], na.rm = T)), error = function(e) NA)); return(stats)}
  quantile_preds <- t(apply(preds, 2, p_stats))
  iqr_to_range <- tryCatch((quantile_preds[, "75%"] - quantile_preds[, "25%"])/(quantile_preds[, "max"] - quantile_preds[, "min"]), error = function(e) NA)
  risk_ratio <- tryCatch((quantile_preds[, "max"] - quantile_preds[, "50%"])/(quantile_preds[, "50%"] - quantile_preds[, "min"]), error = function(e) NA)
  growths <- mapply(function(m, s) rnorm(n_trials, m, s), m = quantile_preds[, "mean"], s = quantile_preds[, "sd"])
  upside_prob <- tryCatch(c(NA, colMeans(apply(growths[,-1]/growths[,-ncol(growths)], 2, function(x) x > 1))), error = function(e) NA)
  pvalues <- mapply(function(m, s) pnorm(seq(min(quantile_preds[, "min"]), max(quantile_preds[, "max"]), length.out = n_trials), m, s), m = quantile_preds[, "mean"], s = quantile_preds[, "sd"])
  divergence <- tryCatch(c(NA, apply(pvalues[,-1] - pvalues[,-ncol(growths)], 2, function(x) abs(max(x, na.rm = T)))), error = function(e) NA)
  quantile_preds <- round(cbind(quantile_preds, iqr_to_range = iqr_to_range, risk_ratio = risk_ratio, upside_prob = upside_prob, divergence = divergence), 3)

  if(is.null(dates)){hist_dates <- 1:length(orig); forcat_dates <- (length(orig) + 1):(length(orig) + seq_len); rownames(quantile_preds) <- paste0("t", 1:seq_len)}
  if(!is.null(dates)){hist_dates <- tail(dates, length(orig)); forcat_dates <- seq.Date(tail(dates, 1), tail(dates, 1) + seq_len * mean(diff(dates)), length.out = seq_len); rownames(quantile_preds) <- as.character(forcat_dates)}
  x_lab <- paste0("Forecasting Horizon for sequence n = ", seq_len)
  y_lab <- paste0("Forecasting Values for ", feat_name)

  plot <- ts_graph(x_hist = hist_dates, y_hist = orig, x_forcat = forcat_dates, y_forcat = quantile_preds[, "50%"], lower = quantile_preds[, "10%"], upper = quantile_preds[, "90%"],
                   label_x = x_lab, label_y = y_lab)

  outcome <- list(errors = errors, quantile_preds = quantile_preds, plot = plot)
  return(outcome)
}


###########
spectral_transformation <- function(ts)
{
  complex <- fft(ts)
  real <- Re(complex)
  imaginary <- Im(complex)
  out <- list(complex = complex, real = real, imaginary = imaginary)
  return(out)
}

########
spectral_inversion <- function(complex = NULL, real, imaginary)
{
  if(is.null(complex)){complex <- complex(real = real, imaginary = imaginary)}
  out <- Re(fft(complex, inverse = T)/length(complex))
  return(out)
}

########
ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                     forcat_band = "darkorange", forcat_line = "darkorange", hist_line = "gray43", label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
{
  all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = c(y_hist, y_forcat))
  forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = y_forcat)

  if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- lower; forcat_data$upper <- upper}

  plot <- ggplot()+geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
  if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes_string(x = "x_forcat", ymin = "lower", ymax = "upper"), alpha = 0.3, fill = forcat_band)}
  plot <- plot + geom_line(data = forcat_data, aes_string(x = "x_forcat", y = "y_forcat"), color = forcat_line, size = line_size)
  if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
  if(is.null(dbreak)){plot <- plot + xlab(label_x)}
  plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)
  plot <- plot + ylab(label_y)  + theme_bw()
  plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

  return(plot)
}


#######
best_deriv <- function(ts, max_diff = 3, thresh = 0.001)
{
  pvalues <- vector(mode = "double", length = as.integer(max_diff))

  for(d in 1:(max_diff + 1))
  {
    model <- lm(ts ~ t, data.frame(ts, t = 1:length(ts)))
    pvalues[d] <- with(summary(model), pf(fstatistic[1], fstatistic[2], fstatistic[3],lower.tail=FALSE))
    ts <- diff(ts)
  }

  best <- tail(cumsum(pvalues < thresh), 1)

  return(best)
}

#####
recursive_diff <- function(vector, deriv)
{
  head_value <- vector("numeric", deriv)
  tail_value <- vector("numeric", deriv)
  if(deriv==0){head_value = NULL; tail_value = NULL}
  if(deriv > 0){for(i in 1:deriv){head_value[i] <- head(vector, 1); tail_value[i] <- tail(vector, 1); vector <- diff(vector)}}
  outcome <- list(vector = vector, head_value = head_value, tail_value = tail_value)
  return(outcome)
}

#####
invdiff <- function(vector, heads)
{
  if(is.null(heads)){return(vector)}
  for(d in length(heads):1){vector <- cumsum(c(heads[d], vector))}
  return(vector)
}

####

jacknife_index <- function(v, n) {ifelse(n <= 0, n <- 1, n); n <- n - 1; purrr::map(1:(length(v) - n), ~ v[- (.x:(.x+n))])}


####

my_metrics <- function(holdout, forecast, actuals, error_scale, error_benchmark)
{
  scale <- switch(error_scale, "deviation" = sd(actuals), "naive" = mean(abs(diff(actuals))))
  benchmark <- switch(error_benchmark, "average" = rep(mean(forecast), length(forecast)), "naive" = rep(tail(actuals, 1), length(forecast)))

  me <- ME(holdout, forecast, na.rm = TRUE)
  mae <- MAE(holdout, forecast, na.rm = TRUE)
  mse <- MSE(holdout, forecast, na.rm = TRUE)
  rmsse <- RMSSE(holdout, forecast, scale, na.rm = TRUE)
  mre <- MRE(holdout, forecast, na.rm = TRUE)
  mpe <- MPE(holdout, forecast, na.rm = TRUE)
  mape <- MAPE(holdout, forecast, na.rm = TRUE)
  rmae <- rMAE(holdout, forecast, benchmark, na.rm = TRUE)
  rrmse <- rRMSE(holdout, forecast, benchmark, na.rm = TRUE)
  rame <- rAME(holdout, forecast, benchmark, na.rm = TRUE)
  mase <- MASE(holdout, forecast, scale, na.rm = TRUE)
  smse <- sMSE(holdout, forecast, scale, na.rm = TRUE)
  sce <- sCE(holdout, forecast, scale, na.rm = TRUE)
  gmrae <- GMRAE(holdout, forecast, benchmark, na.rm = TRUE)

  out <- round(c(me = me, mae = mae, mse = mse, rmsse = rmsse, mpe = mpe, mape = mape, rmae = rmae, rrmse = rrmse, rame = rame, mase = mase, smse = smse, sce = sce, gmrae = gmrae), 3)
  return(out)
}
