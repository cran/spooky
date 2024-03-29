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
#@param bounds list of numeric vector, with minimum and maximum bounds for each numeric features. Default: NULL
#' @param seed Positive integer. Random seed. Default: 42.
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @return This function returns a list including:
#' \itemize{
#' \item exploration: list of all not-null models, complete with predictions, test metrics, prediction stats and plot
#' \item history: a table with the sampled models, hyper-parameters, validation errors
#' \item best_model: results for the best selected model according to the weighted average rank, including:
#' \itemize{
#' \item testing_errors: testing errors for each time feature for the best selected model (for continuous variables: me, mae, mse, rmsse, mpe, mape, rmae, rrmse, rame, mase, smse, sce, gmrae; for factor variables: czekanowski, tanimoto, cosine, hassebrook, jaccard, dice, canberra, gower, lorentzian, clark)
#' \item preds: for continuous variables, min, max, q25, q50, q75, quantiles at selected ci, mean, sd, mode, skewness, kurtosis, IQR to range, risk ratio, upside probability and divergence for each point fo predicted sequences; for factor variables, min, max, q25, q50, q75, quantiles at selected ci, proportions, difformity (deviation of proportions normalized over the maximum possible deviation), entropy, upgrade probability and divergence for each point fo predicted sequences
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
#' @importFrom stats weighted.mean ecdf
#' @importFrom imputeTS na_kalman
#' @importFrom fANCOVA loess.as
#' @importFrom modeest mlv1
#' @import ggplot2
#' @importFrom moments skewness kurtosis
#' @importFrom stats quantile sd lm rnorm pnorm fft runif na.omit
#' @importFrom scales number
#' @importFrom utils tail head
#' @import greybox
#' @importFrom philentropy distance
#' @importFrom entropy entropy
#' @import fastDummies


#'@examples
#'spooky(time_features, seq_len = c(10, 30), lno = c(1, 30), n_samp = 1)
#'
#'

spooky <- function(df, seq_len = NULL, lno = NULL, n_samp = 30, n_windows = 3, ci = 0.8, smoother = FALSE, dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42)
{

  tic.clearlog()
  tic("time")

  set.seed(seed)

  n_length <- nrow(df)

  class_index <- any(map_lgl(df, ~ is.factor(.x) | is.character(.x)))
  all_classes <- all(class_index)
  numeric_index <- map_lgl(df, ~ is.integer(.x) | is.numeric(.x))
  all_numerics <- all(numeric_index)
  if(!(all_classes | all_numerics)){stop("only all numerics or all classes, not both")}

  if(all_classes){df <- dummy_cols(df, select_columns = NULL, remove_first_dummy = FALSE, remove_most_frequent_dummy = TRUE, ignore_na = FALSE, split = NULL, remove_selected_columns = TRUE); binary_class <- rep(T, ncol(df))}
  if(all_numerics){binary_class <- rep(F, ncol(df))}

  if(anyNA(df) & all_numerics){df <- as.data.frame(na_kalman(df)); message("kalman imputation on time features\n")}
  if(anyNA(df) & all_classes){df <- floor(as.data.frame(na_kalman(df))); message("kalman imputation on time features\n")}
  if(smoother == TRUE & all_numerics){df <- as.data.frame(purrr::map(df, ~ suppressWarnings(loess.as(x=1:n_length, y=.x)$fitted))); message("performing optimal smoothing\n")}

  n_feats <- ncol(df)
  feat_names <- colnames(df)


  if(length(seq_len)==1 & length(lno)==1){n_samp <- 1}
  if(is.null(seq_len)){seq_len <- c(1, round(sqrt(n_length)))}
  if(is.null(lno)){lno <- c(1, round(sqrt(n_length)))}
  sqln_set <- round(runif(n_samp, min(seq_len), max(seq_len)))
  lno_set <- round(runif(n_samp, min(lno), max(lno)))

  exploration <- map2(sqln_set, lno_set, ~ windower(df, seq_len = .x, lno = .y, n_windows, ci, dates, error_scale, error_benchmark, binary_class, seed))
  errors <- round(Reduce(rbind, map(exploration, ~ .x$avg_errors)), 4)
  if(n_samp == 1){errors <- t(as.data.frame(errors))}
  history <- as.data.frame(cbind(seq_len = sqln_set, lno = lno_set, errors))
  rownames(history) <- NULL

  if(n_samp > 1){

    if(all_numerics == TRUE){history <- ranker(history, focus = -c(1, 2), inverse = NULL, absolute = c("me", "mpe", "sce"), reverse = FALSE)}
    if(all_classes == TRUE){history <- ranker(history, focus = -c(1, 2), inverse = NULL, absolute = NULL, reverse = FALSE)}

    best_index <- as.numeric(rownames(history[1,]))
    best_model <- exploration[[best_index]]
  }

  if(n_samp == 1){best_model <- exploration[[1]]}

  if(n_feats > 1){testing_errors <- Reduce(rbind, best_model$errors)}
  if(n_feats == 1){testing_errors <- t(as.data.frame(best_model$errors))}
  rownames(testing_errors) <- feat_names
  preds <- map(best_model$model, ~ .x$quantile_preds)
  names(preds) <- feat_names
  plots <- map(best_model$model, ~ .x$plot)
  names(plots) <- feat_names
  pred_scores <- best_model$pred_scores
  names(pred_scores) <- feat_names

  preds <- map2(preds, pred_scores, ~ cbind(.x, pred_scores = .y))

  best_model <- list(testing_errors = testing_errors, preds = preds, plots = plots)

  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best_model = best_model, time_log = time_log)
  return(outcome)
}


#########
windower <- function(df, seq_len, lno = 1, n_windows = 5, ci = 0.8, dates = NULL, error_scale, error_benchmark, binary_class, seed)
{
  n_length <- nrow(df)
  n_feats <- ncol(df)
  feat_names <- colnames(df)
  too_short <- floor(n_length/(n_windows + 1)) <= seq_len
  if(too_short){stop("not enough data for the validation windows")}
  idx <- c(rep(1, n_length%%(n_windows + 1)), rep(1:(n_windows + 1), each = n_length/(n_windows + 1)))
  models <- sapply(1:n_feats, function(f) map(1:n_windows, ~ engine(df[idx <= .x, f, drop = TRUE], seq_len, lno, testing = head(df[idx == .x + 1, f, drop = TRUE], seq_len), ci, dates, feat_names[f], error_scale, error_benchmark, collected_errors = NULL, binary_class[f], seed)), simplify = FALSE)
  errors <- map(map_depth(models, 2, ~ .x$errors), ~ round(colMeans(Reduce(rbind, .x)), 4))
  pred_scores <- map(map_depth(models, 2, ~ .x$pred_scores), ~ round(colMeans(Reduce(rbind, .x)), 4))
  collected_errors <- map(map_depth(models, 2, ~ .x$raw_errors), ~ as.data.frame(Reduce(rbind, .x)))
  if(n_feats > 1){avg_errors <- apply(Reduce(rbind, errors), 2, mean)}
  if(n_feats == 1){avg_errors <- unlist(errors)}

  model <- sapply(1:n_feats, function(f) engine(df[, f, drop = TRUE], seq_len, lno, testing = NULL, ci, dates, feat_names[f], error_scale, error_benchmark, collected_errors = collected_errors[[f]], binary_class[f], seed), simplify = FALSE)
  outcome <- list(errors = errors, avg_errors = avg_errors, pred_scores = pred_scores, model = model)
  return(outcome)
}

#########
engine <- function(ts, seq_len, lno = 1, testing = NULL, ci = 0.8, dates = NULL, feat_name, error_scale, error_benchmark, collected_errors, binary_class, seed)
{
  orig <- ts

  diff_model <- recursive_diff(ts, best_deriv(ts))
  ts <- diff_model$vector

  n_length <- length(ts)
  transformed <- spectral_transformation(ts)
  jacknife_index <- function(v, n) {ifelse(n <= 0, n <- 1, n); n <- n - 1; map(1:(length(v) - n), ~ v[- (.x:(.x+n))])}
  idx_samples <- jacknife_index(1:length(ts), lno)
  n_trials <- length(idx_samples)
  cmplx <- transformed$complex
  preds <- Reduce(rbind, map(idx_samples, ~ tail(spectral_inversion(c(cmplx[.x], vector("complex", seq_len))), seq_len)))
  preds <- t(apply(preds, 1, function(x) tail(invdiff(x, diff_model$tail_value), seq_len)))

  errors <- NULL
  raw_errors <- NULL
  pred_scores <- NULL
  plot <- NULL
  quantile_preds <- NULL
  if(!is.null(testing)){
    errors <- custom_metrics(testing, colMeans(preds), ts, error_scale, error_benchmark, binary_class)
    raw_errors <- t(replicate(nrow(preds), testing)) - preds
    pred_scores <- prediction_score(preds, testing)
  }

  if(!is.null(collected_errors)){

    raw_pred <- prediction_integration(seeds = colMeans(preds), raw_errors = collected_errors)
    q_model <- qpred(raw_pred, holdout_truth = NULL, ts = orig, ci, error_scale, error_benchmark, binary_class, dates, seed)
    quantile_preds <- q_model$quant_pred

    plot <- plotter(quantile_preds, ci, ts = orig, dates, feat_name)
  }

  outcome <- list(errors = errors, raw_errors = raw_errors, pred_scores = pred_scores, quantile_preds = quantile_preds, plot = plot)
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
  out <- Re(fft(complex, inverse = TRUE)/length(complex))
  return(out)
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


####
jacknife_index <- function(v, n) {ifelse(n <= 0, n <- 1, n); n <- n - 1; map(1:(length(v) - n), ~ v[- (.x:(.x+n))])}


  ###
  recursive_diff <- function(vector, deriv)
  {
    vector <- unlist(vector)
    head_value <- vector("numeric", deriv)
    tail_value <- vector("numeric", deriv)
    if(deriv==0){head_value = NULL; tail_value = NULL}
    if(deriv > 0){for(i in 1:deriv){head_value[i] <- head(vector, 1); tail_value[i] <- tail(vector, 1); vector <- diff(vector)}}
    outcome <- list(vector = vector, head_value = head_value, tail_value = tail_value)
    return(outcome)
  }

  ###
  invdiff <- function(vector, heads, add = FALSE)
  {
    vector <- unlist(vector)
    if(is.null(heads)){return(vector)}
    for(d in length(heads):1){vector <- cumsum(c(heads[d], vector))}
    if(add == FALSE){return(vector[-c(1:length(heads))])} else {return(vector)}
  }


  ###
  prediction_score <- function(integrated_preds, ground_truth)
  {
    pfuns <- apply(integrated_preds, 2, ecdf)
    pvalues <- map2_dbl(pfuns, ground_truth, ~ .x(.y))
    scores <- 1 - 2 * abs(pvalues - 0.5)
    return(scores)
  }


  ###
  ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                       forcat_band = "seagreen2", forcat_line = "seagreen4", hist_line = "gray43", label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
  {
    if(is.character(y_hist)){y_hist <- as.factor(y_hist)}
    if(is.character(y_forcat)){y_forcat <- factor(y_forcat, levels = levels(y_hist))}
    if(is.character(lower)){lower <- factor(lower, levels = levels(y_hist))}
    if(is.character(upper)){upper <- factor(upper, levels = levels(y_hist))}

    n_class <- NULL
    if(is.factor(y_hist)){class_levels <- levels(y_hist); n_class <- length(class_levels)}

    all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = as.numeric(c(y_hist, y_forcat)))
    forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = as.numeric(y_forcat))

    if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- as.numeric(lower); forcat_data$upper <- as.numeric(upper)}

    plot <- ggplot()+ geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
    if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes_string(x = "x_forcat", ymin = "lower", ymax = "upper"), alpha = 0.3, fill = forcat_band)}
    plot <- plot + geom_line(data = forcat_data, aes_string(x = "x_forcat", y = "y_forcat"), color = forcat_line, size = line_size)
    if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
    if(is.null(dbreak)){plot <- plot + xlab(label_x)}
    if(is.null(n_class)){plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)}
    if(is.numeric(n_class)){plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), breaks = 1:n_class, labels = class_levels)}
    plot <- plot + ylab(label_y)  + theme_bw()
    plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

    return(plot)
  }

  ###
  qpred <- function(raw_pred, holdout_truth = NULL, ts, ci, error_scale = "naive", error_benchmark = "naive", binary_class = F, dates, seed = 42)
  {
    set.seed(seed)

    raw_pred <- doxa_filter(ts, raw_pred, binary_class)
    quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))

    if(binary_class == F)
    {
      p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), kurtosis = suppressWarnings(kurtosis(x[is.finite(x)], na.rm = TRUE)), skewness = suppressWarnings(skewness(x[is.finite(x)], na.rm = TRUE)))}
      quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
      p_value <- apply(raw_pred, 2, function(x) ecdf(x)(seq(min(raw_pred), max(raw_pred), length.out = 1000)))
      divergence <- c(max(p_value[,1] - seq(0, 1, length.out = 1000)), apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
      upside_prob <- c(mean((raw_pred[,1]/tail(ts, 1)) > 1, na.rm = T), apply(apply(raw_pred[,-1, drop = FALSE]/raw_pred[,-ncol(raw_pred), drop = FALSE], 2, function(x) x > 1), 2, mean, na.rm = T))
      iqr_to_range <- (quant_pred[, "75%"] - quant_pred[, "25%"])/(quant_pred[, "max"] - quant_pred[, "min"])
      above_to_below_range <- (quant_pred[, "max"] - quant_pred[, "50%"])/(quant_pred[, "50%"] - quant_pred[, "min"])
      quant_pred <- round(cbind(quant_pred, iqr_to_range, above_to_below_range, upside_prob, divergence), 4)
    }

    if(binary_class == T)
    {
      p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), prop = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), entropy = entropy(x))}
      quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
      p_value <- apply(raw_pred, 2, function(x) ecdf(x)(c(0, 1)))
      divergence <- c(max(p_value[,1] - c(0, 1)), apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
      upgrade_prob <- c(mean(((raw_pred[,1] + 1)/tail(ts + 1, 1)) > 1, na.rm = T), apply(apply((raw_pred[,-1, drop = FALSE] + 1)/(raw_pred[,-ncol(raw_pred), drop = FALSE] + 1), 2, function(x) x > 1), 2, mean, na.rm = T))
      quant_pred <- round(cbind(quant_pred, upgrade_prob = upgrade_prob, divergence = divergence), 4)
    }

    testing_error <- NULL
    if(!is.null(holdout_truth))
    {
      mean_pred <- colMeans(raw_pred)
      testing_error <- custom_metrics(holdout_truth, mean_pred, ts, error_scale, error_benchmark, binary_class)
      pred_scores <- round(prediction_score(raw_pred, holdout_truth), 4)
      quant_pred <- cbind(quant_pred, pred_scores = pred_scores)
    }

    if(is.Date(dates))
    {
      new_dates<- seq.Date(tail(dates, 1), tail(dates, 1) + nrow(quant_pred) * mean(diff(dates)), length.out = nrow(quant_pred))
      rownames(quant_pred) <- as.character(new_dates)
    }
    else
    {
      rownames(quant_pred) <- paste0("t", 1:nrow(quant_pred))
    }

    outcome <- list(quant_pred = quant_pred, testing_error = testing_error)
    return(outcome)
  }

  ###
  plotter <- function(quant_pred, ci, ts, dates = NULL, feat_name)
  {
    seq_len <- nrow(quant_pred)
    n_ts <- length(ts)

    if(is.Date(dates))
    {
      new_dates<- seq.Date(tail(dates, 1), tail(dates, 1) + seq_len * mean(diff(dates)), length.out = seq_len)
      x_hist <- dates
      x_forcat <- new_dates
      rownames(quant_pred) <- as.character(new_dates)
    }
    else
    {
      x_hist <- 1:n_ts
      x_forcat <- (n_ts + 1):(n_ts + seq_len)
      rownames(quant_pred) <- paste0("t", 1:seq_len)
    }

    quant_pred <- as.data.frame(quant_pred)
    x_lab <- paste0("Forecasting Horizon for sequence n = ", seq_len)
    y_lab <- paste0("Forecasting Values for ", feat_name)

    ###if(is.numeric(n_class) & !is.null(level_names)){ts <- level_names[ts + 1]}
    lower_b <- paste0((1-ci)/2 * 100, "%")
    upper_b <- paste0((ci+(1-ci)/2) * 100, "%")

    plot <- ts_graph(x_hist = x_hist, y_hist = ts, x_forcat = x_forcat, y_forcat = quant_pred[, "50%"], lower = quant_pred[, lower_b], upper = quant_pred[, upper_b], label_x = x_lab, label_y = y_lab)
    return(plot)
  }

  ###
  doxa_filter <- function(ts, mat, binary_class = F)
  {
    discrete_check <- all(ts%%1 == 0)
    all_positive_check <- all(ts >= 0)
    all_negative_check <- all(ts <= 0)
    monotonic_increase_check <- all(diff(ts) >= 0)
    monotonic_decrease_check <- all(diff(ts) <= 0)

    monotonic_fixer <- function(x, mode)
    {
      model <- recursive_diff(x, 1)
      vect <- model$vector
      if(mode == 0){vect[vect < 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
      if(mode == 1){vect[vect > 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
      return(vect)
    }

    if(all_positive_check){mat[mat < 0] <- 0}
    if(all_negative_check){mat[mat > 0] <- 0}
    if(discrete_check){mat <- round(mat)}
    if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}
    if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}

    if(binary_class == T){mat[mat > 1] <- 1; mat[mat < 1] <- 0}
    mat <- na.omit(mat)

    return(mat)
  }


  ###
  custom_metrics <- function(holdout, forecast, actuals, error_scale = "naive", error_benchmark = "naive", binary_class = F)
  {

    if(binary_class == F)
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
    }

    if(binary_class == T)
    {
      dice <- suppressMessages(distance(rbind(holdout, forecast), method = "dice"))
      jaccard <- suppressMessages(distance(rbind(holdout, forecast), method = "jaccard"))
      cosine <- suppressMessages(distance(rbind(holdout, forecast), method = "cosine"))
      canberra <- suppressMessages(distance(rbind(holdout, forecast), method = "canberra"))
      gower <- suppressMessages(distance(rbind(holdout, forecast), method = "gower"))
      tanimoto <- suppressMessages(distance(rbind(holdout, forecast), method = "tanimoto"))
      hassebrook <- 1 - suppressMessages(distance(rbind(holdout, forecast), method = "hassebrook"))
      taneja <- suppressMessages(distance(rbind(holdout, forecast), method = "taneja"))
      lorentzian <- suppressMessages(distance(rbind(holdout, forecast), method = "lorentzian"))
      clark <- suppressMessages(distance(rbind(holdout, forecast), method = "clark"))
      sorensen <- suppressMessages(distance(rbind(holdout, forecast), method = "sorensen"))
      harmonic_mean <- suppressMessages(distance(rbind(holdout, forecast), method = "harmonic_mean"))
      avg <- suppressMessages(distance(rbind(holdout, forecast), method = "avg"))

      out <- round(c(dice, jaccard, cosine, canberra, gower, tanimoto, hassebrook, taneja, lorentzian, clark, sorensen, harmonic_mean, avg), 4)
    }

    return(out)
  }

  ###
  prediction_integration <- function(seeds, raw_errors){as.matrix(as.data.frame(map2(seeds, raw_errors, ~ .x + sample(.y, size = 1000, replace = TRUE))))}

  ###
  ranker <- function(df, focus, inverse = NULL, absolute = NULL, reverse = FALSE)
  {
    rank_set <- df[, focus, drop = FALSE]
    if(!is.null(inverse)){rank_set[, inverse] <- - rank_set[, inverse]}###INVERSION BY COL NAMES
    if(!is.null(absolute)){rank_set[, absolute] <- abs(rank_set[, absolute])}###ABS BY COL NAMES
    index <- apply(scale(rank_set), 1, mean, na.rm = TRUE)
    if(reverse == FALSE){df <- df[order(index),]}
    if(reverse == TRUE){df <- df[order(-index),]}
    return(df)
  }

