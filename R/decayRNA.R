
#' Fit Exponential and Decay-of-Decay Models to RNA Decay Data.
#'
#' This function fits two models to RNA decay data: a simple exponential decay model
#' and a more complex "decay-of-decay" model, where the decay rate itself may vary over time.
#' It returns model parameters, estimated and observed expression at time 0, estimated half-lives,
#' AIC values, and model statistics. See decayRNA.docx for formalization.
#'
#' @param time A numeric vector representing time points (e.g., hours after transcription inhibition).
#' @param abundance A numeric vector of RNA abundance measurements corresponding to each time point.
#'
#' @return A named list containing fitted model outputs:
#' \describe{
#'   \item{exponential}{Results from the exponential decay model, including:
#'     \itemize{
#'       \item \code{params}: Named coefficients \code{y0} and \code{k}
#'       \item \code{initial_expression_estimate}: Fitted \code{y0}
#'       \item \code{initial_expression_observed}: Observed abundance at \code{t == 0}
#'       \item \code{half_life}: Estimated RNA half-life (\code{log(2)/k})
#'       \item \code{AIC}: Akaike Information Criterion for model
#'       \item \code{stats}: Model fit statistics from \code{\link[broom]{glance}}
#'     }
#'   }
#'   \item{decay_of_decay}{Results from the decay-of-decay model, including:
#'     \itemize{
#'       \item \code{params}: Named coefficients \code{y0}, \code{k1}, and \code{k2}
#'       \item \code{initial_expression_estimate}: Fitted \code{y0}
#'       \item \code{initial_expression_observed}: Observed abundance at \code{t == 0}
#'       \item \code{half_life}: Estimated RNA half-life (numerical estimate)
#'       \item \code{AIC}: Akaike Information Criterion for model
#'       \item \code{stats}: Model fit statistics from \code{\link[broom]{glance}}
#'     }
#'   }
#'   \item{model_comparison}{(If both models converge) A list including:
#'     \itemize{
#'       \item \code{anova}: Model comparison using ANOVA
#'       \item \code{AIC_exp}: AIC of the exponential model
#'       \item \code{AIC_decay_decay}: AIC of the decay-of-decay model
#'       \item \code{best_model}: Name of the model with lower AIC
#'     }
#'   }
#'   \item{observed_expression_at_t0}{The observed value(s) of abundance at time = 0}
#' }
#'
#' @importFrom broom glance
#' @importFrom dplyr %>%
#' @importFrom stats nls coef AIC anova uniroot
#' @examples
#' time <- c(0, 1, 2, 3, 4, 5)
#' abundance <- c(100, 60, 38, 25, 18, 12)
#' fit <- fit_decay_models(time, abundance)
#' fit$exponential$initial_expression_estimate
#' fit$exponential$initial_expression_observed
#' fit$observed_expression_at_t0
#'
#' @export
fit_decay2 <- function(df, expression_col="cpm", time_col="time") {
  library(broom)
  library(dplyr)

  # Check if the required columns exist in the data
  if (!expression_col %in% colnames(df)) {
    stop(paste("Error: ", expression_col, "not in the data."))
  }
  if (!time_col %in% colnames(df)) {
    stop(paste("Error: ", time_col, "not in the data."))
  }


  df$time_num <- as.numeric(df[[time_col]])
  df <- df[order(df$time_num), ]
  df$expraw <- df[[expression_col]]

  ##
  cols <- c("ntime", "nozero", "observed_t0","changed_t0",  "k", "y0_m1", "HF_m1", "AIC_m1",
            "k1", "k2", "y0_m2", "HF_m2", "AIC_m2", "bestModel", "anovP", "m1m2Sig")
  #
  result <- as.data.frame(setNames(as.list(rep(NA, length(cols))), cols),
                          stringsAsFactors = FALSE)

  ## Data issues feeding in single time
  if (nrow(df) <= 1 || is.na(df$expraw[1]) || df$expraw[1] == 0) {
    return(result)
  }

  # Check decrease in t0 to t1 to handle experiment effect
  # Check the first pair has a drop.
  # After this other irregular is not meaningful still
  while (
    nrow(df) >= 2 &&
    !is.na(df$expraw[2]) &&
    (df$expraw[1] < df$expraw[2])) {
    df <- df[-1, ]
    }


  # Done here to account for changed
  result$observed_t0 <- df$expraw[[1]]

  # Check time and zero
  result$ntime <- length(unique(df$time_num[!is.na(df$time_num)]))
  result$nozero <- sum(df$expraw > 0, na.rm = TRUE)

  # Early return
  if (result$ntime <= 2 || result$nozero <= 2) {
    return(result)
  }

  # If not continue
  df$normExp <- df$expraw / as.numeric(result$observed_t0)


  ## Note Normalized exp - fix y0 to 1, - modelling easier
  # Model 1: Exponential Decay
  exp_model <- try(nls(normExp ~ y0 * exp(-k * time_num),
                       data = df,
                       start = list(y0 = 1, k = 0.1)), # y0 = 1 due to norm
                   silent = TRUE)

  # Model 2: Decay of Decay
  decay_decay_model <- try(nls(normExp ~ y0 * exp(-k1 * time_num - k2 * time_num^2),
                               data = df,
                               start = list(y0 = 1, k1 = 0.1, k2 = 0.01)), # y0 = 1 due to norm
                           silent = TRUE)

  if (!inherits(exp_model, "try-error")) {
    exp_params <- coef(exp_model)
    result$k <- exp_params["k"]
    result$y0_m1 <- exp_params["y0"]
    result$HF_m1 <- log(2) / result$k

    result$AIC_m1 <- AIC(exp_model)
    # exp_stats <- glance(exp_model) # Activate for checks
  }

  if (!inherits(decay_decay_model, "try-error")) {
    dd_params <- coef(decay_decay_model)
    result$k1 <- dd_params["k1"]
    result$k2 <- dd_params["k2"]
    result$y0_m2 <- dd_params["y0"]

    t_hf_est <- try(uniroot(function(t) exp(-result$k1*t - result$k2*t^2) - 0.5,
                                   c(0, max(time))), silent = TRUE)
    result$HF_m2 <- if (!inherits(t_hf_est, "try-error")) t_hf_est$root else NA
    result$AIC_m2 <- AIC(decay_decay_model)

    # dd_stats <- glance(decay_decay_model) # Activate for checks
  }

  if (!inherits(exp_model, "try-error") && !inherits(decay_decay_model, "try-error")) {
    ano_res <- anova(exp_model, decay_decay_model)
    result$bestModel <- ifelse(result$AIC_m1 < result$AIC_m2, "Exponential", "Decay_of_decay")
    result$anovP <- ano_res$`Pr(>F)`[2]
    result$m1m2Sig <- !is.na(result$anovP) && result$anovP < 0.05

  }

  return(result)
}




#' Fit exponential decay model to RNA expression data
#'
#' @description
#' Models RNA decay kinetics following transcriptional arrest (e.g., Actinomycin D treatment)
#' using a log-linear regression approach. Assumes that mRNA degradation follows first-order
#' kinetics, such that the rate of change in RNA concentration over time is proportional to its abundance.
#'
#' The model estimates the decay constant (k_decay) by fitting:
#'     log(C_t + 1) = -k_decay * t + log(C_0 + 1)
#' where C_t is RNA abundance at time t, and log1p is used for stability with low values.
#'
#' Additionally, estimates the persistence of an RNA species over time following transcriptional shutoff (e.g., actinomycin D treatment)
#' using the area under the curve (AUC) of expression versus time. The AUC is computed using the trapezoid rule by default.
#'
#' @param df A data frame containing RNA expression and timepoint information.
#' @param expression_col Name of the column in `df` that contains RNA expression values (e.g., "cpm").
#' @param time_col Name of the column in `df` that contains time information (e.g., "time").
#' @param method AUC calculation method: "trapezoid" (default) or "sum"
#' @param type Method for fitting the decay model nls (default) or "linear
#'
#' @return A named vector with:
#'   - half_life: Estimated half-life of the RNA (in the same units as time_col).
#'   - decay_rate: Slope of the linear fit, corresponding to -k_decay.
#'     If expression is not detectable at initial time, returns "not expressed".
#'
#' @examples
#' fit_decay(df, expression_col = "cpm", time_col = "time")
fit_decay <- function(df, expression_col="cpm", time_col="time", method = "trapezoid", type="nls") {
  # Check if the required columns exist in the data
  if (!expression_col %in% colnames(df)) {
    stop(paste("Error: ", expression_col, "not in the data."))
  }
  if (!time_col %in% colnames(df)) {
    stop(paste("Error: ", time_col, "not in the data."))
  }

  # time numeric
  df$time_num <- as.numeric(df[[time_col]])
  ftime <- min(df$time_num, na.rm = TRUE)
  fexp <- df[df$time_num == ftime, expression_col][[1]]
  df$expraw <- df[[expression_col]]
  df[[expression_col]] <- df[[expression_col]] / as.numeric(fexp)


  # Check if there are any missing or invalid time values
  if (any(is.na(df$time_num))) {
    stop("Error: Time column contains missing or invalid values.")
  }

  # Arrange by time
  df <- df[order(df$time_num), ]


  ## Do AUC first
  agg_df <- aggregate(df[[expression_col]],
                      by = list(time_num = df$time_num),
                      FUN = mean, na.rm = TRUE) |>
    setNames(c("time_num", "mean_exp")) |>
    (\(x) x[order(x$time_num), ])()


  x <- agg_df$time_num
  y <- agg_df$mean_exp

  # Check validity
  if (length(x) < 2 || all(is.na(y))) {
    auc = NA
  }

  # AUC calculation
  if (method == "trapezoid") {
    auc <- sum(diff(agg_df$time_num) * (head(agg_df$mean_exp, -1) +
                                          tail(agg_df$mean_exp, -1)) / 2, na.rm = TRUE)
  } else {
    auc <- sum(agg_df$mean_exp, na.rm = TRUE)
  }


  # Require at least one positive expression value at time 0
  if (all(fexp <= 0, na.rm = TRUE)) {
    return(c("half_life" = NA, "decay_rate" = 0, "auc"=auc, "t0.exp"=fexp, "n.time"=nrow(df)))
  }



  # Check for sufficient data points and valid expression values
  if (nrow(df) < 3 || any(df[[expression_col]] < 0, na.rm = TRUE)) {
    return(c("half_life" = NA, "decay_rate" = NA, "auc"=auc, "t0.exp"=fexp, "n.time"=nrow(df)))
  }

  if(type == "nls") {
    tryCatch({
      # Fit the non-linear exponential decay model
      fit <- nls(expraw ~ a * exp(-b * time_num),
                 data = df,
                 start = list(a = 1, b = 0.1),
                 algorithm = "port",
                 lower = c(a = 0, b = 0))

      # Extract decay rate and initial expression
      decay_rate <- coef(fit)["b"]
      t0.exp <- as.numeric(coef(fit)["a"])

      # Compute half-life
      half_life <- log(2) / decay_rate

      # Predict fitted values for AUC
      df$fitted <- predict(fit)
      auc <- sum(diff(df$time_num) * (head(df$fitted, -1) + tail(df$fitted, -1)) / 2)

      # Return results
      return(setNames(c(half_life, decay_rate, auc, t0.exp, nrow(df)),
                      c("half_life", "decay_rate", "auc", "t0.exp", "n.time")))

    }, error = function(e) {
      # Return NAs if the model fails
      return(setNames(c(NA, NA, NA, NA, nrow(df)),
                      c("half_life", "decay_rate", "auc", "t0.exp", "n.time")))
    })

  } else {
    # Fit the decay model using natural logarithm (log base e)
    tryCatch({
      fit <- lm(log1p(df[[expression_col]]) ~ df$time_num)
      # Extract the decay rate (negative of the coefficient)
      decay_rate <- as.numeric(-coef(fit)[2])

      if (decay_rate < 0) {
        # negative rate to NA
        return(c("half_life" = NA, "decay_rate" = decay_rate, "auc"=auc, "t0.exp"=fexp, "n.time"=nrow(df)))
      }

      # This still uses the natural log because of exponential decay model
      half_life <- log(2) / decay_rate
      return(c("half_life" = half_life, "decay_rate" = decay_rate, "auc"=auc, "t0.exp"=fexp, "n.time"=nrow(df)))
    }, error = function(e) {
      return(c("half_life" = NA, "decay_rate" = NA, "auc"=auc, "t0.exp"=fexp, "n.time"=nrow(df)))
    })
  }

}
