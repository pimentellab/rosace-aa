#' @import dplyr
#' @import ggplot2
#' @importFrom stats optim 
#' @importFrom tidyr pivot_wider unite separate
#' @importFrom MASS fitdistr
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Plot histogram of scores in Rosette object
#'
#' Plot histograms with variant group label and/or mutant group label
#'
#' @param object An Rosette object
#' @param test.group If TRUE, plot histograms with test group label
#'
#' @return A ggplot object
#'
#' @export
#'
PlotScoreHist <- function(object, test.group = FALSE) {

  df <- object@score.df

  if (test.group) {
    p <- ggplot(df, aes(.data$score, fill = .data$test)) +
      geom_histogram(color = "grey", bins = 30,
                     alpha=0.5, position="identity") +
      theme_bw() +
      facet_wrap(vars(.data$test))
  } else {
    p <- ggplot(df, aes(.data$score)) +
      geom_histogram(color = "grey", bins = 30) +
      theme_bw()
  }

  return(p)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param ctrl.label Boolean vector indicating control variants
#'
#' @rdname PlotDisp
#' @method PlotDisp AssayGrowth
#' @export
#'
PlotDisp.AssayGrowth <- function(object, t, ctrl.label, ...) {
  CheckDots(...)

  mat <- object@counts[ctrl.label, ]

  n_arr <- seq(10, nrow(mat), length.out = 50)
  if (missing(t)) {
    phi_arr <- sapply(n_arr, function(n) inferPerVarPhi(mat, n))
  } else {
    phi_arr <- sapply(n_arr, function(n) inferPerVarPhi(mat, n, t))
  }

  data <- data.frame(n_arr, phi_arr)
  p <- ggplot(data, aes(x = n_arr, y = phi_arr)) +
    geom_line() +
    xlab("number of negative control variants randomly chosen") +
    ylab("per variant phi")

  return(p)
}

#' @param name Name of AssayGrowth
#' @param ctrl.col Name of column in score.df that contains type of variants
#' @param ctrl.name Label of control variants in type.col
#'
#' @rdname PlotDisp
#' @method PlotDisp Rosace
#' @export
#'
# main function to use
PlotDisp.Rosace <- function(object, t, name, ctrl.col, ctrl.name, ...) {
  CheckDots(...)

  var.data <- ExtractVarAssay(object, name = name, norm = FALSE)
  ctrl.label <- (var.data[[ctrl.col]] == ctrl.name)
  assay <- ExtractAssay(object, name = name)

  if (missing(t)) {
    p <- PlotDisp(object = assay, ctrl.label = ctrl.label)
  } else {
    p <- PlotDisp(object = assay, t = t, ctrl.label = ctrl.label)
  }

  return(p)
}

#' @param ctrl.label Boolean vector indicating control variants
#'
#' @rdname EstimateDisp
#' @method EstimateDisp AssayGrowth
#' @export
EstimateDisp.AssayGrowth <- function(object, ctrl.label, ...) {
  CheckDots(...)

  mat <- object@counts[ctrl.label, ]
  return(inferPerVarPhi(mat))
}

#' @param name Name of AssayGrowth
#' @param ctrl.col Name of column in score.df that contains type of variants
#' @param ctrl.name Label of control variants in ctrl.col
#'
#' @rdname EstimateDisp
#' @method EstimateDisp Rosace
#' @export
#'
EstimateDisp.Rosace <- function(object, name, ctrl.col, ctrl.name, ...) {
  CheckDots(...)

  var.data <- ExtractVarAssay(object, name = name, norm = FALSE)
  ctrl.label <- var.data[[ctrl.col]] == ctrl.name
  assay <- ExtractAssay(object, name = name)

  if (!isa(assay, "AssayGrowth")) {
    stop("EstimateDisp only supports AssayGrowth object.")
  }

  return(EstimateDisp(object = assay, ctrl.label = ctrl.label))
}

#' @rdname EstimateDispStart
#' @method EstimateDispStart AssayGrowth
#' @export
EstimateDispStart.AssayGrowth <- function(object, ...) {
  CheckDots(...)

  count0 <- object@counts[1, ]
  return(inferDispLibrary(count0))
}

#' @param name Name of AssayGrowth
#' @rdname EstimateDispStart
#' @method EstimateDispStart Rosace
#' @export
#'
EstimateDispStart.Rosace <- function(object, name, ...) {
  CheckDots(...)

  assay <- ExtractAssay(object, name = name)

  if (!isa(assay, "AssayGrowth")) {
    stop("EstimateDisp only supports AssayGrowth object.")
  }

  return(EstimateDispStart(object = assay))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Log likelihood function for dirichlet multinomial
dirichlet_multinomial_log_pdf <- function(x, Phi, p) {
  alpha <- Phi * p
  n <- sum(x)
  alpha_0 <- sum(alpha)
  res <- log(n)
  res <- res + lbeta(alpha_0, n)

  for (i in 1:length(x)) {
    x_k <- x[i]
    alpha_k <- alpha[i]
    if (x_k > 0) {
      res <- res - log(x_k)
      res <- res - lbeta(alpha_k, x_k)
    }
  }
  return(res)
}

# Optim wrapper for negative dirichlet multinomial log likelihood
dirichlet_multinomial_negative_ll_optim_wrapper <- function(x, p) {
  func <- function(Phi) {

    if (is.null(ncol(x)) || (ncol(x) == 1)) {
      return(-dirichlet_multinomial_log_pdf(x, Phi, p))
    }
    else {
      sum <- 0
      for (i in 1:ncol(x)) {
        sum = sum - dirichlet_multinomial_log_pdf(x[, i], Phi, p)
      }
      return(sum)
    }
  }
  return(func)
}

# Infer per variant phi (dispersion) for counts
# Dirichlet alpha = phi * p
inferPerVarPhi <- function(counts, sample, t) {

  # add pseudo-count
  counts <- counts + 0.5

  if (missing(sample)) {
    n <- nrow(counts)
    p <- counts[, 1] / sum(counts[, 1])
    if (missing(t)) {
      x <- counts[, -1]
    } else {
      if (t < 2) {
        stop("t must be greater than 1.")
      }
      x <- counts[, t]
    }
  } else {
    n <- as.integer(sample)
    indices <- sample(seq(1, nrow(counts)), size = n)
    p <- counts[indices, 1] / sum(counts[indices, 1])
    if (missing(t)) {
      x <- counts[indices, -1]
    } else {
      if (t < 2) {
        stop("t must be greater than 1.")
      }
      x <- counts[indices, t]
    }
  }
  Phi_init <- 2.5 * n
  res <- stats::optim(Phi_init, dirichlet_multinomial_negative_ll_optim_wrapper(x, p), method = "BFGS")
  return(res$par[1]/n)
}

# count0: vector input (starting count)
inferDispLibrary <- function(count0) {
  n <- length(count0)
  Phi_init <- 2.5 * n
  p <- rep(1/n, n)
  res <- stats::optim(Phi_init, dirichlet_multinomial_negative_ll_optim_wrapper(count0, p), method = "BFGS")
  return(res$par[1]/n)
}

# fit variant label distribution
fit_label_dist <- function(df) {

  var.dist <- data.frame(label = c("ctrl", unique(df$test)), mean = 0, sd = 0, count = 0)

  # control variants
  score_ctrl <- df$score[df$ctrl == TRUE]
  fit <- MASS::fitdistr(score_ctrl, "normal")
  var.dist[1, 2] <- fit$estimate[1]
  var.dist[1, 3] <- fit$estimate[2]
  var.dist[1, 4] <- length(score_ctrl)

  # all other variants
  vec_label <- unique(df$test)
  for (i in 1:length(vec_label)) {
    score_label <- df$score[(df$test == vec_label[i]) & (df$ctrl == FALSE)]
    fit <- MASS::fitdistr(score_label, "normal")
    var.dist[i + 1, 2] <- fit$estimate[1]
    var.dist[i + 1, 3] <- fit$estimate[2]
    var.dist[i + 1, 4] <- length(score_label)
  }
  
  return(var.dist)
}
