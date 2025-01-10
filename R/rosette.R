#' @import dplyr
#' @importFrom readr write_tsv
#' @importFrom tidyr pivot_wider
#' @importFrom methods setClass new
#' @importFrom stats runif quantile
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Rosette Class
#'
#' A Rosette object is a container for inferring summary statistics from a
#' AssayGrowth of a Rosace experiment.
#'
#' @slot score.df A data.frame with columns index, score, pos, wt, mut, ctrl, test, blosum
#' @slot pos.df A data.frame with columns pos, act3, act4, act5
#' @slot var.label label for simulation ("Neutral" + (), ("Neg"), ("Pos"), or ("Neg", "Pos")
#' @slot var.dist A data.frame with columns label, mean, sd, count
#' @slot disp Dispersion of raw count (sequencing step)
#' @slot disp.start Dispersion of starting variant library
#' @slot rounds Number of rounds
#' @slot project.name A character string indicating the name of the project
#'
#' @exportClass Rosette
#'
Rosette <- methods::setClass(
  Class = 'Rosette',
  slots = c(
    score.df = 'data.frame', # index, score, pos, wt, mut, ctrl, test, blosum
    pos.df = 'data.frame', # pos, act3, act4, act5
    var.label = 'vector', #
    var.dist = 'data.frame', # var_group, mean, sd, count
    disp = 'numeric',
    disp.start = "numeric",
    rounds = "numeric",
    project.name = 'character'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param score.name A character string indicating the name of the score
#' @param assay.name If the score is computed from an Assayset, use this assay to compute dispersion
#' @param pos.col A character string indicating the name of the position column in var.data
#' @param ctrl.col A character string indicating the name of the control column in var.data
#' @param ctrl.name A character string indicating the label of the control in the control column
#' @param wt.col For Growth screen, optional, for blosum grouping
#' @param mut.col For Growth screen, optional, for blosum grouping
#' @param aa.code "single" or "triple" amino acid coding, for blosum grouping
#' @param var.label label for simulation ((), ("Neg"), ("Pos"), or ("Neg", "Pos").
#' "Neutral" is included by default.
#'
#' @rdname CreateRosetteObject
#' @method CreateRosetteObject Rosace
#'
#' @export
#'
CreateRosetteObject.Rosace <- function(object, project.name, score.name, assay.name,
                                       pos.col, ctrl.col, ctrl.name, wt.col, mut.col, aa.code = "single",
                                       var.label = c(),
                                       ...) {
  CheckDots(...)

  score <- ExtractScore(object, name = score.name)
  if (!startsWith(score@method, "ROSACE")) {
    stop("The new rosette only supports simulation from Rosace Score object.")
  }

  # estimate dispersion
  if (score@type == "AssayGrowth") {
    assay <- ExtractAssay(object, name = score@assay.name)
    var.data <- ExtractVarAssay(object, name = score@assay.name, norm = FALSE)
    rounds <- score@misc$rounds
  } else if (score@type == "AssaySetGrowth") {
    if (missing(assay.name)) {
      stop("Rosette supports score from type 'AssaySetGrowth' only with given assay.name for dispersion computation.")
    }
    assay <- ExtractAssay(object, name = assay.name)
    var.data <- ExtractVarAssay(object, name = assay.name, norm = FALSE)
    rounds <- assay@rounds
  } else {
    stop("Rosette currently only support score from type 'AssayGrowth' or 'AssaySetGrowth'.")
  }
  ctrl.label <- var.data[[ctrl.col]] == ctrl.name
  disp <- EstimateDisp(object = assay, ctrl.label = ctrl.label)
  disp.start <- EstimateDispStart(object = assay)

  # extract variant info
  var.data <- ExtractVarScore(object, name = score.name)

  # extract hypothesis testing label
  df.score <- OutputScore(object = score, sig.test = 0.05, pos.info = FALSE, blosum.info = FALSE, pos.act.info = FALSE)
  test.label <- df.score$label
  var.label <- sort(var.label)
  if (length(var.label) == 0) {
    warning("Only simulate neutral label! Please check and proceed.")
    test.label <- "Neutral"
  } else if (length(var.label) == 2) {
    if (var.label[1] != "Neg" || var.label[2] != "Pos") {
      stop("Only accepts label 'Neg' and 'Pos'.")
    }
  } else if (length(var.label) == 1 && var.label == "Pos") {
    test.label <- ifelse(test.label == "Neg", "Neutral", test.label)
  } else if (length(var.label) == 1 && var.label == "Neg") {
    test.label <- ifelse(test.label == "Pos", "Neutral", test.label)
  } else {
    stop("Invalid var.label. Accepts {}, {'Neg'}, {'Pos'}, or {'Neg', 'Pos'}. ")
  }
  var.label <- c(var.label, "Neutral")

  # call score-level function to create rosette object
  rosette <- CreateRosetteObject(object = score,
                                 project.name = project.name,
                                 pos.label = var.data[[pos.col]],
                                 ctrl.label = var.data[[ctrl.col]] == ctrl.name,
                                 wt.label = var.data[[wt.col]],
                                 mut.label = var.data[[mut.col]],
                                 test.label = test.label,
                                 var.label = var.label,
                                 aa.code = aa.code,
                                 disp = disp,
                                 disp.start = disp.start,
                                 rounds = rounds)

  return(rosette)
}

#' @param pos.label A vector of position labels
#' @param ctrl.label A vector of control labels (boolean: ctrl or not)
#' @param wt.label A vector of wildtype labels
#' @param mut.label A vector of mutant labels
#' @param test.label A vector of hypothesis testing label (Neg, Neutral, Pos)
#' @param aa.code "single" or "triple" amino acid coding, for blosum grouping
#' @param var.label label for simulation ("Neutral" + (), ("Neg"), ("Pos"), or ("Neg", "Pos")
#' @param disp Dispersion of raw count
#' @param disp.start Dispersion of variant library
#' @param rounds Number of rounds in the orignal assay/assayset
#'
#' @rdname CreateRosetteObject
#' @method CreateRosetteObject Score
#' @export
#'
CreateRosetteObject.Score <- function(object, project.name,
                                      pos.label, ctrl.label, wt.label, mut.label, test.label, aa.code = "single",
                                      var.label,
                                      disp, disp.start,
                                      rounds, ...) {
  CheckDots(...)

  ### score.df vanilla ###
  # score, pos, wt, mut, ctrl, test, blosum
  df <- data.frame(score = object@score[[2]], pos = pos.label, wt = wt.label,
                   mut = mut.label, ctrl = ctrl.label, test = test.label)
  # map blosum group
  df$blosum <- MapBlosumScore(wt.vec = wt.label, mut.vec = mut.label, aa.code = aa.code)

  ### var.dist ###
  # label, mean, sd, count
  var.dist <- fit_label_dist(df = df)
  prop.sig <- sum((df$test != "Neutral") & (df$ctrl == FALSE)) / sum(df$ctrl == FALSE)

  ### pos.df ###
  # pos, act3, act4, act5
  df_pos <- df %>%
    dplyr::group_by(.data$pos, .data$wt) %>%
    dplyr::summarise(mean = mean(.data$score, na.rm = FALSE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(pos = as.numeric(.data$pos)) %>%
    dplyr::arrange(.data$pos)
  # df_pos$act5 <- (abs(df_pos$mean) >= stats::quantile(abs(df_pos$mean), 1 - prop.sig))

  ### blosum.df ###
  df_blosum <- df %>%
    dplyr::select(.data$pos, .data$wt, .data$mut) %>%
    dplyr::mutate(pos = as.numeric(.data$pos)) %>%
    tidyr::pivot_wider(names_from = .data$mut, values_from = .data$wt) %>%
    dplyr::arrange(.data$pos)
  # impute with row value
  for (j in 2:ncol(df_blosum)) {
    df_blosum[[j]] <- MapBlosumScore(wt.vec = df_pos$wt,
                                     mut.vec = rep(colnames(df_blosum)[j], nrow(df_blosum)),
                                     aa.code = aa.code)
  }

  ### score.df ###
  var_info <- df_blosum %>%
    tidyr::pivot_longer(!.data$pos, names_to = "mut", values_to = "blosum") %>%
    dplyr::left_join(df) %>%
    dplyr::select(-.data$wt) %>%
    dplyr::left_join(df_pos[, 1:2])
  # impute missing info
  var_info <- var_info %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ctrl = ifelse(is.na(.data$ctrl) && (.data$mut == .data$wt), TRUE, .data$ctrl),
                  ctrl = ifelse(is.na(.data$ctrl) && (.data$mut != .data$wt), FALSE, .data$ctrl)) %>%
    dplyr::ungroup()
  var_info <- var_info %>%
    tidyr::unite(col ="index", .data$wt, .data$pos, .data$mut, remove = FALSE) %>%
    select(.data$index, .data$score, .data$pos, .data$wt, .data$mut,
          .data$ctrl, .data$test, .data$blosum) %>%
    arrange(.data$pos, .data$mut)

  ### create object ###
  rosette <- methods::new(
    Class = 'Rosette',
    score.df = var_info,
    pos.df = df_pos,
    rounds = rounds,
    project.name = project.name,
    disp = disp,
    disp.start = disp.start,
    var.label = var.label,
    var.dist = var.dist
  )

  return(rosette)
}

