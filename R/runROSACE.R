#' @import cmdstanr
#' @import dplyr
#' @importFrom readr write_tsv
#' @importFrom stringr str_detect
#' @importFrom stats sd quantile
#' @importFrom posterior default_convergence_measures
#' @importFrom tidyr drop_na
#' @importFrom impute impute.knn
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check Stan Setup
#'
#' @param mc.cores Number of cores to use for parallel builds
#' @param install.update whether to update CmdStan
#'
#' @return None
#'
#' @export
#'
CheckStanSetup <- function(mc.cores, install.update = TRUE) {
  check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
  if (install.update) {
    install_cmdstan(cores = mc.cores, quiet = TRUE)
  }
  cmdstan_path()
  cmdstan_version()
}

#' Compile Stan model
#'
#' @param file Stan model file path
#' @param print whether to print model out model
#'
#' @return CmdStanModel object
#'
#' @export
#'
CompileModel <- function(file, print = TRUE) {
  mod <- cmdstanr::cmdstan_model(file)
  if (print) {
    mod$print()
  }
  mod$exe_file()

  return(mod)
}

#' Run Stan MCMC with input list and compiled model
#'
#' @param input list of input data
#' @param mod compiled model
#' @param seed random seed
#' @param refresh number of iterations between progress updates
#'
#' @return CmdStanMCMC object
#'
#' @export
#'
MCMCRunStan <- function(input, mod, seed = 100, refresh = 100) {
  fit <- mod$sample(
    data = input,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    refresh = refresh,
    save_cmdstan_config=TRUE
  )

  diagnostics <- fit$diagnostic_summary()
  print(diagnostics)

  return(fit)
}

#' Extract MCMC diagnostics from CmdStanMCMC object
#'
#' @param fit CmdStanMCMC object
#' @param sampler whether to get diagnostics summary or
#' extract detailed sampler diagnostics
#'
#' @return list of diagnostics summary or data.frame of sampler diagnostics
#'
#' @export
#'
MCMCDiagnostics <- function(fit, sampler = FALSE) {

  if (!sampler) {
    diags_summary <- fit$diagnostic_summary()
    return(diags_summary)
  } else {
    diags_sampler <- fit$sampler_diagnostics(format = "df")
    return(diags_sampler)
  }

}

#' Extract "functional score" posterior distribution from CmdStanMCMC object
#'
#' @param fit CmdStanMCMC object
#' @param param.key parameter key to extract
#' @param param.post posterior distribution of parameters (optional)
#' If not given, extract from fit object.
#' @param savefile If given, save tsv output to file
#' @param output.lfsr whether to compute lfsr
#'
#' @return data.frame of "functional score" posterior distribution
#'
#' @export
#'
MCMCScoreDf <- function(fit, param.key, param.post, savefile, output.lfsr = TRUE){

  if (missing(param.post)) {
    # warning("No input of 'param.post'. By default calling summary measures
    #         through 'cmdstanr' summary. Might be slow.",
    #         call. = FALSE,
    #         immediate. = TRUE)

    param.name <- fit$metadata()$variables
    param.name <- param.name[stringr::str_detect(param.name,
                                                 paste("^", param.key, "\\[", sep = ""))]
    if (length(param.name) == 0) {
      param.name <- fit$metadata()$variables
      param.name <- param.name[stringr::str_detect(param.name, paste("^", param.key, sep = ""))]
    }
    if (length(param.name) == 0) {
      stop("No parameter found. Check the spelling of param.key.")
    }

    df <-
      fit$summary(param.name, mean, stats::sd,
                  ~stats::quantile(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    colnames(df)[3] <- "sd"

  } else {
    df <- param.post %>%
      dplyr::filter(stringr::str_detect(.data$variable,
                                        paste("^", param.key, "\\[", sep = "")))
    if (nrow(df) == 0) {
      stop("'param.post' does not contain variables with 'param.key'.")
    }
  }

  if (output.lfsr) {
    lfsr <- MCMCLfsr(fit, param.key = param.key)
    df$lfsr <- lfsr
  }

  if (!missing(savefile)) {
    readr::write_tsv(df, file = savefile)
  }

  return(df)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param pos.label vector of true position of variants
#' @param ctrl.label vector of whether variant is control
#' @param stop.label vector of whether variant is stop/nonsense
#' @param blosum.label vector of blosum score of variants
#' @param thred integer, threshold for number of variants per position label (index)
#' Passed to function "varPosIndexMap"
#‘
#' @rdname GenRosaceInput
#' @method GenRosaceInput AssayGrowth
#' @export
GenRosaceInput.AssayGrowth <- function(object, save.input, pos.label, ctrl.label, stop.label, blosum.label = NA, thred = 10, ...) {
  CheckDots(...)

  # generate variants - mean count group mapping
  # TODO: change 25 heuristics here
  row_idx <- apply(array(object@norm.var.names), 1, FUN = function(x) {which(x == object@var.names)})
  raw.counts <- object@counts[row_idx, ]
  vMAPm <- ceiling(rank(rowSums(raw.counts, na.rm = TRUE))/25)

  if (!is.na(pos.label[1]) && !is.na(blosum.label[1])) {

    # generate variants - position mapping
    df_map <-
      varPosIndexMap(var.names = object@norm.var.names,
                     pos.label = pos.label,
                     ctrl.label = ctrl.label,
                     stop.label = stop.label,
                     blosum.label = blosum.label,
                     thred = thred)

    P_syn <- 0
    if (!is.na(ctrl.label[1])) {
      P_syn <- P_syn + length(unique(df_map[df_map$ctrl, ]$index))
    }
    if (!is.na(stop.label[1])) {
      P_syn <- P_syn + length(unique(df_map[df_map$stop, ]$index))
    }

    if (max(df_map$index) != length(unique(df_map$index))) {
      stop("Error when generating position index.")
    }

    blosum_count <- as.numeric(table(df_map$blosum))

    # generate input list
    input <- list(m = object@norm.counts,
                  T = object@rounds + 1,
                  t = seq(0, object@rounds)/object@rounds,
                  V = length(object@norm.var.names),
                  vMAPp = df_map$index,
                  P = length(unique(df_map$index)),
                  vMAPm = vMAPm,
                  vMAPb = df_map$blosum,
                  B = max(df_map$blosum),
                  M = max(vMAPm),
                  blosum_count = blosum_count,
                  P_syn = P_syn)

  } else if (!is.na(pos.label[1])) {
    # generate variants - position mapping
    df_map <-
      varPosIndexMap(var.names = object@norm.var.names,
                     pos.label = pos.label,
                     ctrl.label = ctrl.label,
                     stop.label = stop.label,
                     blosum.label = blosum.label,
                     thred = thred)
    if (max(df_map$index) != length(unique(df_map$index))) {
      stop("Error when generating position index.")
    }

    # generate input list
    input <- list(m = object@norm.counts,
                  T = object@rounds + 1,
                  t = seq(0, object@rounds)/object@rounds,
                  V = length(object@norm.var.names),
                  vMAPp = df_map$index,
                  P = length(unique(df_map$index)),
                  vMAPm = vMAPm,
                  M = max(vMAPm))

  } else {
    # generate variants "map"
    df_map <- data.frame(variants = object@norm.var.names)

    # generate input list
    input <- list(m = object@norm.counts,
                  T = object@rounds + 1,
                  t = seq(0, object@rounds)/object@rounds,
                  V = length(object@norm.var.names),
                  vMAPm = vMAPm,
                  M = max(vMAPm))
  }

  # save and return
  if (!missing(save.input)) {
    save(input, df_map, file = save.input)
  }
  return(list(input = input, df_map = df_map))

}

#' @param pos.label vector of true position of variants
#' @param ctrl.label vector of whether variant is control
#' @param stop.label vector of whether variant is stop/nonsense
#' @param blosum.label vector of blosum score of variants
#' @param thred integer, threshold for number of variants per position label (index)
#' Passed to function "varPosIndexMap"
#' @rdname GenRosaceInput
#' @method GenRosaceInput AssaySetGrowth
#' @export
GenRosaceInput.AssaySetGrowth <- function(object, save.input, pos.label, ctrl.label, stop.label, blosum.label = NA, thred = 10, ...) {
  CheckDots(...)

  # generate variants - mean count group mapping
  # TODO: change 25 heuristics here
  vMAPm <- ceiling(rank(rowSums(object@raw.counts, na.rm = TRUE))/25)

  # impute the combined.counts matrix if not complete
  impute.output <- imputeAssaysCountKNN(object@combined.counts, object@rounds)
  counts <- impute.output$counts
  rounds <- impute.output$rounds

   if (!is.na(pos.label[1]) && !is.na(blosum.label[1])) {

    # generate variants - position mapping
    df_map <-
      varPosIndexMap(var.names = object@var.names,
                     pos.label = pos.label,
                     ctrl.label = ctrl.label,
                     stop.label = stop.label,
                     blosum.label = blosum.label,
                     thred = thred)

    P_syn <- 0
    if (!is.na(ctrl.label[1])) {
      P_syn <- P_syn + length(unique(df_map[df_map$ctrl, ]$index))
    }
    if (!is.na(stop.label[1])) {
      P_syn <- P_syn + length(unique(df_map[df_map$stop, ]$index))
    }

    if (max(df_map$index) != length(unique(df_map$index))) {
      stop("Error when generating position index.")
    }

    blosum_count <- as.numeric(table(df_map$blosum))

    # generate input list
    input <- list(m = counts,
                  T = sum(rounds + 1),
                  t = unlist(lapply(rounds, function(x) seq(0, x)/max(rounds))),
                  V = length(object@var.names),
                  vMAPp = df_map$index,
                  P = length(unique(df_map$index)),
                  vMAPm = vMAPm,
                  vMAPb = df_map$blosum,
                  B = max(df_map$blosum),
                  M = max(vMAPm),
                  blosum_count = blosum_count,
                  P_syn = P_syn)

  } else if (!is.na(pos.label[1])) {
    # generate variants - position mapping
    df_map <-
      varPosIndexMap(var.names = object@var.names,
                    pos.label = pos.label,
                    ctrl.label = ctrl.label,
                    stop.label = stop.label,
                    blosum.label = blosum.label,
                    thred = thred)

    if (max(df_map$index) != length(unique(df_map$index))) {
      stop("Error when generating position index.")
    }

    # generate input list
    input <- list(m = counts,
                  T = sum(rounds + 1),
                  t = unlist(lapply(rounds, function(x) seq(0, x)/max(rounds))),
                  V = length(object@var.names),
                  vMAPp = df_map$index,
                  P = length(unique(df_map$index)),
                  vMAPm = vMAPm,
                  M = max(vMAPm))
  } else {
    df_map <- data.frame(variants = object@var.names)

    # generate input list
    input <- list(m = counts,
                  T = sum(rounds + 1),
                  t = unlist(lapply(rounds, function(x) seq(0, x)/max(rounds))),
                  V = length(object@var.names),
                  vMAPm = vMAPm,
                  M = max(vMAPm))
  }

  if (!missing(save.input)) {
    save(input, df_map, file = save.input)
  }
  return(list(input = input, df_map = df_map))

}

#' @rdname MCMCCreateScore
#' @method MCMCCreateScore Assay
#' @export
MCMCCreateScore.Assay <- function(object, main.score,
                                  param.post, diags) { # optional, can be missing

  # score_all <- cbind(var.map, main.score)
  score_all <- main.score

  score <- score_all %>%
    dplyr::select(.data$variants, .data$mean, .data$sd, .data$lfsr)
  optional.score <- score_all %>%
    dplyr::select(-.data$variants, -.data$mean, -.data$sd, -.data$lfsr)

  misc <- list()
  if (!missing(param.post)) {
    misc <- append(misc, list(param.post = param.post))
  }
  if (!missing(diags)) {
    misc <- append(misc, list(diags = diags))
  }
  if (isa(object, "AssayGrowth")) { # for Rosette Object
    misc <- append(misc, list(rounds = object@rounds))
  }

  score <- CreateScoreObject(method = "ROSACE",
                             type = class(object)[1],
                             assay.name = names(object),
                             score = score,
                             optional.score = optional.score,
                             misc = misc)
}

#' @rdname MCMCCreateScore
#' @method MCMCCreateScore AssaySet
#' @export
MCMCCreateScore.AssaySet <- function(object, main.score,
                                     param.post, diags) { # optional, can be missing

  # score_all <- cbind(var.map, main.score)
  score_all <- main.score

  score <- score_all %>% dplyr::select(.data$variants, .data$mean, .data$sd, .data$lfsr)
  optional.score <- score_all %>% dplyr::select(-.data$variants, -.data$mean, -.data$sd, -.data$lfsr)

  misc <- list()
  if (!missing(param.post)) {
    misc <- append(misc, list(param.post = param.post))
  }
  if (!missing(diags)) {
    misc <- append(misc, list(diags = diags))
  }
  if (isa(object, "AssaySetGrowth")) { # for Rosette Object
    misc <- append(misc, list(rounds = object@rounds))
  }

  score <- CreateScoreObject(method = "ROSACE",
                             type = class(object)[1],
                             assay.name = names(object),
                             score = score,
                             optional.score = optional.score,
                             misc = misc)
}

#' @param pos.label vector of true position of variants
#' @param ctrl.label vector of whether variant is in the control group (NA if none provided)
#' @param stop.label vector of whether variant is in the stop/nonsense group (NA if none provided)
#' @param blosum.label vector of blosum score of variants
#' @param pos.act For Growth screen, optional, boolean, whether to fit changepoint position activation model
#'
#' @rdname RunRosace
#' @method RunRosace AssayGrowth
#' @export
#'
RunRosace.AssayGrowth <- function(object, savedir, mc.cores = 4, debug = FALSE, install = TRUE,
                                  pos.label, ctrl.label, stop.label, blosum.label, pos.act = FALSE, ...) {
  CheckDots(..., args = "thred")
  return(helperRunRosaceGrowth(object = object,
                               savedir = savedir,
                               mc.cores = mc.cores,
                               pos.label = pos.label,
                               ctrl.label = ctrl.label,
                               stop.label = stop.label,
                               blosum.label = blosum.label,
                               pos.act = pos.act,
                               debug = debug,
                               install = install,
                               ...))
}

#' @param pos.label vector of true position of variants
#' @param ctrl.label vector of whether variant is in the control group (NA if none provided)
#' @param stop.label vector of whether variant is in the stop/nonsense group (NA if none provided)
#' @param blosum.label vector of blosum score of variants
#' @param pos.act For Growth screen, optional, boolean, whether to fit changepoint position activation model
#'
#' @rdname RunRosace
#' @method RunRosace AssaySetGrowth
#' @export
#'
RunRosace.AssaySetGrowth <- function(object, savedir, mc.cores = 4, debug = FALSE, install = TRUE,
                                     pos.label, ctrl.label, stop.label, blosum.label, pos.act = FALSE, ...) {
  CheckDots(..., args = "thred")
  return(helperRunRosaceGrowth(object = object,
                               savedir = savedir,
                               mc.cores = mc.cores,
                               pos.label = pos.label,
                               ctrl.label = ctrl.label,
                               stop.label = stop.label,
                               blosum.label = blosum.label,
                               pos.act = pos.act,
                               debug = debug,
                               install = install,
                               ...))
}

#' @param name Name of the object to be analyzed
#' @param type "Assay" or "AssaySet"
#' @param pos.col For Growth screen, the column name for position in the var.data
#' (optional in no_pos mode)
#' @param ctrl.col For Growth screen, optional for control to have separate position indexes
#' @param ctrl.name For Growth screen, optional, the name of the control type
#' @param stop.col For Growth screen, optional for stop/nonsense mutations to have separate position indexes
#' @param stop.name For Growth screen, optional, the name of the stop/nonsense type
#' @param wt.col For Growth screen, optional, for blosum grouping
#' @param mut.col For Growth screen, optional, for blosum grouping
#' @param aa.code "single" or "triple" amino acid coding, for blosum grouping
#' @param pos.act For Growth screen, optional, boolean, whether to fit changepoint position activation model
#'
#' @rdname RunRosace
#' @method RunRosace Rosace
#' @export
#'
RunRosace.Rosace <- function(object, savedir, mc.cores = 4, debug = FALSE, install = TRUE,
                             name, type,
                             pos.col, ctrl.col, ctrl.name, stop.col, stop.name,
                             wt.col, mut.col, aa.code, pos.act = FALSE, ...) {

  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }
  
  # Extract Assay
  if (type == "Assay") {
    sub_object <- ExtractAssay(object, name)
  } else if (type == "AssaySet") {
    sub_object <- ExtractAssaySet(object, name)
  } else {
    stop("Unsupported class type. Provide Assay or AssaySet.")
  }

  # Generate Score Object
  if (isa(sub_object, "AssayGrowth") || isa(sub_object, "AssaySetGrowth")) {
    # Growth screen
    # pos.col optional
    if (missing(pos.col)) {
      warning("position column (pos.col) not provided. run the no-position model.")
      pos.label <- NA
    } else {
      if (type == "Assay") {
        pos.label <- ExtractVarAssay(object, name, norm = TRUE)[[pos.col]]
      } else {
        pos.label <- ExtractVarAssaySet(object, name)[[pos.col]]
      }
    }

    # ctrl.col optional
    if (missing(ctrl.col) || missing(ctrl.name)) {
      warnings("control column (ctrl.col) or name (ctrl.name) not provided.")
      ctrl.label <- NA
    } else {
      if (type == "Assay") {
        ctrl.label <- ExtractVarAssay(object, name, norm = TRUE)[[ctrl.col]] == ctrl.name
      } else {
        ctrl.label <- ExtractVarAssaySet(object, name)[[ctrl.col]] == ctrl.name
      }
    }

    # stop.col optional
    if (missing(stop.col) || missing(stop.name)) {
      warnings("stop column or name not provided.")
      stop.label <- NA
    } else {
      if (type == "Assay") {
        stop.label <- ExtractVarAssay(object, name, norm = TRUE)[[stop.col]] == stop.name
      } else {
        stop.label <- ExtractVarAssaySet(object, name)[[stop.col]] == stop.name
      }
    }

    # blosum.col optional
    if (missing(wt.col) || missing(mut.col) || missing(aa.code)) {
      warnings("wildtype aa column (wt.col), mutation aa column (mut.col), or amino acid coding scheme (aa.code) not provided.")
      blosum.label <- NA
      pos.act <- FALSE
    } else {
      if (type == "Assay") {
        wt.label <- ExtractVarAssay(object, name, norm = TRUE)[[wt.col]]
        mut.label <- ExtractVarAssay(object, name, norm = TRUE)[[mut.col]]
        blosum.label <- MapBlosumScore(wt.vec = wt.label, mut.vec = mut.label, aa.code = aa.code)
      } else {
        wt.label <- ExtractVarAssaySet(object, name)[[wt.col]]
        mut.label <- ExtractVarAssaySet(object, name)[[mut.col]]
        blosum.label <- MapBlosumScore(wt.vec = wt.label, mut.vec = mut.label, aa.code = aa.code)
      }
    }

    # run Rosace
    # AssayGrowth/AssaySetGrowth: pos.label (could be NA), ctrl.label (could be NA), stop.label (could be NA), thred (optional)
    score <- RunRosace(object = sub_object,
              savedir = savedir,
              mc.cores = mc.cores,
              debug = debug,
              install = install,
              pos.label = pos.label,
              ctrl.label = ctrl.label,
              stop.label = stop.label,
              blosum.label = blosum.label,
              pos.act = pos.act,
              ...)
  } else {
    stop("THIS IS VERY WRONG. Check ExtractAssay and ExtractAssaySet.")
  }

  if (debug) {
    fit <- score$fit
    score <- score$score
    return(list(fit = fit, score = score))

  } else {
    # Add Score Object to Rosace if not debugging
    object <- AddScoreData(object = object, score = score)
    return(object)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ...: thred (optional)
helperRunRosaceGrowth <- function(object, savedir, mc.cores,
                                  pos.label, ctrl.label, stop.label, blosum.label, pos.act = FALSE,
                                  debug = FALSE, install = TRUE, ...) {
  # create directory if not exists
  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }

  # stan check
  CheckStanSetup(mc.cores = mc.cores, install.update = install)

  # model
  if (!is.na(pos.label[1]) && !is.na(blosum.label[1])) {

    # processing control label
    if (is.na(ctrl.label[1])) {
      ctrl.label <- (blosum.label == 5)  # hard code synonymous variant is 5
    }
    if (sum(ctrl.label) == 0) {
      ctrl.label <- NA
    } else if (sum(ctrl.label) < 5) {
      warning("Number of control variants is less than 5. Might cause divergence at the control phi.")
    }
    # processing stop label
    if (!is.na(stop.label[1])) {
      if (sum(stop.label) == 0) {
        stop.label <- NA
      } else if (sum(stop.label) < 5) {
        warning("Number of nonsense variants is less than 5. Might cause divergence at the nonsense phi.")
      }
    }

    # write model
    if (pos.act) {
      if ((!is.na(ctrl.label[1])) || (!is.na(stop.label[1]))) {
        mod.file <- WriteStanModel(type = "growth_pos_blosum_act")
      } else {
        mod.file <- WriteStanModel(type = "growth_pos_blosum_act_nosyn")
      }
      method <- "ROSACE3"
    } else {
      if ((!is.na(ctrl.label[1])) || (!is.na(stop.label[1]))) {
        mod.file <- WriteStanModel(type = "growth_pos_blosum")
      } else {
        mod.file <- WriteStanModel(type = "growth_pos_blosum_nosyn")
      }
      method <- "ROSACE2"
      warnings("Running blosum model without position activation.")
    }
  } else if (!is.na(pos.label[1])) {
    # write pos-aware Rosace model
    mod.file <- WriteStanModel(type = "growth_pos")
    method <- "ROSACE1"
    warning("Running position model. No blosum label input.")
  } else {
    # write pos-unaware Rosace model
    mod.file <- WriteStanModel(type = "growth_nopos")
    method <- "ROSACE0"
    warning("Running no-position model. Ignore control and blosum label.")
  }
  mod <- CompileModel(file = mod.file, print = FALSE)

  # growth: pos.label, thred in ...
  input.list <- GenRosaceInput(object = object,
                               save.input =
                                 paste(savedir, "/input_", names(object),
                                       ".RData", sep = ""),
                               pos.label = pos.label,
                               ctrl.label = ctrl.label,
                               stop.label = stop.label,
                               blosum.label = blosum.label,
                               ...)
  input <- input.list$input
  df_map <- input.list$df_map

  # MCMC sampling
  # WARNING: FIT is a temporary environment!!!
  fit <- MCMCRunStan(input, mod, seed = 100, refresh = 10)

  # MCMC diagnostics
  diags <- MCMCDiagnostics(fit, sampler = FALSE)

  # MCMC score
  main.score <- MCMCScoreDf(fit, param.key = "beta",
                            savefile = paste(savedir, "/beta.tsv", sep = ""),
                            output.lfsr = TRUE)
  main.score <- cbind(df_map, main.score)

  epsilon2 <-  MCMCScoreDf(fit, param.key = "epsilon2",
                           savefile = paste(savedir, "/epsilon2.tsv", sep = ""),
                           output.lfsr = FALSE)

  if (!is.na(pos.label[1]) && !is.na(blosum.label[1])) {
    nu <- MCMCScoreDf(fit, param.key = "nu",
                      savefile = paste(savedir, "/nu.tsv", sep = ""),
                      output.lfsr = FALSE)

    if (pos.act) {
      rho <- MCMCScoreDf(fit, param.key = "rho",
                         savefile = paste(savedir, "/rho.tsv", sep = ""),
                         output.lfsr = FALSE)
    }
  }

  if (!is.na(pos.label[1])) {
    sigma2 <-  MCMCScoreDf(fit, param.key = "sigma2",
                          savefile = paste(savedir, "/sigma2.tsv", sep = ""),
                          output.lfsr = FALSE)
    phi <-  MCMCScoreDf(fit, param.key = "phi",
                        savefile = paste(savedir, "/phi.tsv", sep = ""),
                        output.lfsr = FALSE)

    ##### mapping the phi and sigma (rho) to main.score
    phi <- phi %>%
      dplyr::mutate(index = 1:nrow(phi)) %>%
      dplyr::select(.data$index, phi_mean = .data$mean, phi_sd = .data$sd)
    sigma2 <- sigma2 %>% dplyr::select(sigma2_mean = .data$mean, sigma2_sd = .data$sd)
    df_pos_index <- cbind(phi, sigma2)
    if (pos.act) {
      rho <- rho %>% dplyr::select(rho_mean = .data$mean, rho_sd = .data$sd)
      df_pos_index <- cbind(df_pos_index, rho)
    }
    main.score <- main.score %>%
      dplyr::left_join(df_pos_index, by = c("index" = "index"))

    if (!is.na(blosum.label[1])) {
      ##### mapping nu to main.score
      nu <- nu %>%
        dplyr::mutate(blosum = 1:nrow(nu)) %>%
        dplyr::select(.data$blosum, nu_mean = .data$mean, nu_sd = .data$sd)
      main.score <- main.score %>%
        dplyr::left_join(nu, by = c("blosum" = "blosum"))
    }
  }

  # Create Score Object
  score <- MCMCCreateScore(object = object, main.score = main.score, diags = diags)
  score@method <- method

  if (debug) {
    return(list(score = score, fit = fit))
  } else {
    return(score)
  }
}

#' Compute lfsr for variables from CmdStanMCMC object
#'
#' @param fit CmdStanMCMC object
#' @param param.key character string to match parameter names (functinal score)
#'
#' @return numeric vector of lfsr for each parameter
#'
MCMCLfsr <- function(fit, param.key = "beta") {

  param.name <- fit$metadata()$variables
  param.name <- param.name[stringr::str_detect(param.name,
                                               paste("^", param.key, sep = ""))]
  param.draw <- fit$draws(variables = param.name, format = "draws_matrix")
  lfsr <- apply(param.draw, MARGIN = 2,
                function(x) {min(mean(x > .Machine$double.eps),
                                 mean(x < .Machine$double.eps))})
  return(lfsr)
}

#' Map AA substitution to BLOSUM90 matrix
#'
#' Stop codon is *. If no AA substitution found code as -7.
#'
#' @param wt.vec vector of wildtype AA
#' @param mut.vec vector of mutation AA
#' @param aa.code "single" or "triple" amino acid coding
#'
#' @return numeric vector blosum90 score
#'
MapBlosumScore <- function(wt.vec, mut.vec, aa.code = "single") {

  # blosum90 is internal data in R folder
  blosum_map <- function(wt, mut) {
    if ((wt %in% rownames(blosum90)) & (mut %in% colnames(blosum90))) {
      # 3 is the max for human AA missense
      # human AA synonymous starts at 5
      return(min(blosum90[wt, mut], 5))
    } else if (mut == "DEL") {
      return(-7)
    } else if (startsWith(mut, "INS")) {
      return(-8)
    } else {
      warning(paste("No blosum mapping for ", wt, " and ", mut, sep = ""))
      return(-9)
    }
  }
  v_blosum_map <- Vectorize(blosum_map)

  # aa_table is internal data in R folder
  aa_table$triple <- toupper(aa_table$triple)
  if (aa.code == "triple") {
    wt.vec <- as.character(mapply(function(x) {aa_table$single[aa_table$triple == x]}, x = toupper(wt.vec)))
    mut.vec <- as.character(mapply(function(x) {aa_table$single[aa_table$triple == x]}, x = toupper(mut.vec)))
  } else if (aa.code == "single") {
    wt.vec <- toupper(wt.vec)
    mut.vec <- toupper(mut.vec)
  } else {
    stop("Invalid entry of AA code: single or triple")
  }

  return(as.numeric(v_blosum_map(wt.vec, mut.vec)))
}


#' Map each variant to its position label (index)
#'
#' Threshold is tricky to choose.
#'
#' @param var.names character vector of variant names
#' @param pos.label character or numeric vector of true positions
#' @param ctrl.label vector of whether variant is control
#' @param stop.label vector of whether variant is stop/nonsense
#' @param blosum.label vector of blosum score of variants
#' @param thred integer, threshold for number of variants per position label (index)
#'
#' @return data.frame with columns 'variant', 'pos', 'index'
#'
varPosIndexMap <- function(var.names, pos.label, ctrl.label, stop.label, blosum.label = NA, thred = 10) {

  df_map <- data.frame(variants = var.names, pos = pos.label)
  n_pos <- df_map %>%
    dplyr::group_by(.data$pos) %>%
    dplyr::summarise(n_pos =  n())
  n_pos$index <- 0

  # group position into index
  n_pos <- n_pos[order(n_pos$pos), ] # order by position
  curr_index <- 1
  counter <- 0
  for (i in 1:nrow(n_pos)) {
    n_pos$index[i] <- curr_index
    counter <- counter + n_pos$n_pos[i]
    if (counter >= thred) {
      # Reset counter and increment curr_index
      counter <- 0
      curr_index <- curr_index + 1
    }
  }
  if (counter < thred) {
    n_pos$index[n_pos$index == curr_index] <- curr_index - 1
    counter <- 0
  }
  df_map <- df_map %>% dplyr::left_join(n_pos)

  # map synonymous mutation index
  n_syn_group <- max(n_pos$n_pos) - 1
  if (!is.na(ctrl.label[1]) && sum(ctrl.label) > 0) {
    df_map$ctrl <- ctrl.label

    counter <- 0
    for (i in which(df_map$ctrl)[order(df_map$pos[df_map$ctrl])]) {
      df_map$index[i] <- curr_index
      counter <- counter + 1
      if (counter >= n_syn_group) {
        counter <- 0
        curr_index <- curr_index + 1
      }
    }
    if ((counter < thred) && (sum(df_map$ctrl) >= thred)) {
      df_map$index[df_map$index == curr_index] <- curr_index - 1
    } else if (sum(df_map$ctrl) != 0) {
      curr_index <- curr_index + 1
    }
    # df_map <- df_map %>% dplyr::mutate(index = ifelse(.data$ctrl, curr_idx, .data$index))
  } else {
    df_map$ctrl <- FALSE
  }

  # map stop/nonsense mutation index
  n_syn_group <- max(n_pos$n_pos) - 1
  if (!is.na(stop.label[1]) && sum(stop.label) > 0) {
    df_map$stop <- stop.label

    counter <- 0
    for (i in which(df_map$stop)[order(df_map$pos[df_map$stop])]) {
      df_map$index[i] <- curr_index
      counter <- counter + 1
      if (counter >= n_syn_group) {
        counter <- 0
        curr_index <- curr_index + 1
      }
    }
    if ((counter < thred) && (sum(df_map$stop) >= thred)) {
      df_map$index[df_map$index == curr_index] <- curr_index - 1
    } else if (sum(df_map$stop) != 0) {
      curr_index <- curr_index + 1
    }
  } else {
    df_map$stop <- FALSE
  }

  if (!is.na(blosum.label[1])) {
    df_map$blosum_score <- blosum.label
    df_map$blosum <- blosum.label

    # change stop codon blosum to the synonymous one (5) for inference
    # stop codon should not have position-dependent aa substitution score
    df_map$blosum[df_map$stop] <- 5

    # combine blosum group with less than 20% position representation
    cove <- 0.2
    compute_coverage <- function(blo, pos, label) {
      return(mean(table(blo, pos)[as.character(label), ] > 0))
    }
    label <- sort(unique(df_map$blosum))
    label_rm <- c()
    for (i in 1:(length(label) - 1)) {
      if (compute_coverage(blo = df_map$blosum, pos = df_map$index, label[i]) < cove) {
        if ((i + 1) != length(label)) {
          df_map$blosum[df_map$blosum == label[i]] <- label[i + 1]
          label_rm <- c(label_rm, i)
        } else {
          label_rm <- c(label_rm, i)
          label_new <- label[-label_rm]
          df_map$blosum[df_map$blosum == label[i]] <- label_new[length(label_new) - 1]
        }
      }
    }
    # print(rowMeans(table(df_map$blosum, df_map$index) > 0))
    df_map$blosum <- dplyr::dense_rank(df_map$blosum) # make the labels non negative
    # print(rowMeans(table(df_map$blosum, df_map$index) > 0))
  }

  return(df_map)
}

#' Impute Assays count
#'
#' If replicates are with different number of rounds,
#' truncate counts beyond the miminum number of rounds.
#' Otherwise, impute with mean counts of replicates.
#'
#' @param counts A matrix of counts
#' @param rounds A vector of number of rounds
#'
#' @return A list of imputed counts and rounds
#'
imputeAssaysCount <- function(counts, rounds) {

  if (sum(is.na(counts)) == 0) {
    print("No count imputation needed for combined.assays.")
    return(list(counts = counts,
                rounds = rounds))
  }

  # if replicates are with different number of rounds
  if (length(unique(rounds)) >= 2) {
    warnings("Impute assays counts with different number of rounds.
             Truncate counts beyond miminum number of rounds.
             Ex: [2, 3, 3] to [2, 2, 2].",
             immediate. = TRUE)

    end <- cumsum(rounds + 1)
    n_del <- rounds - min(rounds)
    idx_del <- c()
    for (i in 1:length(rounds)) {
      if (n_del[i] != 0) {
        idx_del <- c(idx_del, (end[i] - n_del[i] + 1):end[i])
      }
    }
    rm(end, n_del)

    counts <- counts[, -idx_del]
    rounds <- rep(min(rounds), length(rounds))

    if (ncol(counts) != (sum(rounds + 1))) {
      stop("imputeAssaysCount: error in count column removal.")
    }
  }

  # impute with new counts and rounds
  miss_idx <- which(rowSums(is.na(counts)) > 0)
  for (i in miss_idx) {
    row <- counts[i, ]
    # row to matrix
    data <- matrix(row, ncol = max(rounds + 1), byrow = TRUE)
    # column mean input
    for(j in 1:ncol(data)){
      data[is.na(data[, j]), j] <- mean(data[, j], na.rm = TRUE)
    }
    # matrix to row
    counts[i, ] <- c(t(data))
  }

  if (sum(is.na(counts))) {
    stop("imputeAssaysCount: error in imputation. NA still exists.")
  }

  return(list(counts = counts,
              rounds = rounds))
}

#' Impute Assays count (KNN Package)
#'
#' @param counts A matrix of counts
#' @param rounds A vector of number of rounds
#'
#' @return A list of imputed counts and rounds
#'
imputeAssaysCountKNN <- function(counts, rounds) {

  # TODO: change row max heristics 0.99 here
  mat_impute <- impute::impute.knn(
    counts,
    k = 10,
    rowmax = 0.99,
    maxp = nrow(counts),
    rng.seed = 362436069
  )
  mat_impute <- mat_impute$data

  return(list(counts = mat_impute,
              rounds = rounds))
}

