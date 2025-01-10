#' @import dplyr
#' @importFrom tidyr separate
#' @importFrom readr write_tsv
#' @importFrom stats rmultinom rnorm rpois rbinom
#' @importFrom compositions rDirichlet.rcomp
#' @importFrom rjson toJSON
#' @importFrom utils write.table
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Simulate data from a config file
#'
#' @param config A list of config
#' @param save.tsv A boolean indicating whether to save the data in tsv format
#' @param save.rosace A boolean indicating whether to save the Rosace object
#' @param save.enrich2 A boolean indicating whether to save the data for enrich2
#'
#' @return None
#'
#' @export
#'
runRosette <- function(config,
                       save.tsv = TRUE,
                       save.rosace = TRUE,
                       save.enrich2 = FALSE){

  # create save directory
  save.dir <- config[['sim']][['save.sim']]
  save.dir <- file.path(save.dir,
                        paste(config[["sim"]][["type.sim"]],
                              "_rep",  config[["exp"]][["n.rep"]],
                              "_rd", config[["exp"]][["n.round"]],
                              "_", config[["sim"]][["mode.sim"]],
                              sep = ""))
  if (!dir.exists(save.dir)) {
    dir.create(save.dir, recursive = TRUE)
  }

  # start simulation
  sim_rounds <- config[["sim"]][["n.sim"]]
  for(i in 1:sim_rounds){
    set.seed(i)
    print(paste("Starting simulation round", i))

    score_df <- GenerateEffect(cfg = config) # expected_label, expected_effect, true_effect
    counts_list <- GenerateCount(cfg = config, effects = score_df$expected_effect,
                                 var_index = score_df$index) # counts, sequencing

    sub.save.dir <- file.path(save.dir, paste("sim", i, sep = ""))
    if (!dir.exists(sub.save.dir)) {
      dir.create(sub.save.dir, recursive = TRUE)
    }

    output_tsv(score_df = score_df, counts_list = counts_list, save.dir = sub.save.dir)
    output_rosace(score_df = score_df, counts_list = counts_list,
                  Nrep = config[["exp"]][["n.rep"]], Nround = config[["exp"]][["n.round"]],
                  save.dir = sub.save.dir)
    output_enrich2(counts_list = counts_list,
                   Nrep = config[["exp"]][["n.rep"]], Nround = config[["exp"]][["n.round"]],
                   mode.sim = config[['sim']][["mode.sim"]], save.dir = sub.save.dir)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GenerateEffect <- function(cfg) {

  mode_sim <- cfg[["sim"]][["mode.sim"]]
  type_exp <- cfg[["sim"]][["type.sim"]]
  Nround <- cfg[["exp"]][["n.round"]]
  Nrep <- cfg[["exp"]][["n.rep"]]

  # expected effect distribution
  expected_effect_list <-
    expected_effect_distribution_list(var_para = cfg[["effect"]][["var.dist"]],
                                      var_shrink_normal = cfg[["effect"]][["var.shrink"]],
                                      effect_scale = Nround/cfg[["effect"]][["rounds"]],
                                      null_group = "Neutral")

  # generate expected effect
  score_df <- generate_expected_effect(expected_effect_list = expected_effect_list,
                                       score_df = cfg[["effect"]][["score.df"]],
                                       pos_df = cfg[["effect"]][["pos.df"]],
                                       var_dist = cfg[["effect"]][["var.dist"]],
                                       mode_sim = mode_sim)

  # calculate true effect depending on experiments
  score_df$true_effect <-
    generate_true_effect(expected_effect = score_df$expected_effect,
                         type_exp = type_exp,
                         wt_effect = cfg[["effect"]][["wt.effect"]],
                         Nround = Nround)

  return(score_df)
}

GenerateCount <- function(cfg, effects, var_index) {

  Nround <- cfg[["exp"]][["n.round"]]
  Nrep <- cfg[["exp"]][["n.rep"]]
  Ncell_pop <- cfg[["pop"]][["pop.size"]] # per variant
  type_exp <- cfg[["sim"]][["type.sim"]]
  Nvar <- length(effects)
  Ncell_pop <- Ncell_pop * Nvar

  # starting count distribution
  dist_start <- starting_count_distribution(Ncell_pop = Ncell_pop, Nvar = Nvar,
                                            disp_start = cfg[["pop"]][["lib.disp"]] * cfg[["pop"]][["lib.shrink"]])

  # generate starting counts
  counts <- generate_starting_counts(dist_start = dist_start,
                                     Nround = Nround, Nrep = Nrep, Nvar = Nvar,
                                     replicate_mode = cfg[["exp"]][["mode.rep"]],
                                     var_index = var_index)

  # run experiment
  if(type_exp == 'binding') {
    # TODO: test binding screen
    stop("Function not tested yet.")
    # counts <- run_binding(counts = counts,
    #                       effects = effects,
    #                       Ncell_pop = Ncell_pop, Nround = Nround, Nrep = Nrep,
    #                       messages = TRUE)
  } else if (type_exp == 'growth') {
    counts <- run_growth(counts = counts,
                         effects = effects,
                         Ncell_pop = Ncell_pop, Nround = Nround, Nrep = Nrep,
                         wt = cfg[["effect"]][["wt.effect"]],
                         messages = TRUE)
  } else{
    stop(paste('Invalid assay mode', type_exp))
  }

  # generate sequencing count
  sequencing <- generate_sequencing_counts(counts = counts,
                                           depth = cfg[["seq"]][["seq.depth"]],
                                           dispersion = cfg[["seq"]][["seq.disp"]] * cfg[["seq"]][["seq.shrink"]],
                                           Nrep = Nrep,
                                           Nround = Nround)

  return(list(counts = counts,
              sequencing = sequencing))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Generate column name for counts/effects
# Example: rep1.c_0
rep_index <- function(Nround, Nrep){
  index <- paste0(rep(paste0('rep', seq(1:Nrep)), each = Nround + 1),
                  rep(paste0('.c_', seq(0, Nround)), Nrep))
  return(index)
}


# A multinomial distribution generator for starting cell population
starting_count_distribution <- function(Ncell_pop, Nvar, disp_start){
  gen <- function() {
    p <- as.numeric(compositions::rDirichlet.rcomp(n = 1, alpha = rep(1/Nvar, Nvar) * disp_start * Nvar))
    distln <- stats::rmultinom(n = 1, size = Ncell_pop, prob = p)
    return(distln)
  }

  return(gen)
}

# Initialize cell count data.frame and generate starting cell population
generate_starting_counts <- function(dist_start,
                                     Nround, Nrep, Nvar,
                                     replicate_mode,
                                     var_index){
  counts_df <- data.frame(matrix(0, nrow = Nvar, ncol = (Nround + 1) * Nrep))
  colnames(counts_df) <- rep_index(Nround, Nrep)
  row.names(counts_df) <- var_index

  if (replicate_mode == 'bio'){
    # independent draws for each replicate
    counts_df[, grep("c_0$", colnames(counts_df))] <-
      replicate(Nrep, as.numeric(dist_start()), simplify = FALSE)
  } else if(replicate_mode == 'tech'){
    # same draw for all replicates
    counts_df[, grep("c_0$", colnames(counts_df))] <- as.numeric(dist_start())
  } else{
    stop(paste("Invalid replicate mode", replicate_mode))
  }

  return(counts_df)
}


# A normal distribution generator for expected effect
expected_effect_normal <- function(loc, scale) {
  gen <- function(size) {
    distn <- stats::rnorm(size, mean = loc, sd = scale)
    return(distn)
  }
  return(gen)
}
# zero generator for expected effect
expected_effect_zero <- function() {
  gen <- function(size) {
    return(rep(0, size))
  }
  return(gen)
}

# A list of distributions for ctrl and each variant group
expected_effect_distribution_list <- function(var_para,
                                              var_shrink_normal,
                                              effect_scale,
                                              null_group){
  syn_mean <- var_para[var_para$label == 'ctrl',]$mean

  expected_effect <- list()

  expected_effect <- lapply(1:nrow(var_para), function(i){
    row <- var_para[i, ]
    if (row$label == null_group || row$label == "ctrl") {
      expected_effect_zero()
    } else {
      expected_effect_normal(loc = (row$mean - syn_mean) * effect_scale,
                             scale = var_shrink_normal * row$sd * effect_scale^2)
    }
  })

  names(expected_effect) <- var_para$label
  return(expected_effect)
}

# Generate expected effect based on mutant distribution
generate_expected_effect <- function(expected_effect_list, score_df, pos_df,
                                     var_dist, mode_sim) {

  score_df$expected_label <- ""
  score_df$expected_effect <- 0
  score_df$expected_label[score_df$ctrl] <- "ctrl"

  var_dist_sorted <- var_dist[-1, ] %>% dplyr::arrange(.data$label)
  var_para <- var_dist_sorted$label
  var_prop <- var_dist_sorted$count
  var_prop <- var_prop / sum(var_prop)

  point2label <- function(point, ctrl, var_para, var_prop) {
    # take advantage of the order after sorting
    # "Neg" "Neutral" "Pos"
    cutoff <- quantile(point[!ctrl], probs = cumsum(var_prop[-1]))
    if (length(cutoff) == 1) {
      label_ <- dplyr::case_when(
        point <= cutoff ~ var_para[1],
        .default = var_para[2]
      )
    } else if (length(cutoff) == 2) {
      label_ <- dplyr::case_when(
        point <= cutoff[1] ~ var_para[1],
        point > cutoff[2] ~ var_para[3],
        .default = var_para[2]
      )
    } else {
      stop("Invalid variant distribution label. Please provide 'Pos', 'Neutral' and 'Neg'.")
    }
    label_[ctrl] <- "ctrl"
    return(label_)
  }

  ### generate label ###
  if (length(var_para) == 1) {
    warning("All the variants have the same label.")
    score_df$expected_label[!score_df$ctrl] <- var_para
  } else if (mode_sim == 1){
    # random variant label
    score_df$expected_label[!score_df$ctrl] <-
      sample(x = var_para, size = sum(!score_df$ctrl), prob = var_prop, replace = TRUE)
  } else if (mode_sim == 2) {
    # random position label
    pos_df$pos_label <- sample(x = var_para, size = nrow(pos_df), prob = var_prop, replace = TRUE)
    score_df <- score_df %>%
      dplyr::left_join(pos_df %>% dplyr::select(.data$pos, .data$pos_label)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(expected_label = ifelse(.data$expected_label == "", .data$pos_label, .data$expected_label)) %>%
      dplyr::ungroup()
  } else if (mode_sim == 3) {
    # position + global blosum
    score_df <- score_df %>%
      dplyr::group_by(.data$blosum) %>%
      dplyr::mutate(point_blosum = mean(.data$score, na.rm = TRUE)) %>%
      dplyr::left_join(pos_df %>% dplyr::select(.data$pos, point_pos = .data$mean)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(point = .data$point_blosum + .data$point_pos,
                    expected_label = point2label(.data$point, .data$ctrl, var_para, var_prop))
  } else if (mode_sim == 4) {
    # position + position-activated blosum (act4, random)
    # Warnings: proportion of position activation is 50%.
    pos_df$act4 <- sample(x = c(TRUE, FALSE), size = nrow(pos_df), prob = c(0.5, 0.5), replace = TRUE)
    score_df <- score_df %>%
      dplyr::left_join(pos_df %>% dplyr::select(.data$pos, point_pos = .data$mean, .data$act4)) %>%
      dplyr::group_by(.data$blosum, .data$act4) %>%
      dplyr::mutate(point_blosum = mean(.data$score, na.rm = TRUE)) %>%
      dplyr::ungroup() %>% 
      dplyr::rowwise() %>%
      dplyr::mutate(point_blosum = ifelse(.data$act4, .data$point_blosum, 0)) %>%
      dplyr::ungroup() %>% 
      dplyr::mutate(point = .data$point_blosum + .data$point_pos,
                    expected_label = point2label(.data$point, .data$ctrl, var_para, var_prop))
  } else if (mode_sim == 5) {
    # position activation has structure based on phi
    cutoff <- quantile(pos_df$mean, probs = c(0.25, 0.5, 0.25))
    pos_df <- pos_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(act5 = ifelse(.data$mean > cutoff[2] || .data$mean < cutoff[1], FALSE, TRUE)) %>%
      dplyr::ungroup() 
    score_df <- score_df %>%
      dplyr::left_join(pos_df %>% dplyr::select(.data$pos, point_pos = .data$mean, .data$act5)) %>%
      dplyr::group_by(.data$blosum, .data$act5) %>%
      dplyr::mutate(point_blosum = mean(.data$score, na.rm = TRUE)) %>%
      dplyr::ungroup() %>% 
      dplyr::rowwise() %>%
      dplyr::mutate(point_blosum = ifelse(.data$act5, .data$point_blosum, 0)) %>%
      dplyr::ungroup() %>% 
      dplyr::mutate(point = .data$point_blosum + .data$point_pos,
                    expected_label = point2label(.data$point, .data$ctrl, var_para, var_prop))
  } else {
    stop("Invalid mode of simulation. Enter san integer from 1 to 5.")
  }

  ### from label to effect size ###
  score_df <- score_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(expected_effect = expected_effect_list[[.data$expected_label]](size = 1)) %>%
    dplyr::ungroup()

  return(score_df)
}

# Generate true effects (helper function)
# probability being selected/growth rate depending on simulation type
calc_true_effect_from_expected <- function(expected_effect, type_exp, wt, Nround){

  if (type_exp == 'binding'){
    stop("No binding experiment simulation implemented.")
    # true_effect <- exp(expected_effect/Nround + log(wt))
    # if(sum(true_effect > 1) > 0){
    #   print(paste(sum(true_effect > 1), "variants have binding effect more than 1. Shrink to 1."))
    #   true_effect[true_effect > 1] <- 1
    # }
  }
  else if (type_exp == 'growth'){
    if (wt < 0) {
      wt_effect <- -1
      true_effect <- wt_effect * (1 + (expected_effect/log(2)/Nround/wt))
      if(sum(true_effect > 0) > 0){
        warning(sum(true_effect > 0), " variants have growth effect larger than 0. ",
                "Maximum true effect is ", max(true_effect),
                call. = FALSE, immediate. = TRUE)
        # true_effect[true_effect > 0] <- 0
      }
    }
    else if (wt > 0) {
      wt_effect <- 1
      true_effect = wt_effect * (1 + (expected_effect/log(2)/Nround/wt))
      if(sum(true_effect < 0) > 0){
        warning(sum(true_effect < 0), "variants have growth effect smaller than 0.",
                "Minimum true effect is ", min(true_effect),
                call. = FALSE, immediate. = TRUE)
        # true_effect[true_effect < 0] <- 0
      }
    }
    else{
      stop(paste("Invalid wild type doubling rate", wt))
    }
  }
  else{
    stop(paste("Invalid assay mode ", type_exp))
  }

  return(true_effect)
}

# Genrerate true effects: dataframe
generate_true_effect <- function(expected_effect,
                                 type_exp,
                                 wt_effect,
                                 Nround){
  true_effect <-
    calc_true_effect_from_expected(expected_effect = expected_effect,
                                   type_exp = type_exp,
                                   wt = wt_effect,
                                   Nround = Nround)
  return(true_effect)
}

# Resample cells to initial cell population
resample_counts <- function(count, depth, dispersion=0){

  if (dispersion > 0){
    alpha <- nrow(count) * dispersion * count / sum(count)
    var_zero <- subset(alpha, alpha[, 1] == 0)
    var_nonzero <- subset(alpha, alpha[, 1] != 0)
    p <- c(compositions::rDirichlet.rcomp(n = 1, alpha = var_nonzero[,1]))
    p <- data.frame(p, row.names = row.names(var_nonzero))
    colnames(p) <- colnames(var_zero)
    p <- rbind(p, var_zero)
    p <- p[rownames(count), 1]
  } else if(dispersion == 0){
    p <- count / sum(count)
    p <- p[, 1]
  } else{
    stop(paste("Invalid dispersion value", dispersion))
  }

  sample_counts <- stats::rmultinom(n = 1, size = length(p) * as.integer(depth), p)
  return(unname(sample_counts))
}

# # TODO: NOT TESTED!!!
# # Run binding screen
# run_binding <- function(counts, effects, Ncell_pop, Nround, Nrep, messages = TRUE){
#   for(rep in 1:Nrep){
#     for(round in 1:Nround){
#       for(row in 1:nrow(counts)){
#         Nprev <- counts[row, paste0('rep', rep, '.c_', round - 1)]
#         p <- effects[row, paste0('rep', rep)]
#         counts[row, paste0('rep', rep, '.c_', round)] <- stats::rbinom(1, Nprev, p)
#       }

#       # grow cells to original population
#       counts[paste0('rep', rep, '.c_', round)] <-
#         resample_counts(counts[paste0('rep', rep, '.c_', round)], depth = Ncell_pop / nrow(counts))

#       if(messages){
#         print(paste('Finished rep', rep, "round", round))
#       }
#     }
#   }

#   return(counts)
# }

# Run growth screen
run_growth <- function(counts, effects, Ncell_pop, Nround, Nrep, wt, messages = TRUE){

  wt_effect <- sign(wt)
  if (wt == 0) {
    stop(paste("Invalid wild type doubling rate", wt))
  }
  t <- wt * log(2) / wt_effect

  for(rep in 1:Nrep){
    for(round in 1:Nround){
      for(row in 1:nrow(counts)){
        c <- counts[row, paste0('rep', rep, '.c_', round - 1)]
        if(c > 0){
          counts[row, paste0('rep', rep, '.c_', round)] <-
            stats::rpois(1, c * exp(t * effects[row]))
          # if(wt > 0){
          #   counts[row, paste0('rep',rep,'.c_',round)] <-
          #     np.random.negative_binomial(c, exp(-1 * t * effects[row])) + c
          # }
          # if(wt < 0){
          #   counts[row, paste0('rep',rep,'.c_',round)] <-
          #     rpois(1, c * exp(t * effects[row]))
          # }
        } else{
          counts[row, paste0('rep', rep, '.c_', round)] <- 0
        }
      }

      # if exceeding initial cell population
      if(sum(counts[paste0('rep', rep, '.c_', round)]) > Ncell_pop){
        counts[paste0('rep', rep, '.c_', round)] <-
          as.numeric(resample_counts(counts[paste0('rep', rep, '.c_', round)],
                          depth = Ncell_pop / nrow(counts)))
      } else if(sum(counts[paste0('rep', rep, '.c_', round)] == 0) > 0.8 * nrow(counts)){
        stop(paste('Too many rounds of cell growth. 80% variants die at round', round))
      }

      if(messages){
        print(paste('Finished rep', rep, "round", round))
      }
    }
  }

  return(counts)
}

# Generate sequencing count with dispersion
generate_sequencing_counts <- function(counts, depth, dispersion, Nrep, Nround){

  seq_counts <- data.frame(matrix(0, nrow = nrow(counts), ncol = ncol(counts)))
  colnames(seq_counts) <- colnames(counts)
  row.names(seq_counts) <- row.names(counts)

  for(rep in 1:Nrep){
    for(round in 0:Nround){
      seq_counts[paste0('rep',rep,'.c_',round)] =
        as.numeric(resample_counts(counts[paste0('rep', rep, '.c_', round)], depth, dispersion))
    }
  }

  return(seq_counts)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal: output
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_tsv <- function(score_df, counts_list, save.dir) {

  sub.save.dir <- file.path(save.dir, "tsv")
  if (!dir.exists(sub.save.dir)) {
    dir.create(sub.save.dir, recursive = TRUE)
  }

  # output effects
  readr::write_tsv(score_df, file = file.path(sub.save.dir, "effects.tsv"))

  # output counts
  counts <- counts_list$counts
  sequencing <- counts_list$sequencing
  readr::write_tsv(counts, file = file.path(sub.save.dir, "cell_counts.tsv"))
  readr::write_tsv(sequencing, file = file.path(sub.save.dir, "sequencing_counts.tsv"))

}

output_rosace <- function(score_df, counts_list, Nrep, Nround, save.dir) {

  # saving directory
  sub.save.dir <- file.path(save.dir, "rosace")
  if (!dir.exists(sub.save.dir)) {
    dir.create(sub.save.dir, recursive = TRUE)
  }

  # create rosace object
  for (i in 1:Nrep) {
    idx_start <- (Nround + 1) * i - Nround
    idx_end <- (Nround + 1) * i
    assay <- CreateAssayObject(counts = as.matrix(counts_list$sequencing[idx_start:idx_end]),
                               var.names = rownames(counts_list$sequencing),
                               key = "simulation", rep = i, type = "growth")
    if (i == 1) {
      rosace <- CreateRosaceObject(object = assay)
    } else {
      rosace <- AddAssayData(object = rosace, assay = assay)
    }
  }

  # # add ground truth data
  # effects <- data.frame(variants = rownames(effects_list$effects),
  #                       expected_effect = effects_list$expected_effect,
  #                       expected_label = effects_list$expected_label)
  # rosace@misc <- list(effects = effects)

  # process var.data
  rosace@var.data <- score_df

  # save rosace object
  saveRDS(rosace, file = file.path(sub.save.dir, "rosace.rds"))

}

output_enrich2 <- function(counts_list, Nrep, Nround, mode.sim, save.dir) {

  # saving directory
  sub.save.dir <- file.path(save.dir, "enrich2")
  if (!dir.exists(sub.save.dir)) {
    dir.create(sub.save.dir, recursive = TRUE)
  }
  data.dir <- file.path(sub.save.dir, "data")
  if (!dir.exists(data.dir)) {
    dir.create(data.dir, recursive = TRUE)
  }
  result.dir <- file.path(sub.save.dir, "results")
  if (!dir.exists(result.dir)) {
    dir.create(result.dir, recursive = TRUE)
  }

  # add wild-type info
  ctrl_label <- endsWith(rownames(counts_list$sequencing), 'ctrl')

  # output counts
  for (i in 1:Nrep) {
    for (j in 0:Nround) {
      count <- counts_list$sequencing[(i - 1) * (Nround + 1) + j + 1]
      colnames(count) <- "count"
      count[nrow(count) + 1, 1] <- sum(count[[1]][ctrl_label])
      rownames(count)[nrow(count)] <- "_wt"
      utils::write.table(count, quote = FALSE, sep="\t",
                         file = file.path(data.dir, paste("count_rep", i, "_c", j, ".tsv", sep = "")))
    }
  }

  # output json file
  experiment <- list('name' = 'simulation', 'output directory' = result.dir, 'conditions' = list())
  condition <- list('name' = mode.sim, selections = list())
  for (i in 1:Nrep) {
    selection <- list('name' = paste("rep", i, sep = ""), 'libraries' = list())
    for (j in 0:Nround) {
      seqlib = list('counts file' = file.path(data.dir, paste("count_rep", i, "_c", j, ".tsv", sep = "")),
                    'identifiers' = list(),
                    'name' = paste("rep", i, "_c", j, sep = ""),
                    'report filtered read' = FALSE,
                    'timepoint' = j)
      selection$libraries[[j + 1]] <- seqlib
    }
    condition$selections[[i]] <- selection
  }
  experiment$conditions[[1]] <- condition

  # save json file
  jsonData <- rjson::toJSON(experiment, indent = 2)
  write(jsonData, file = file.path(sub.save.dir, "config.json"))
}


