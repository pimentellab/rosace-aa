#' @importFrom readr read_tsv
#' @importFrom jsonlite fromJSON
#' @importFrom rjson toJSON
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Write the config list in json format
#'
#' @param cfg A list of config
#' @param file File name
#' @param savedir Directory to save the file
#'
#' @return None
#'
#' @export
#'
WriteConfig <- function(cfg, file, savedir) {
  jsonData <- rjson::toJSON(cfg, indent = 2)
  write(jsonData, file = paste0(file.path(savedir, file),".json"))
}

#' Read a full config file
#' (with score.df, pos.df, var.dist as data frame)
#'
#' @param file Config file path
#'
#' @return A list of config
#'
#' @rdname ReadConfig
#' @export
#'
ReadFullConfig <- function(file) {
  cfg <- jsonlite::fromJSON(file)
  return(cfg)
}

#' Create the config list for simulation from a Rosette object
#'
#' @param object A Rosette object
#' @param n.sim Number of simulations
#' @param save.sim Directory to save the simulation data
#' @param type.sim Type of simulation: growth or binding
#' @param mode.sim Simulation model code (interger 1 to 5)
#' @param n.rep Number of replicates
#' @param wt.effect Wild-type effect (binding) or doubling rate (growth)
#' @param n.round Number of rounds
#' @param mode.rep Replicate mode: "bio" (default) or "tech"
#' @param seq.shrink Shrinkage factor for sequencing dispersion (times seq.disp)
#' @param seq.depth Sequencing depth, by default 200
#' @param pop.size cell population size per variant
#' @param lib.shrink Shrinkage factor for library dispersion
#' @param var.shrink Shrinkage factor for variance (var.dist)
#'
#' @return A list of config
#'
#' @export
#'
# (1): input
# (2): optional. missing okay. can be inferred from rosette.
# (3): optional. have default value.
# (4): no input. directy from rosette.
CreateConfig <-
  function(object,
           n.sim, save.sim, type.sim, mode.sim, # sim (1)
           n.rep, # exp (1)
           wt.effect, # effect (1)
           n.round, # exp: optional from rosette (2)
           mode.rep = "bio", # exp (3)
           seq.shrink = 2, seq.depth = 200,  # seq (3)
           pop.size = 100, lib.shrink = 2, # pop (3)
           var.shrink = 0.8 # effect (3)
           # seq.disp # seq: from rosette (4)
           # lib.disp # pop: from rosette (4)
           # score.df, pos.df, var.dist, rounds # effect: from rosette (4)
           ) {

    # sim config
    sim <- listSimConfig(n.sim = n.sim, save.sim = save.sim, type.sim = type.sim, mode.sim = mode.sim) # sim (1)

    # exp config
    if (missing(n.round)) {
      n.round <- object@rounds
    }
    exp <- listExpConfig(n.rep = n.rep, # exp (1)
                         mode.rep = mode.rep, # exp (3)
                         n.round = n.round) # exp: optional from rosette (2)

    # pop config 
    lib.disp <- object@disp.start
    pop <- listPopConfig(pop.size = pop.size, # pop (3)
                         lib.disp = lib.disp, # pop (4)
                         lib.shrink = lib.shrink) # pop (3)

    # seq config
    if (length(object@disp) == 0) {
      stop("Object has empty slot 'disp'. Run 'AddDisp' first.")
    }
    seq.disp <- object@disp
    seq <- listSeqConfig(seq.disp = seq.disp,  # seq: from rosette (4)
                         seq.shrink = seq.shrink, seq.depth = seq.depth) # seq (3)

    # effect config
    effect <-
      listEffectConfig(wt.effect = wt.effect, # effect (1)
                       score.df = object@score.df, pos.df = object@pos.df,
                       var.dist = object@var.dist, rounds = object@rounds, # effect: from rosette (4)
                       var.shrink = var.shrink) # effect (3)

    cfg <- listConfig(list.sim = sim,
                      list.exp = exp,
                      list.seq = seq,
                      list.pop = pop,
                      list.effect = effect)

    return(cfg)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

listConfig <- function(list.sim, list.exp, list.seq,
                       list.pop, list.effect) {
  return(list(sim = list.sim,
              exp = list.exp,
              seq = list.seq,
              pop = list.pop,
              effect = list.effect))
}

listSimConfig <- function(n.sim, save.sim, type.sim, mode.sim) {

  return(list(n.sim = n.sim, # number of simulation
              save.sim = save.sim, # saving directory
              type.sim = type.sim, # type of experiment: growth  
              mode.sim = mode.sim)) # mode of experiment: numeric 1-5

}

listExpConfig <- function(n.rep,
                          mode.rep = "bio",
                          n.round = 3) {

  return(list(n.rep = n.rep, # number of replicates
              mode.rep = mode.rep, # replicate mode
              n.round = n.round)) # number of rounds
}

listSeqConfig <- function(seq.disp, seq.shrink = 1.5, seq.depth = 200) {

  if (seq.shrink < 1) {
    stop("Seq config: sequencing dispersion shrinkage has to be no less than 1.")
  }

  return(list(seq.disp = seq.disp, # sequencing dispersion ***import from rosette
              seq.shrink = seq.shrink, # shrink dispersion
              seq.depth = seq.depth)) # sequencing depth
}

listPopConfig <- function(pop.size, lib.disp, lib.shrink) {

  if (lib.shrink < 1) {
    stop("Pop config: library dispersion shrinkage has to be no less than 1.")
  }

  return(list(pop.size = pop.size,
              lib.disp = lib.disp,
              lib.shrink = lib.shrink)) # population size

}

listEffectConfig <- function(var.dist, score.df, pos.df, rounds,  
                             wt.effect, var.shrink = 0.8) {


  if (!'data.frame' %in% class(var.dist)) {
    var.dist <- readr::read_tsv(var.dist)
  }

  if (var.shrink > 1) {
    stop("Effect config: var.shrink has to be no more than 1.")
  } else if (var.shrink < 0.5) {
    warnings("Effect config: var.shink is lower than 0.5.")
  }

  return(list(rounds = rounds, # max rounds of exp ***import from rosette
              wt.effect = wt.effect, # wt effect
              var.shrink = var.shrink,
              var.dist = var.dist, # DATAFRAME ***import from rosette
              score.df = score.df, # DATAFRAME ***import from rosette
              pos.df = pos.df)) # DATAFRAME ***import from rosette
}
