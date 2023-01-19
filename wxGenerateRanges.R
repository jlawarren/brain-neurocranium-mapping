# wxGenerateRanges
# This function creates index ranges to determine species (Homo sapiens or Pan troglodytes)
# and morphological module (brain or neurocranium)
# José Luis Alatorre Warren, University of Zurich, 2017

wxGenerateRanges <- function(procrustesSets,side,submodule) {

  # Create ranges from the main dataset (hpbn)
  # H (Homo), P (Pan), B (Brain), N (Neurocranium)
  # Homo sapiens: specimens from 1 to 41
  # Pan troglodytes: specimens from 42 to 65
  # Brain landmarks (full dataset): from 1 to 1573
  # Neurocranium landmarks (full dataset): 1574 to 3258
  ranges <- list()
  
  # Species range
  ranges[["H"]] <-  1:41;
  ranges[["P"]] <- 42:65;
  
  # Module range
  
  if (side == 'rl' && submodule == 'all'){
    
    ranges[["B"]] <- 1:dim(procrustesSets$rl$b$all)[1]
    ranges[["N"]] <- (dim(procrustesSets$rl$b$all)[1]+1):(dim(procrustesSets$rl$bn$all)[1])

  } else if (side == 'rl' && submodule == 'par'){
  
    ranges[["B"]] <- 1:dim(procrustesSets$rl$b$par)[1]
    ranges[["N"]] <- (dim(procrustesSets$rl$b$par)[1]+1):(dim(procrustesSets$rl$bn$par)[1])
    
  } else if (side == 'r' && submodule == 'all'){

    ranges[["B"]] <- 1:dim(procrustesSets$r$b$all)[1]
    ranges[["N"]] <- (dim(procrustesSets$r$b$all)[1]+1):(dim(procrustesSets$r$bn$all)[1])
    
  } else if (side == 'r' && submodule == 'par'){
    
    ranges[["B"]] <- 1:dim(procrustesSets$r$b$par)[1]
    ranges[["N"]] <- (dim(procrustesSets$r$b$par)[1]+1):(dim(procrustesSets$r$bn$par)[1])

  } else if (side == 'l' && submodule == 'all'){
    
    ranges[["B"]] <- 1:dim(procrustesSets$l$b$all)[1]
    ranges[["N"]] <- (dim(procrustesSets$l$b$all)[1]+1):(dim(procrustesSets$l$bn$all)[1])
    
  } else if (side == 'l' && submodule == 'par'){
    
    ranges[["B"]] <- 1:dim(procrustesSets$l$b$par)[1]
    ranges[["N"]] <- (dim(procrustesSets$l$b$par)[1]+1):(dim(procrustesSets$l$bn$par)[1])
    
  } else {
    stop("Error in variable 'side'")
  }

  return(ranges)

}