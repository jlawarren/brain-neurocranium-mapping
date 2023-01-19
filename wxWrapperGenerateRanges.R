# wxWrapperGenerateRanges
# This functions is a wrapper of wxGenerateRanges, which creates
# index ranges to determine species (Homo sapiens or Pan troglodytes)
# and morphological module (brain or neurocranium)
# José Luis Alatorre Warren, University of Zurich, 2017

wxWrapperGenerateRanges <- function(procrustesSets) {
  
  # Create variable list
  ranges <- list()
  
  # Generate ranges from Procrustes sets
  ranges$dense$rl$all    <- wxGenerateRanges(procrustesSets$dense,   side='rl',submodule='all')
  ranges$dense$rl$par    <- wxGenerateRanges(procrustesSets$dense,   side='rl',submodule='par')
  ranges$dense$r$all     <- wxGenerateRanges(procrustesSets$dense,   side='r', submodule='all')
  ranges$dense$r$par     <- wxGenerateRanges(procrustesSets$dense,   side='r', submodule='par')
  ranges$dense$l$all     <- wxGenerateRanges(procrustesSets$dense,   side='l', submodule='all')
  ranges$dense$l$par     <- wxGenerateRanges(procrustesSets$dense,   side='l', submodule='par')
  ranges$sparse$rl$all   <- wxGenerateRanges(procrustesSets$sparse,  side='rl',submodule='all')
  ranges$sparse$rl$par   <- wxGenerateRanges(procrustesSets$sparse,  side='rl',submodule='par')
  ranges$sparse$r$all    <- wxGenerateRanges(procrustesSets$sparse,  side='r', submodule='all')
  ranges$sparse$r$par    <- wxGenerateRanges(procrustesSets$sparse,  side='r', submodule='par')
  ranges$sparse$l$all    <- wxGenerateRanges(procrustesSets$sparse,  side='l', submodule='all')
  ranges$sparse$l$par    <- wxGenerateRanges(procrustesSets$sparse,  side='l', submodule='par')
  ranges$extremes$rl$all <- wxGenerateRanges(procrustesSets$extremes,side='rl',submodule='all')
  ranges$extremes$rl$par <- wxGenerateRanges(procrustesSets$extremes,side='rl',submodule='par')
  ranges$extremes$r$all  <- wxGenerateRanges(procrustesSets$extremes,side='r', submodule='all')
  ranges$extremes$r$par  <- wxGenerateRanges(procrustesSets$extremes,side='r', submodule='par')
  ranges$extremes$l$all  <- wxGenerateRanges(procrustesSets$extremes,side='l', submodule='all')
  ranges$extremes$l$par  <- wxGenerateRanges(procrustesSets$extremes,side='l', submodule='par')
  
  return(ranges)
  
  }