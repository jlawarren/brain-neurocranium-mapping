# wxWrapperTreeOfAnalyses
# This functions is a wrapper of wxTreeOfAnalyses 
# The function wxTreeOfAnalyses computes Principal Component Analysis (PCA),
# Multivariate Multiple Regression (MMR), and Two-Block
# Partial Least Squares (PLS) regression on 2 populations
# with 2 morphological modules each.
# José Luis Alatorre Warren, University of Zurich, 2017

wxWrapperTreeOfAnalyses <- function(alignedData,ranges) {

  # Create list
  treeOfResults <- list()
  
  # Generate tree of results
  treeOfResults$dense$rl$all    <- wxTreeOfAnalyses(procrustesSets$dense$rl$bn$all,   ranges$dense$rl$all,   boolSimplified = FALSE)
  treeOfResults$dense$rl$par    <- wxTreeOfAnalyses(procrustesSets$dense$rl$bn$par,   ranges$dense$rl$par,   boolSimplified = FALSE) 
  treeOfResults$dense$r$all     <- wxTreeOfAnalyses(procrustesSets$dense$r$bn$all,    ranges$dense$r$all,    boolSimplified = FALSE) 
  treeOfResults$dense$r$par     <- wxTreeOfAnalyses(procrustesSets$dense$r$bn$par,    ranges$dense$r$par,    boolSimplified = FALSE) 
  treeOfResults$dense$l$all     <- wxTreeOfAnalyses(procrustesSets$dense$l$bn$all,    ranges$dense$l$all,    boolSimplified = FALSE) 
  treeOfResults$dense$l$par     <- wxTreeOfAnalyses(procrustesSets$dense$l$bn$par,    ranges$dense$l$par,    boolSimplified = FALSE)
  treeOfResults$sparse$rl$all   <- wxTreeOfAnalyses(procrustesSets$sparse$rl$bn$all,  ranges$sparse$rl$all,  boolSimplified = TRUE)
  treeOfResults$sparse$rl$par   <- wxTreeOfAnalyses(procrustesSets$sparse$rl$bn$par,  ranges$sparse$rl$par,  boolSimplified = TRUE) 
  treeOfResults$sparse$r$all    <- wxTreeOfAnalyses(procrustesSets$sparse$r$bn$all,   ranges$sparse$r$all,   boolSimplified = TRUE) 
  treeOfResults$sparse$r$par    <- wxTreeOfAnalyses(procrustesSets$sparse$r$bn$par,   ranges$sparse$r$par,   boolSimplified = TRUE) 
  treeOfResults$sparse$l$all    <- wxTreeOfAnalyses(procrustesSets$sparse$l$bn$all,   ranges$sparse$l$all,   boolSimplified = TRUE) 
  treeOfResults$sparse$l$par    <- wxTreeOfAnalyses(procrustesSets$sparse$l$bn$par,   ranges$sparse$l$par,   boolSimplified = TRUE) 
  treeOfResults$extremes$rl$all <- wxTreeOfAnalyses(procrustesSets$extremes$rl$bn$all,ranges$extremes$rl$all,boolSimplified = TRUE)
  treeOfResults$extremes$rl$par <- wxTreeOfAnalyses(procrustesSets$extremes$rl$bn$par,ranges$extremes$rl$par,boolSimplified = TRUE) 
  treeOfResults$extremes$r$all  <- wxTreeOfAnalyses(procrustesSets$extremes$r$bn$all, ranges$extremes$r$all, boolSimplified = TRUE) 
  treeOfResults$extremes$r$par  <- wxTreeOfAnalyses(procrustesSets$extremes$r$bn$par, ranges$extremes$r$par, boolSimplified = TRUE) 
  treeOfResults$extremes$l$all  <- wxTreeOfAnalyses(procrustesSets$extremes$l$bn$all, ranges$extremes$l$all, boolSimplified = TRUE) 
  treeOfResults$extremes$l$par  <- wxTreeOfAnalyses(procrustesSets$extremes$l$bn$par, ranges$extremes$l$par, boolSimplified = TRUE)
  
  return(treeOfResults)
  
}