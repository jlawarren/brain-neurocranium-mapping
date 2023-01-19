# wxGenerateProcrustesSet
# This function creates multiple morphometric sets from 2 species
# (Homo sapiens and Pan troglodytes) with 2 morphometric modules
# each (brain and neurocranium).
# José Luis Alatorre Warren, University of Zurich, 2017

wxGenerateProcrustesSets <- function(anatomicalFeatures,tableOfSubmodules) {

  # Create a list of brain and neurocranium submodules
  listOfSubmodules <- list()
  
  # Create boolean vectors to generate R and L hemispheres
  boolVectorLeft           <- ((tableOfSubmodules[,"Side"]=="L")|(tableOfSubmodules[,"Side"]=="M"))
  boolVectorRight          <- ((tableOfSubmodules[,"Side"]=="R")|(tableOfSubmodules[,"Side"]=="M"))
  
  # Create basic RL modules
  listOfSubmodules$rl$b$all <- anatomicalFeatures[which(tableOfSubmodules[,"b_all"]==1)]
  listOfSubmodules$rl$b$par <- anatomicalFeatures[which(tableOfSubmodules[,"b_par"]==1)]
  listOfSubmodules$rl$n$int <- anatomicalFeatures[which(tableOfSubmodules[,"n_int"]==1)]  
  listOfSubmodules$rl$n$ext <- anatomicalFeatures[which(tableOfSubmodules[,"n_ext"]==1)]  
  listOfSubmodules$rl$n$par <- anatomicalFeatures[which(tableOfSubmodules[,"n_par"]==1)]
  
  # Create basic L modules
  listOfSubmodules$l$b$all  <- anatomicalFeatures[which((tableOfSubmodules[,"b_all"]==1)&boolVectorLeft)]
  listOfSubmodules$l$b$par  <- anatomicalFeatures[which((tableOfSubmodules[,"b_par"]==1)&boolVectorLeft)]
  listOfSubmodules$l$n$int  <- anatomicalFeatures[which((tableOfSubmodules[,"n_int"]==1)&boolVectorLeft)]  
  listOfSubmodules$l$n$ext  <- anatomicalFeatures[which((tableOfSubmodules[,"n_ext"]==1)&boolVectorLeft)]  
  listOfSubmodules$l$n$par  <- anatomicalFeatures[which((tableOfSubmodules[,"n_par"]==1)&boolVectorLeft)]  
  
  # Create basic R modules
  listOfSubmodules$r$b$all  <- anatomicalFeatures[which((tableOfSubmodules[,"b_all"]==1)&boolVectorRight)]
  listOfSubmodules$r$b$par  <- anatomicalFeatures[which((tableOfSubmodules[,"b_par"]==1)&boolVectorRight)]
  listOfSubmodules$r$n$int  <- anatomicalFeatures[which((tableOfSubmodules[,"n_int"]==1)&boolVectorRight)]  
  listOfSubmodules$r$n$ext  <- anatomicalFeatures[which((tableOfSubmodules[,"n_ext"]==1)&boolVectorRight)]  
  listOfSubmodules$r$n$par  <- anatomicalFeatures[which((tableOfSubmodules[,"n_par"]==1)&boolVectorRight)]  
  
  # Create additional RL modules
  listOfSubmodules$rl$n$all  <- c(listOfSubmodules$rl$n$int,listOfSubmodules$rl$n$ext)
  listOfSubmodules$rl$bn$all <- c(listOfSubmodules$rl$b$all,listOfSubmodules$rl$n$all)
  listOfSubmodules$rl$bn$par <- c(listOfSubmodules$rl$b$par,listOfSubmodules$rl$n$par)
  
  # Create additional L modules
  listOfSubmodules$l$n$all   <- c(listOfSubmodules$l$n$int,listOfSubmodules$l$n$ext)
  listOfSubmodules$l$bn$all  <- c(listOfSubmodules$l$b$all,listOfSubmodules$l$n$all)
  listOfSubmodules$l$bn$par  <- c(listOfSubmodules$l$b$par,listOfSubmodules$l$n$par)
  
  # Create additional R modules
  listOfSubmodules$r$n$all   <- c(listOfSubmodules$r$n$int,listOfSubmodules$r$n$ext)
  listOfSubmodules$r$bn$all  <- c(listOfSubmodules$r$b$all,listOfSubmodules$r$n$all)
  listOfSubmodules$r$bn$par  <- c(listOfSubmodules$r$b$par,listOfSubmodules$r$n$par)
  
  # Convert list of submodules as a node and print the result
  print(as.Node(listOfSubmodules),limit=10000)
  
  # Create list with 3D arrays
  procrustesSets <- list()

  for (ii in 1:3){
    
    if (ii == 1){
      densityLabel <- "dense"
    }
    if (ii == 2){
      densityLabel <- "sparse"
    }
    if (ii == 3){
      densityLabel <- "extremes"
    }
    
    procrustesSets[[densityLabel]]$rl$bn$all  <- wxRearrengeSubmodulesTopology(listOfSubmodules$rl$bn$all)[[densityLabel]]
    procrustesSets[[densityLabel]]$rl$bn$par  <- wxRearrengeSubmodulesTopology(listOfSubmodules$rl$bn$par)[[densityLabel]]
    procrustesSets[[densityLabel]]$rl$b$all   <- wxRearrengeSubmodulesTopology(listOfSubmodules$rl$b$all)[[densityLabel]]
    procrustesSets[[densityLabel]]$rl$b$par   <- wxRearrengeSubmodulesTopology(listOfSubmodules$rl$b$par)[[densityLabel]]
    procrustesSets[[densityLabel]]$rl$n$all   <- wxRearrengeSubmodulesTopology(listOfSubmodules$rl$n$all)[[densityLabel]]
    procrustesSets[[densityLabel]]$rl$n$int   <- wxRearrengeSubmodulesTopology(listOfSubmodules$rl$n$int)[[densityLabel]]
    procrustesSets[[densityLabel]]$rl$n$ext   <- wxRearrengeSubmodulesTopology(listOfSubmodules$rl$n$ext)[[densityLabel]]
    procrustesSets[[densityLabel]]$rl$n$par   <- wxRearrengeSubmodulesTopology(listOfSubmodules$rl$n$par)[[densityLabel]]
    
    procrustesSets[[densityLabel]]$r$bn$all   <- wxRearrengeSubmodulesTopology(listOfSubmodules$r$bn$all)[[densityLabel]]
    procrustesSets[[densityLabel]]$r$bn$par   <- wxRearrengeSubmodulesTopology(listOfSubmodules$r$bn$par)[[densityLabel]]
    procrustesSets[[densityLabel]]$r$b$all    <- wxRearrengeSubmodulesTopology(listOfSubmodules$r$b$all)[[densityLabel]]
    procrustesSets[[densityLabel]]$r$b$par    <- wxRearrengeSubmodulesTopology(listOfSubmodules$r$b$par)[[densityLabel]]
    procrustesSets[[densityLabel]]$r$n$all    <- wxRearrengeSubmodulesTopology(listOfSubmodules$r$n$all)[[densityLabel]]
    procrustesSets[[densityLabel]]$r$n$int    <- wxRearrengeSubmodulesTopology(listOfSubmodules$r$n$int)[[densityLabel]]
    procrustesSets[[densityLabel]]$r$n$ext    <- wxRearrengeSubmodulesTopology(listOfSubmodules$r$n$ext)[[densityLabel]]
    procrustesSets[[densityLabel]]$r$n$par    <- wxRearrengeSubmodulesTopology(listOfSubmodules$r$n$par)[[densityLabel]]
    
    procrustesSets[[densityLabel]]$l$bn$all   <- wxRearrengeSubmodulesTopology(listOfSubmodules$l$bn$all)[[densityLabel]]
    procrustesSets[[densityLabel]]$l$bn$par   <- wxRearrengeSubmodulesTopology(listOfSubmodules$l$bn$par)[[densityLabel]]
    procrustesSets[[densityLabel]]$l$b$all    <- wxRearrengeSubmodulesTopology(listOfSubmodules$l$b$all)[[densityLabel]]
    procrustesSets[[densityLabel]]$l$b$par    <- wxRearrengeSubmodulesTopology(listOfSubmodules$l$b$par)[[densityLabel]]
    procrustesSets[[densityLabel]]$l$n$all    <- wxRearrengeSubmodulesTopology(listOfSubmodules$l$n$all)[[densityLabel]]
    procrustesSets[[densityLabel]]$l$n$int    <- wxRearrengeSubmodulesTopology(listOfSubmodules$l$n$int)[[densityLabel]]
    procrustesSets[[densityLabel]]$l$n$ext    <- wxRearrengeSubmodulesTopology(listOfSubmodules$l$n$ext)[[densityLabel]]
    procrustesSets[[densityLabel]]$l$n$par    <- wxRearrengeSubmodulesTopology(listOfSubmodules$l$n$par)[[densityLabel]]

  }

  return(procrustesSets)
  
}