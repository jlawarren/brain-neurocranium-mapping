# wxRearrengeSubmodulesTopology
# This function rearranges data such that the nodes "sparse" and "dense"
# are the first parents and not the last children
# José Luis Alatorre Warren, University of Zurich, 2017

wxRearrengeSubmodulesTopology <- function(listOfSubmoduleFeatures) {
  
  numberOfSubmoduleFeatures <- length(listOfSubmoduleFeatures)

  for (ii in 1:numberOfSubmoduleFeatures){

    currentFeatureDense    <- listOfSubmoduleFeatures[[ii]][["dense"]]
    currentFeatureSparse   <- listOfSubmoduleFeatures[[ii]][["sparse"]]
    currentFeatureExtremes <- listOfSubmoduleFeatures[[ii]][["extremes"]]
    
    if (ii == 1){
      dense    <- currentFeatureDense
      sparse   <- currentFeatureSparse
      extremes <- currentFeatureExtremes
    } else {
      dense    <- abind(dense,  currentFeatureDense,  along=1)
      sparse   <- abind(sparse, currentFeatureSparse, along=1)
      extremes <- abind(extremes, currentFeatureExtremes, along=1)
    }

  }
  
  procrustesSets          <- list()
  procrustesSets$dense    <- dense
  procrustesSets$sparse   <- sparse
  procrustesSets$extremes <- extremes
  return(procrustesSets)

}