# wxWrapperVariationPerFeature
# This function is a wrapper of wxVariationFeature, which
# computes the amount of variation per point and per feature
# José Luis Alatorre Warren, University of Zurich, 2018

wxWrapperVariationPerFeature <- function(anatomicalFeatures,ranges) {
  
  # Create a list by copying anatomicalFeatures
  results <- anatomicalFeatures
  
  # Total number of features
  numberOfFeatures <- length(anatomicalFeatures)

  # Now that we have the same names, empty each field
  for (ii in 1:numberOfFeatures){
    results[[ii]] <- list()
  }
  
  # Create ranges to index members of each species (Hs and Pt)
  rangesOfSpecies    <- list()
  rangesOfSpecies$h  <- ranges$dense$rl$all$H
  rangesOfSpecies$p  <- ranges$dense$rl$all$P  

  for (ii in 1:numberOfFeatures){
    
    results[[ii]]$hp$dense    <- wxVariationPerFeature(anatomicalFeatures[[ii]]$dense)
    results[[ii]]$hp$sparse   <- wxVariationPerFeature(anatomicalFeatures[[ii]]$sparse)
    results[[ii]]$hp$extremes <- wxVariationPerFeature(anatomicalFeatures[[ii]]$extremes)
    results[[ii]]$h $dense    <- wxVariationPerFeature(anatomicalFeatures[[ii]]$dense   [,,rangesOfSpecies$h])
    results[[ii]]$h $sparse   <- wxVariationPerFeature(anatomicalFeatures[[ii]]$sparse  [,,rangesOfSpecies$h])
    results[[ii]]$h $extremes <- wxVariationPerFeature(anatomicalFeatures[[ii]]$extremes[,,rangesOfSpecies$h])
    results[[ii]]$p $dense    <- wxVariationPerFeature(anatomicalFeatures[[ii]]$dense   [,,rangesOfSpecies$p])
    results[[ii]]$p $sparse   <- wxVariationPerFeature(anatomicalFeatures[[ii]]$sparse  [,,rangesOfSpecies$p])
    results[[ii]]$p $extremes <- wxVariationPerFeature(anatomicalFeatures[[ii]]$extremes[,,rangesOfSpecies$p])
    
  }
  
  return(results)
    
}