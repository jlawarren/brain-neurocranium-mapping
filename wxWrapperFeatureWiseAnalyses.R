# wxWrapperFeatureWiseAnalyses
# This function computes is a wrapper of wxFeatureWiseAnalyses,
# which computes Principal Component Analysis (PCA),
# Multivariate Multiple Regression (MMR),
# Two-Block Partial Least Squares (PLS) regression,
# and Adams's CR coefficient on a feature wise level
# José Luis Alatorre Warren, University of Zurich, 2017

wxWrapperFeatureWiseAnalyses <- function(anatomicalFeatures,ranges) {
  
  # Create list
  treeOfResults <- list()
  
  # Create ranges to index members of each species (Hs and Pt)
  rangesOfSpecies <- list()
  rangesOfSpecies$hp <- 1:(length(ranges$dense$rl$all$H)+length(ranges$dense$rl$all$P))
  rangesOfSpecies$h  <- ranges$dense$rl$all$H
  rangesOfSpecies$p  <- ranges$dense$rl$all$P
  
  # Create feature-wise analyses for both species together and separately
  treeOfResults$hp <- wxFeatureWiseAnalyses(anatomicalFeatures,rangesOfSpecies$hp)
  treeOfResults$h  <- wxFeatureWiseAnalyses(anatomicalFeatures,rangesOfSpecies$h)
  treeOfResults$p  <- wxFeatureWiseAnalyses(anatomicalFeatures,rangesOfSpecies$p)
  
  return(treeOfResults)
  
}