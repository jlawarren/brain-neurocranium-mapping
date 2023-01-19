# wxVariationOfPairsOfFeatures
# This function is a wrapper of wxVariationFeature, which
# computes the amount of variation per point and per feature
# José Luis Alatorre Warren, University of Zurich, 2018

wxVariationOfPairsOfFeatures <- function(anatomicalFeatures,variationResults) {
 
  # Generate all combinations of n elements taken k at a time
  numberOfFeatures     <- length(anatomicalFeatures)
  matrixOfCombinations <- wxGenerateCombinations(nElements=numberOfFeatures,kGroupSize=2)  
  numberOfCombinations <- dim(matrixOfCombinations)[1]
  
  # Preallocate space in matrix
  matrixOfVariationOfPairsOfFeatures <- matrix(0,nrow=numberOfCombinations,ncol=14)
  
  for (ii in 1:numberOfCombinations){
    
    # Get current feature indices
    currentIndexFeature1 <- as.numeric(matrixOfCombinations[ii,1])
    currentIndexFeature2 <- as.numeric(matrixOfCombinations[ii,2])
    
    # Fill column with feature number codes
    matrixOfVariationOfPairsOfFeatures[ii,1] <- currentIndexFeature1
    matrixOfVariationOfPairsOfFeatures[ii,2] <- currentIndexFeature2
    
    # Fill variation (standard deviation) per feature pair
    matrixOfVariationOfPairsOfFeatures[ii, 3] <- variationResults[[currentIndexFeature1]]$hp$sparse  $standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii, 4] <- variationResults[[currentIndexFeature1]]$h $sparse  $standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii, 5] <- variationResults[[currentIndexFeature1]]$p $sparse  $standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii, 6] <- variationResults[[currentIndexFeature2]]$hp$sparse  $standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii, 7] <- variationResults[[currentIndexFeature2]]$h $sparse  $standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii, 8] <- variationResults[[currentIndexFeature2]]$p $sparse  $standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii, 9] <- variationResults[[currentIndexFeature1]]$hp$extremes$standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii,10] <- variationResults[[currentIndexFeature1]]$h $extremes$standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii,11] <- variationResults[[currentIndexFeature1]]$p $extremes$standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii,12] <- variationResults[[currentIndexFeature2]]$hp$extremes$standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii,13] <- variationResults[[currentIndexFeature2]]$h $extremes$standardDeviationFeature
    matrixOfVariationOfPairsOfFeatures[ii,14] <- variationResults[[currentIndexFeature2]]$p $extremes$standardDeviationFeature    
    
  }
  
  return(matrixOfVariationOfPairsOfFeatures)
   
}