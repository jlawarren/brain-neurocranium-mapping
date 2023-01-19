# wxVariationPerFeature
# This function computes the amount of variation per point and per feature
# José Luis Alatorre Warren, University of Zurich, 2018

wxVariationPerFeature <- function(currentFeature) {
  
  # Create a list
  results <- list()
  
  # Parameters of currentFeature
  featureLength       <- dim(currentFeature)[1]
  coordinatesPerPoint <- dim(currentFeature)[2]
  
  # Preallocate space for relevant variables
  results$standardDeviationCoordinates <- matrix(0,nrow=featureLength,ncol=coordinatesPerPoint)

  # Compute standard deviation of each coordinate
  for (jj in 1:featureLength){
    for (kk in 1:coordinatesPerPoint){
      results$standardDeviationCoordinates[jj,kk] <- sd(currentFeature[jj,kk,])
    }
  }
  
  # Compute standard deviation averages per point and per feature
  results$standardDeviationPoints  <- as.numeric(rowMeans(results$standardDeviationCoordinates))
  results$standardDeviationFeature <- as.numeric(mean(results$standardDeviationPoints))
  
  return(results)
  
}