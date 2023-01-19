# wxSimplify
# This function reduces the density of the anatomical features by a
# factor of 10, while preserving their two ends or extremes.
# It can also only grab the two ends of extremes of each feature.
# José Luis Alatorre Warren, University of Zurich, 2017

wxSimplify <- function(originalArray,simplificationType) {
  
  if (simplificationType=="sparse") {
    
    # The features will be simplified by a factor of 10
    reductionFactor <- 10
    
    # Main algorithm
    featureLength         <- dim(originalArray)[1]
    originalIndexVector   <- 1:featureLength
    divisionRemainder     <- featureLength%%reductionFactor
    startingPoint         <- 1 + divisionRemainder%/%2
    simplifiedIndexVector <- originalIndexVector[seq(startingPoint, featureLength, reductionFactor)]
    
    # Correction for remainders equal to 2 
    if(divisionRemainder==2){
      simplifiedIndexVector[1] <- 1
    }
    
    # Append first point of the feature
    if ((simplifiedIndexVector[1]!=1)){
      simplifiedIndexVector <- append(1,simplifiedIndexVector)
    }
    
    # Append last point of the feature
    if (tail(simplifiedIndexVector,n=1)!=featureLength){
      simplifiedIndexVector <- append(simplifiedIndexVector,featureLength)
    }
    
    # Create simplified array
    simplifiedArray <- originalArray[simplifiedIndexVector,,]
  
  } else if (simplificationType=="extremes") {
    
    # Just grab the first and the last points of each specimen
    lastIndex <- dim(originalArray)[1]
    simplifiedArray <- originalArray[c(1,lastIndex),,]
    
  }
  
  return(simplifiedArray)

}