# wxGenerateCombinations
# Generate all combinations of n elements taken k at a time
# José Luis Alatorre Warren, University of Zurich, 2018

wxGenerateCombinations <- function(nElements,kGroupSize) {

# Generate all combinations of n elements taken k at a time
arrayOfNumbers <- as.character(1:nElements)
matrixOfCombinations <- t(combn(arrayOfNumbers,kGroupSize))

return(matrixOfCombinations)

}