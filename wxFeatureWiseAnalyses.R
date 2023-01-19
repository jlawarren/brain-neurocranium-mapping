# wxFeatureWiseAnalyses
# This function computes  
# Principal Component Analysis (PCA),
# Multivariate Multiple Regression (MMR),
# Two-Block Partial Least Squares (PLS) regression,
# and Adams's CR coefficient on a feature wise level
# José Luis Alatorre Warren, University of Zurich, 2017

wxFeatureWiseAnalyses <- function(anatomicalFeatures,rangeOfSpecies) {

  # Generate all combinations of n elements taken k at a time
  numberOfFeatures <- length(anatomicalFeatures)
  matrixOfCombinations <- wxGenerateCombinations(nElements=numberOfFeatures,kGroupSize=2)  
  numberOfCombinations <- dim(matrixOfCombinations)[1]
  
  # Create matrix with results
  matrixOfResults <- matrix(list(),numberOfCombinations,20)
  
  # Fill combinations in main table
  matrixOfResults[,1:2] <- matrixOfCombinations
  
  # Create intial variables
  currentIndices <- matrix(0,2,1)
  numberOfPcaComponentsUsed <- 6
  currentFeatures <- list()
  chainedFeatures <- list()
  partitionVector <- list()
  currentLengths  <- list()
  
  # Create a list containing all results
  treeOfResults          <- list()
  treeOfResults$sparse   <- list()
  treeOfResults$extremes <- list()  

  for (ii in 1:numberOfCombinations){
    
    # Print values to see progress of computation
    print(ii)
    print(length(rangeOfSpecies))

    # Expand list
    treeOfResults$sparse[[ii]]   <- list()
    treeOfResults$extremes[[ii]] <- list()    
    
    # Write names of current features
    treeOfResults$sparse  [[ii]]$featureNames$f1 <- names(anatomicalFeatures[as.numeric(matrixOfCombinations[ii,1])])
    treeOfResults$sparse  [[ii]]$featureNames$f2 <- names(anatomicalFeatures[as.numeric(matrixOfCombinations[ii,2])])
    treeOfResults$extremes[[ii]]$featureNames$f1 <- names(anatomicalFeatures[as.numeric(matrixOfCombinations[ii,1])])
    treeOfResults$extremes[[ii]]$featureNames$f2 <- names(anatomicalFeatures[as.numeric(matrixOfCombinations[ii,2])])
    
    # Update current indices
    currentIndices[1] <- as.numeric(matrixOfCombinations[ii,1])
    currentIndices[2] <- as.numeric(matrixOfCombinations[ii,2])
    
    # Update current features
    currentFeatures[[1]] <- anatomicalFeatures[[currentIndices[1]]]
    currentFeatures[[2]] <- anatomicalFeatures[[currentIndices[2]]]
    
    # Select the relevant specimens (both species, or only humans or only chimps), and discard the others
    currentFeatures[[1]]$dense    <- currentFeatures[[1]]$dense   [,,rangeOfSpecies]
    currentFeatures[[2]]$dense    <- currentFeatures[[2]]$dense   [,,rangeOfSpecies]
    currentFeatures[[1]]$sparse   <- currentFeatures[[1]]$sparse  [,,rangeOfSpecies]
    currentFeatures[[2]]$sparse   <- currentFeatures[[2]]$sparse  [,,rangeOfSpecies]
    currentFeatures[[1]]$extremes <- currentFeatures[[1]]$extremes[,,rangeOfSpecies]
    currentFeatures[[2]]$extremes <- currentFeatures[[2]]$extremes[,,rangeOfSpecies]

    # Glue features
    chainedFeatures$sparse   <- abind(currentFeatures[[1]]$sparse,  currentFeatures[[2]]$sparse,  along=1)
    chainedFeatures$extremes <- abind(currentFeatures[[1]]$extremes,currentFeatures[[2]]$extremes,along=1)
    
    # Compute PCA of each feature separately
    treeOfResults$sparse  [[ii]]$pca$f1 <- plotTangentSpace(currentFeatures[[1]]$sparse,  warpgrids=FALSE)
    treeOfResults$sparse  [[ii]]$pca$f2 <- plotTangentSpace(currentFeatures[[2]]$sparse,  warpgrids=FALSE)
    treeOfResults$extremes[[ii]]$pca$f1 <- plotTangentSpace(currentFeatures[[1]]$extremes,warpgrids=FALSE)
    treeOfResults$extremes[[ii]]$pca$f2 <- plotTangentSpace(currentFeatures[[2]]$extremes,warpgrids=FALSE)
    
    # Dimensionality reduction of PCA components
    treeOfResults$sparse  [[ii]]$pca$f1$pc.scores.reduced <- treeOfResults$sparse  [[ii]]$pca$f1$pc.scores[1:numberOfPcaComponentsUsed]
    treeOfResults$sparse  [[ii]]$pca$f2$pc.scores.reduced <- treeOfResults$sparse  [[ii]]$pca$f2$pc.scores[1:numberOfPcaComponentsUsed]
    treeOfResults$extremes[[ii]]$pca$f1$pc.scores.reduced <- treeOfResults$extremes[[ii]]$pca$f1$pc.scores[1:numberOfPcaComponentsUsed]
    treeOfResults$extremes[[ii]]$pca$f2$pc.scores.reduced <- treeOfResults$extremes[[ii]]$pca$f2$pc.scores[1:numberOfPcaComponentsUsed]
    
    # Multivariate Multiple Regression (MMR) from principal components computed from geomorph's procD.lm 
    # Regressing feature 1 (response) PC scores on feature 2 (predictor) PC scores
    feature1 <- treeOfResults$sparse  [[ii]]$pca$f1$pc.scores.reduced
    feature2 <- treeOfResults$sparse  [[ii]]$pca$f2$pc.scores.reduced
    gdf <- geomorph.data.frame(feature1=feature1,feature2=feature2)
    treeOfResults$sparse[[ii]]$mmr$responseF1 <- procD.lm(f1=feature1~feature2,RRPP=TRUE,effect.type=c("F","SS","cohen"),data=gdf,print.progress=TRUE)
    
    # Multivariate Multiple Regression (MMR) from principal components computed from geomorph's procD.lm 
    # Regressing feature 1 (response) PC scores on feature 2 (predictor) PC scores
    feature1 <- treeOfResults$sparse  [[ii]]$pca$f2$pc.scores.reduced
    feature2 <- treeOfResults$sparse  [[ii]]$pca$f1$pc.scores.reduced
    gdf <- geomorph.data.frame(feature1=feature1,feature2=feature2)
    treeOfResults$sparse[[ii]]$mmr$responseF2 <- procD.lm(f1=feature1~feature2,RRPP=TRUE,effect.type=c("F","SS","cohen"),data=gdf,print.progress=TRUE)
    
    # Multivariate Multiple Regression (MMR) from principal components computed from geomorph's procD.lm 
    # Regressing feature 1 (response) PC scores on feature 2 (predictor) PC scores
    feature1 <- treeOfResults$extremes[[ii]]$pca$f1$pc.scores.reduced
    feature2 <- treeOfResults$extremes[[ii]]$pca$f2$pc.scores.reduced
    gdf <- geomorph.data.frame(feature1=feature1,feature2=feature2)
    treeOfResults$extremes[[ii]]$mmr$responseF1 <- procD.lm(f1=feature1~feature2,RRPP=TRUE,effect.type=c("F","SS","cohen"),data=gdf,print.progress=TRUE)
    
    # Multivariate Multiple Regression (MMR) from principal components computed from geomorph's procD.lm 
    # Regressing feature 1 (response) PC scores on feature 2 (predictor) PC scores
    feature1 <- treeOfResults$extremes[[ii]]$pca$f2$pc.scores.reduced
    feature2 <- treeOfResults$extremes[[ii]]$pca$f1$pc.scores.reduced
    gdf <- geomorph.data.frame(feature1=feature1,feature2=feature2)
    treeOfResults$extremes[[ii]]$mmr$responseF2 <- procD.lm(f1=feature1~feature2,RRPP=TRUE,effect.type=c("F","SS","cohen"),data=gdf,print.progress=TRUE)

    # Compute lenghts of features
    currentLengths$f1$sparse   <- dim(currentFeatures[[1]]$sparse)[1]
    currentLengths$f2$sparse   <- dim(currentFeatures[[2]]$sparse)[1]
    currentLengths$f1$extremes <- dim(currentFeatures[[1]]$extremes)[1]
    currentLengths$f2$extremes <- dim(currentFeatures[[2]]$extremes)[1]
    
    # Preallocate partition vectors
    partitionVector$sparse   <- character(dim(chainedFeatures$sparse)[1])
    partitionVector$extremes <- character(dim(chainedFeatures$extremes)[1])
    
    # Fill partition vectors 
    partitionVector$sparse[1:currentLengths$f1$sparse] <- "A"
    partitionVector$sparse[(currentLengths$f1$sparse+1):dim(chainedFeatures$sparse)[1]] <- "B"
    partitionVector$extremes[1:currentLengths$f1$extremes] <- "A"
    partitionVector$extremes[(currentLengths$f1$extremes+1):dim(chainedFeatures$extremes)[1]] <- "B"

    # Compute morphological integration (RV) using a two-block Partial Least Squares (PLS) analysis
    treeOfResults$sparse  [[ii]]$integration <- integration.test(chainedFeatures$sparse,  partition.gp=partitionVector$sparse,  iter=999,seed=NULL,print.progress=TRUE)
    treeOfResults$extremes[[ii]]$integration <- integration.test(chainedFeatures$extremes,partition.gp=partitionVector$extremes,iter=999,seed=NULL,print.progress=TRUE)    
        
    # Compute morphological modularity test (CR)
    treeOfResults$sparse  [[ii]]$modularity <- modularity.test(chainedFeatures$sparse,  partition.gp=partitionVector$sparse,  iter=999,CI=FALSE,seed=NULL,print.progress=TRUE)
    treeOfResults$extremes[[ii]]$modularity <- modularity.test(chainedFeatures$extremes,partition.gp=partitionVector$extremes,iter=999,CI=FALSE,seed=NULL,print.progress=TRUE)    
    
  }
  
  return(treeOfResults)

}