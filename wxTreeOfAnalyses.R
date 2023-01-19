# wxTreeOfAnalyses
# This function computes Principal Component Analysis (PCA),
# Multivariate Multiple Regression (MMR), Two-Block
# Partial Least Squares (PLS) regression, and Adams' CR
# coefficient on 2 populations with 2 morphological modules each.
# José Luis Alatorre Warren, University of Zurich, 2017

wxTreeOfAnalyses <- function(alignedData,ranges,boolSimplified) {
  
  # Create lists
  # hpbn: Homo, Pan, Brain, Neurocranium
  # hbn:  Homo,      Brain, Neurocranium
  # pbn:        Pan, Brain, Neurocranium
  # hpb:  Homo, Pan, Brain
  # hpn:  Homo, Pan,        Neurocranium
  # hb:   Homo,      Brain
  # hn:   Homo,             Neurocranium
  # pb:        Pan,  Brain
  # pn:        Pan,         Neurocranium
  hpbn <- list()
  hbn  <- list()
  pbn  <- list()
  hpb  <- list()
  hpn  <- list()
  hb   <- list()
  hn   <- list()
  pb   <- list()
  pn   <- list()
  
  # h: Homo sapiens
  # p: Pan troglodytes
  # b: Brain
  # n: Neurocranium
  # Extract the corresponding Procrustes alignments from run110n_morpho_equalized
  hpbn[["procrustes"]] <- alignedData
  hbn[["procrustes"]]  <- hpbn$procrustes[,,ranges$H]
  pbn[["procrustes"]]  <- hpbn$procrustes[,,ranges$P]
  hpb[["procrustes"]]  <- hpbn$procrustes[ranges$B,,]
  hpn[["procrustes"]]  <- hpbn$procrustes[ranges$N,,]
  hb[["procrustes"]]   <-  hbn$procrustes[ranges$B,,]
  hn[["procrustes"]]   <-  hbn$procrustes[ranges$N,,]
  pb[["procrustes"]]   <-  pbn$procrustes[ranges$B,,]
  pn[["procrustes"]]   <-  pbn$procrustes[ranges$N,,]
  
  # Save hpbn$procrustes in disk as CSV file (if it does not exist yet)
  if(!file.exists(names$csvOriginalFile)){
    write.table(hpbn$procrustes,
                file = names$csvOriginalFile,
                append = FALSE,
                quote = FALSE,
                sep = " ",
                eol = "\n",
                na = "NA",
                dec = ".",
                row.names = FALSE,
                col.names = FALSE,
                qmethod = c("escape", "double"),
                fileEncoding = "")
  }
  
  # Save hpbn in disk as .mat file (MATLAB file) (if it does not exist yet)
  if(!file.exists(names$matOriginalFile)){
    writeMat(con = names$matOriginalFile, A=hpbn$procrustes)
  }  
  
  # Perform PCA of shape variation of the aligned data
  # The function plotTangentSpace of the R package "geomorph" performs PCA from aligned data.
  # On the other hand, the R package "Morpho" performs (using the function procSym)
  # the PCA analysis together with Procrustes alignment. Unfortunately, this impedes extracting
  # PCA from subsets of aligned points.
  # Do this for hpbn, hbn, pbn, hpb, hpn, hb, hn, pb and pn
  hpbn[["pca"]] <- plotTangentSpace(hpbn$procrustes,
                                    axis1 = 1,
                                    axis2 = 2,
                                    warpgrids = FALSE,
                                    mesh = NULL,
                                    label = NULL,
                                    groups = NULL,
                                    legend = FALSE)
  hbn[["pca"]]  <- plotTangentSpace(hbn$procrustes,
                                    axis1 = 1,
                                    axis2 = 2,
                                    warpgrids = FALSE,
                                    mesh = NULL,
                                    label = NULL,
                                    groups = NULL,
                                    legend = FALSE)
  pbn[["pca"]]  <- plotTangentSpace(pbn$procrustes,
                                    axis1 = 1,
                                    axis2 = 2,
                                    warpgrids = FALSE,
                                    mesh = NULL,
                                    label = NULL,
                                    groups = NULL,
                                    legend = FALSE)
  hpb[["pca"]]  <- plotTangentSpace(hpb$procrustes,
                                    axis1 = 1,
                                    axis2 = 2,
                                    warpgrids = FALSE,
                                    mesh = NULL,
                                    label = NULL,
                                    groups = NULL,
                                    legend = FALSE)
  hpn[["pca"]]  <- plotTangentSpace(hpn$procrustes,
                                    axis1 = 1,
                                    axis2 = 2,
                                    warpgrids = FALSE,
                                    mesh = NULL,
                                    label = NULL,
                                    groups = NULL,
                                    legend = FALSE)
  hb[["pca"]]   <- plotTangentSpace(hb$procrustes,
                                    axis1 = 1,
                                    axis2 = 2,
                                    warpgrids = FALSE,
                                    mesh = NULL,
                                    label = NULL,
                                    groups = NULL,
                                    legend = FALSE)
  hn[["pca"]]   <- plotTangentSpace(hn$procrustes,
                                    axis1 = 1,
                                    axis2 = 2,
                                    warpgrids = FALSE,
                                    mesh = NULL,
                                    label = NULL,
                                    groups = NULL,
                                    legend = FALSE)
  pb[["pca"]]   <- plotTangentSpace(pb$procrustes,
                                    axis1 = 1,
                                    axis2 = 2,
                                    warpgrids = FALSE,
                                    mesh = NULL,
                                    label = NULL,
                                    groups = NULL,
                                    legend = FALSE)
  pn[["pca"]]   <- plotTangentSpace(pn$procrustes,
                                    axis1 = 1,
                                    axis2 = 2,
                                    warpgrids = FALSE,
                                    mesh = NULL,
                                    label = NULL,
                                    groups = NULL,
                                    legend = FALSE)
  
  # Dimensionaly reduction
  numberOfDimensions <- 10;
  hpbn$pca[["pc.scores.reduced"]] <- hpbn$pca$pc.scores[,1:numberOfDimensions]
  hbn$pca[["pc.scores.reduced"]]  <-  hbn$pca$pc.scores[,1:numberOfDimensions]
  pbn$pca[["pc.scores.reduced"]]  <-  pbn$pca$pc.scores[,1:numberOfDimensions]
  hpb$pca[["pc.scores.reduced"]]  <-  hpb$pca$pc.scores[,1:numberOfDimensions]
  hpn$pca[["pc.scores.reduced"]]  <-  hpn$pca$pc.scores[,1:numberOfDimensions]
  hb$pca[["pc.scores.reduced"]]   <-   hb$pca$pc.scores[,1:numberOfDimensions]
  hn$pca[["pc.scores.reduced"]]   <-   hn$pca$pc.scores[,1:numberOfDimensions]
  pb$pca[["pc.scores.reduced"]]   <-   pb$pca$pc.scores[,1:numberOfDimensions]
  pn$pca[["pc.scores.reduced"]]   <-   pn$pca$pc.scores[,1:numberOfDimensions]
  remove(numberOfDimensions)
  
  # Create empty lists to allocate full and reduced Multivariate Multiple Regressions (MMR)
  hpbn[["mmr"]] <- list()
  hbn[["mmr"]]  <- list()
  pbn[["mmr"]]  <- list()
  
  # Multivariate Multiple Regression (MMR) using all principal components using lm
  # Regressing brain (full) PC scores (response) on neurocranium (full) PC scores (predictor)
  hpbn$mmr$lm[["allpcs"]] <- lm(formula = hpb$pca$pc.scores ~ hpn$pca$pc.scores)
  hbn$mmr$lm[["allpcs"]]  <- lm(formula =  hb$pca$pc.scores ~  hn$pca$pc.scores)
  pbn$mmr$lm[["allpcs"]]  <- lm(formula =  pb$pca$pc.scores ~  pn$pca$pc.scores)
  
  # Multivariate Multiple Regression (MMR) using a reduced set of principal components using lm 
  # Regressing brain (reduced) PC scores (response) on neurocranium (reduced) PC scores (predictor)
  hpbn$mmr$lm[["fewpcs"]] <- lm(formula = hpb$pca$pc.scores.reduced ~ hpn$pca$pc.scores.reduced)
  hbn$mmr$lm[["fewpcs"]]  <- lm(formula =  hb$pca$pc.scores.reduced ~  hn$pca$pc.scores.reduced)
  pbn$mmr$lm[["fewpcs"]]  <- lm(formula =  pb$pca$pc.scores.reduced ~  pn$pca$pc.scores.reduced)
  
  # Multivariate Multiple Regression (MMR) using all principal components using geomorph's procD.lm 
  # Regressing brain (full) PC scores (response) on neurocranium (full) PC scores (predictor)
  brain                         <- hpb$pca$pc.scores
  neurocranium                  <- hpn$pca$pc.scores
  gdf                           <- geomorph.data.frame(brain = brain,neurocranium = neurocranium)
  hpbn$mmr$geomorph[["allpcs"]] <- procD.lm(f1 = brain ~ neurocranium, iter = 999, seed = NULL,
                                            RRPP = TRUE, effect.type = c("F","SS", "cohen"),
                                            int.first = FALSE, data = gdf, print.progress = TRUE)
  brain                         <- hb$pca$pc.scores
  neurocranium                  <- hn$pca$pc.scores
  gdf                           <- geomorph.data.frame(brain = brain,neurocranium = neurocranium)
  hbn$mmr$geomorph[["allpcs"]]  <- procD.lm(f1 = brain ~ neurocranium, iter = 999, seed = NULL,
                                            RRPP = TRUE, effect.type = c("F","SS", "cohen"),
                                            int.first = FALSE, data = gdf, print.progress = TRUE)
  brain                         <- pb$pca$pc.scores
  neurocranium                  <- pn$pca$pc.scores
  gdf                           <- geomorph.data.frame(brain = brain,neurocranium = neurocranium)
  pbn$mmr$geomorph[["allpcs"]]  <- procD.lm(f1 = brain ~ neurocranium, iter = 999, seed = NULL,
                                            RRPP = TRUE, effect.type = c("F","SS", "cohen"),
                                            int.first = FALSE, data = gdf, print.progress = TRUE)
  brain                         <- hpb$pca$pc.scores.reduced
  neurocranium                  <- hpn$pca$pc.scores.reduced
  gdf                           <- geomorph.data.frame(brain = brain,neurocranium = neurocranium)
  hpbn$mmr$geomorph[["fewpcs"]] <- procD.lm(f1 = brain ~ neurocranium, iter = 999, seed = NULL,
                                            RRPP = TRUE, effect.type = c("F","SS", "cohen"),
                                            int.first = FALSE, data = gdf, print.progress = TRUE)
  brain                         <- hb$pca$pc.scores.reduced
  neurocranium                  <- hn$pca$pc.scores.reduced
  gdf                           <- geomorph.data.frame(brain = brain,neurocranium = neurocranium)
  hbn$mmr$geomorph[["fewpcs"]]  <- procD.lm(f1 = brain ~ neurocranium, iter = 999, seed = NULL,
                                            RRPP = TRUE, effect.type = c("F","SS", "cohen"),
                                            int.first = FALSE, data = gdf, print.progress = TRUE)
  brain                         <- pb$pca$pc.scores.reduced
  neurocranium                  <- pn$pca$pc.scores.reduced
  gdf                           <- geomorph.data.frame(brain = brain,neurocranium = neurocranium)
  pbn$mmr$geomorph[["fewpcs"]]  <- procD.lm(f1 = brain ~ neurocranium, iter = 999, seed = NULL,
                                            RRPP = TRUE, effect.type = c("F","SS", "cohen"),
                                            int.first = FALSE, data = gdf, print.progress = TRUE)  
  
  # Create empty lists to allocate Two-Block Partial Least Squares (PLS) regression
  hpbn[["pls"]] <- list()
  hbn[["pls"]]  <- list()
  pbn[["pls"]]  <- list()
  
  # PLS regressions using brain as response and the neurocranium as predictor
  # This is done for 3 cases: HsPt, Hs, Pt
  # And each case is based on 2 options: based on covariation matrix (cov) and on correalation matrix (cor)
  hpbn$pls[["cov"]] <- pls2B(hpb$procrustes,
                             hpn$procrustes,
                             tol = 1e-12,
                             same.config = TRUE,
                             rounds = 2,
                             useCor = FALSE,
                             mc.cores = parallel::detectCores())
  hpbn$pls[["cor"]] <- pls2B(hpb$procrustes,
                             hpn$procrustes,
                             tol = 1e-12,
                             same.config = TRUE,
                             rounds = 2,
                             useCor = TRUE,
                             mc.cores = parallel::detectCores())
  hbn$pls[["cov"]] <- pls2B(hb$procrustes,
                            hn$procrustes,
                            tol = 1e-12,
                            same.config = TRUE,
                            rounds = 2,
                            useCor = FALSE,
                            mc.cores = parallel::detectCores())
  hbn$pls[["cor"]] <- pls2B(hb$procrustes,
                            hn$procrustes,
                            tol = 1e-12,
                            same.config = TRUE,
                            rounds = 2,
                            useCor = TRUE,
                            mc.cores = parallel::detectCores())
  pbn$pls[["cov"]] <- pls2B(pb$procrustes,
                            pn$procrustes,
                            tol = 1e-12,
                            same.config = TRUE,
                            rounds = 2,
                            useCor = FALSE,
                            mc.cores = parallel::detectCores())
  pbn$pls[["cor"]] <- pls2B(pb$procrustes,
                            pn$procrustes,
                            tol = 1e-12,
                            same.config = TRUE,
                            rounds = 2,
                            useCor = TRUE,
                            mc.cores = parallel::detectCores())
  
  # Prepare variables needed to generate brain shapes from the PLS regression
  # The variable numberOfPlsComponents will determine how many brain landmark configurations
  # will be created from the PLS regression.
  # A brain configuration will be generated per PLS component: the 1st brain configuration
  # will correspond to the 1st PLS, the 2nd brain configuration to the 2nd PLS component, and
  # the nth brain configuration to the nth PLS component.
  # The variable numberOfStandardDeviationsUsed is the number of standard deviations that will
  # be used to create the brain configurations (per PLS component). In other words, if
  # the variable numberOfStandardDeviationsUsed is set to 2, a configuration will be generated
  # using standard deviation of 1 and another one with standard deviation of 2
  numberOfBrainPoints <- length(ranges$B)
  numberOfSpatialDimensions <- 3
  numberOfPlsComponents <- 3
  numberOfStandardDeviationsUsed <- 2
  
  # Create empty lists to allocate brain configurations (responses) from the PLS analysis
  # using neurocranial configurations (predictors)
  hpbn$pls[["brainShapes"]][["cov"]] <- list()
  hpbn$pls[["brainShapes"]][["cor"]] <- list()
  hbn$pls[["brainShapes"]][["cov"]]  <- list()
  hbn$pls[["brainShapes"]][["cor"]]  <- list()
  pbn$pls[["brainShapes"]][["cov"]]  <- list()
  pbn$pls[["brainShapes"]][["cor"]]  <- list()
  
  # Generating brain configurations (responses) from the PLS analysis 
  # using  neurocranial configurations (predictors)
  for(ii in seq(1, numberOfPlsComponents, by=1)){
    
    currentNamePls <- paste("pls",toString(ii),sep = "")
    
    for(jj in seq(1, numberOfStandardDeviationsUsed, by=1)){
      
      currentNameSd  <- paste("sd", toString(jj),sep = "")
      
      hpbn$pls$brainShapes$cov[[currentNamePls]][[currentNameSd]] <- plsCoVar(hpbn$pls$cov,i=ii,sdx=jj,sdy=jj) 
      hpbn$pls$brainShapes$cor[[currentNamePls]][[currentNameSd]] <- plsCoVar(hpbn$pls$cor,i=ii,sdx=jj,sdy=jj) 
      hbn$pls$brainShapes$cov[[currentNamePls]][[currentNameSd]]  <- plsCoVar( hbn$pls$cov,i=ii,sdx=jj,sdy=jj) 
      hbn$pls$brainShapes$cor[[currentNamePls]][[currentNameSd]]  <- plsCoVar( hbn$pls$cor,i=ii,sdx=jj,sdy=jj) 
      pbn$pls$brainShapes$cov[[currentNamePls]][[currentNameSd]]  <- plsCoVar( pbn$pls$cov,i=ii,sdx=jj,sdy=jj) 
      pbn$pls$brainShapes$cor[[currentNamePls]][[currentNameSd]]  <- plsCoVar( pbn$pls$cor,i=ii,sdx=jj,sdy=jj) 
      
    }
  }
  remove(ii)
  remove(jj)
  remove(currentNameSd)
  remove(currentNamePls)    
  remove(numberOfBrainPoints)
  remove(numberOfPlsComponents)
  remove(numberOfSpatialDimensions)
  remove(numberOfStandardDeviationsUsed)
  
  # Compute analysis of morphological modularity (CR) if and only if the data has been simplified
  # In other words, only if the data have been transformed from having dense landmarks to having sparse landmarks
  if (boolSimplified == TRUE){

    # Prepare partition vector: label with brain (B) landamarks and with neurocranial (N) landmarks
    # The partition vector will be used to compute the morphological integration (RV) and modularity (CR) tests
    partitionVector <- character(length=dim(hpbn$procrustes)[1])
    partitionVector[ranges$B] <- "B"
    partitionVector[ranges$N] <- "N"

    # Compute morphological integration (RV) using a two-block Partial Least Squares (PLS) analysis
    hpbn$integration <- integration.test(A=hpbn$procrustes,partition.gp=partitionVector,iter=999,seed=NULL,print.progress=TRUE)
    hbn$integration  <- integration.test(A=hbn$procrustes, partition.gp=partitionVector,iter=999,seed=NULL,print.progress=TRUE)
    pbn$integration  <- integration.test(A=pbn$procrustes, partition.gp=partitionVector,iter=999,seed=NULL,print.progress=TRUE)

    # Compute morphological modularity (CR)
    hpbn$modularity <- modularity.test(A=hpbn$procrustes,partition.gp=partitionVector,iter=999,CI=FALSE,seed=NULL,print.progress=TRUE)
    hbn$modularity  <- modularity.test(A=hbn$procrustes, partition.gp=partitionVector,iter=999,CI=FALSE,seed=NULL,print.progress=TRUE)
    pbn$modularity  <- modularity.test(A=pbn$procrustes, partition.gp=partitionVector,iter=999,CI=FALSE,seed=NULL,print.progress=TRUE)

  }
  
  # Pretty-print the structure of the lists. It is similar to the R function
  # The function list.tree is similar to the R function str
  underscoreLength <- 90
  underscoreString <- paste(replicate(underscoreLength, "_"), collapse = "")
  print(underscoreString,quote=FALSE); list.tree(hpbn,numbers=FALSE,attr.print=FALSE,size=FALSE)
  print(underscoreString,quote=FALSE); list.tree( hbn,numbers=FALSE,attr.print=FALSE,size=FALSE)
  print(underscoreString,quote=FALSE); list.tree( pbn,numbers=FALSE,attr.print=FALSE,size=FALSE)
  print(underscoreString,quote=FALSE); list.tree( hpb,numbers=FALSE,attr.print=FALSE,size=FALSE)
  print(underscoreString,quote=FALSE); list.tree( hpn,numbers=FALSE,attr.print=FALSE,size=FALSE)
  print(underscoreString,quote=FALSE); list.tree(  hb,numbers=FALSE,attr.print=FALSE,size=FALSE)
  print(underscoreString,quote=FALSE); list.tree(  hn,numbers=FALSE,attr.print=FALSE,size=FALSE)
  print(underscoreString,quote=FALSE); list.tree(  pb,numbers=FALSE,attr.print=FALSE,size=FALSE)
  print(underscoreString,quote=FALSE); list.tree(  pn,numbers=FALSE,attr.print=FALSE,size=FALSE)
  remove(underscoreLength)
  remove(underscoreString)
  
  # Create list with results and return its value
  results <- list()
  results$hpbn <- hpbn
  results$hbn <- hbn
  results$pbn <- pbn
  results$hpb <- hpb
  results$hpn <- hpn
  results$hb <- hb
  results$hn <- hn
  results$pb <- pb
  results$pn <- pn
  return(results)

}