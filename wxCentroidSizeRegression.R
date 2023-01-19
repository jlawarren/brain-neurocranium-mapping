# wxCentroidSizeRegression
# This functions computes a regression of brain and neurocranial
# on size (log centroid size)
# José Luis Alatorre Warren, University of Zurich, 2017

wxCentroidSizeRegression <- function(rawData,ranges) {
  
  # Create list with results
  results <- list()
  
  # Create ranges for BN configuration
  ranges$sparse$rl$all$BN   <- c(ranges$sparse  $rl$all$B,ranges$sparse  $rl$all$N)
  ranges$extremes$rl$all$BN <- c(ranges$extremes$rl$all$B,ranges$extremes$rl$all$N)
  
  # _______________________________________________________________________________________________
  # Create subsets (species and modules) with raw data
  
  # Use dense semilandmarks
  results$dense$hpbn$rawData    <- rawData
  results$dense$hpb$rawData     <- rawData[ranges$dense$rl$all$B,,]
  results$dense$hpn$rawData     <- rawData[ranges$dense$rl$all$N,,]
  results$dense$hbn$rawData     <- rawData[                     ,,ranges$dense$rl$all$H]
  results$dense$pbn$rawData     <- rawData[                     ,,ranges$dense$rl$all$P]
  results$dense$hb$rawData      <- rawData[ranges$dense$rl$all$B,,ranges$dense$rl$all$H]
  results$dense$hn$rawData      <- rawData[ranges$dense$rl$all$N,,ranges$dense$rl$all$H]
  results$dense$pb$rawData      <- rawData[ranges$dense$rl$all$B,,ranges$dense$rl$all$P]
  results$dense$pn$rawData      <- rawData[ranges$dense$rl$all$N,,ranges$dense$rl$all$P]

  # Use  semilandmarks
  results$sparse$hpbn$rawData   <- rawData[ranges$sparse$rl$all$BN,,]
  results$sparse$hpb$rawData    <- rawData[ranges$sparse$rl$all$B ,,]
  results$sparse$hpn$rawData    <- rawData[ranges$sparse$rl$all$N ,,]
  results$sparse$hbn$rawData    <- rawData[ranges$sparse$rl$all$BN,,ranges$sparse$rl$all$H]
  results$sparse$pbn$rawData    <- rawData[ranges$sparse$rl$all$BN,,ranges$sparse$rl$all$P]
  results$sparse$hb$rawData     <- rawData[ranges$sparse$rl$all$B ,,ranges$sparse$rl$all$H]
  results$sparse$hn$rawData     <- rawData[ranges$sparse$rl$all$N ,,ranges$sparse$rl$all$H]
  results$sparse$pb$rawData     <- rawData[ranges$sparse$rl$all$B ,,ranges$sparse$rl$all$P]
  results$sparse$pn$rawData     <- rawData[ranges$sparse$rl$all$N ,,ranges$sparse$rl$all$P]

  # Use extreme semilandmarks
  results$extremes$hpbn$rawData <- rawData[ranges$extremes$rl$all$BN,,]
  results$extremes$hpb$rawData  <- rawData[ranges$extremes$rl$all$B ,,]
  results$extremes$hpn$rawData  <- rawData[ranges$extremes$rl$all$N ,,]
  results$extremes$hbn$rawData  <- rawData[ranges$extremes$rl$all$BN,,ranges$extremes$rl$all$H]
  results$extremes$pbn$rawData  <- rawData[ranges$extremes$rl$all$BN,,ranges$extremes$rl$all$P]
  results$extremes$hb$rawData   <- rawData[ranges$extremes$rl$all$B ,,ranges$extremes$rl$all$H]
  results$extremes$hn$rawData   <- rawData[ranges$extremes$rl$all$N ,,ranges$extremes$rl$all$H]
  results$extremes$pb$rawData   <- rawData[ranges$extremes$rl$all$B ,,ranges$extremes$rl$all$P]
  results$extremes$pn$rawData   <- rawData[ranges$extremes$rl$all$N ,,ranges$extremes$rl$all$P]

  # _______________________________________________________________________________________________
  # Perform procrustes analysis as implemented in geormorph's gpagen

  # Use dense semilandmarks
  results$dense$hpbn$gpagen <- gpagen(results$dense$hpbn$rawData)
  results$dense$hpb$gpagen  <- gpagen(results$dense$hpb$rawData)
  results$dense$hpn$gpagen  <- gpagen(results$dense$hpn$rawData)
  results$dense$hbn$gpagen  <- gpagen(results$dense$hbn$rawData)
  results$dense$pbn$gpagen  <- gpagen(results$dense$pbn$rawData)
  results$dense$hb$gpagen   <- gpagen(results$dense$hb$rawData)
  results$dense$hn$gpagen   <- gpagen(results$dense$hn$rawData)
  results$dense$pb$gpagen   <- gpagen(results$dense$pb$rawData)
  results$dense$pn$gpagen   <- gpagen(results$dense$pn$rawData)
  
  # Use  semilandmarks
  results$sparse$hpbn$gpagen <- gpagen(results$sparse$hpbn$rawData)
  results$sparse$hpb$gpagen  <- gpagen(results$sparse$hpb$rawData)
  results$sparse$hpn$gpagen  <- gpagen(results$sparse$hpn$rawData)
  results$sparse$hbn$gpagen  <- gpagen(results$sparse$hbn$rawData)
  results$sparse$pbn$gpagen  <- gpagen(results$sparse$pbn$rawData)
  results$sparse$hb$gpagen   <- gpagen(results$sparse$hb$rawData)
  results$sparse$hn$gpagen   <- gpagen(results$sparse$hn$rawData)
  results$sparse$pb$gpagen   <- gpagen(results$sparse$pb$rawData)
  results$sparse$pn$gpagen   <- gpagen(results$sparse$pn$rawData)
  
  # Use extremes semilandmarks
  results$extremes$hpbn$gpagen <- gpagen(results$extremes$hpbn$rawData)
  results$extremes$hpb$gpagen  <- gpagen(results$extremes$hpb$rawData)
  results$extremes$hpn$gpagen  <- gpagen(results$extremes$hpn$rawData)
  results$extremes$hbn$gpagen  <- gpagen(results$extremes$hbn$rawData)
  results$extremes$pbn$gpagen  <- gpagen(results$extremes$pbn$rawData)
  results$extremes$hb$gpagen   <- gpagen(results$extremes$hb$rawData)
  results$extremes$hn$gpagen   <- gpagen(results$extremes$hn$rawData)
  results$extremes$pb$gpagen   <- gpagen(results$extremes$pb$rawData)
  results$extremes$pn$gpagen   <- gpagen(results$extremes$pn$rawData)
  
  # _______________________________________________________________________________________________
  # Create geomorph data frames
  
  # Use dense semilandmarks
  results$dense$hpbn$gdf <- geomorph.data.frame(results$dense$hpbn$gpagen)
  results$dense$hpb$gdf  <- geomorph.data.frame(results$dense$hpb$gpagen)
  results$dense$hpn$gdf  <- geomorph.data.frame(results$dense$hpn$gpagen)
  results$dense$hbn$gdf  <- geomorph.data.frame(results$dense$hbn$gpagen)
  results$dense$pbn$gdf  <- geomorph.data.frame(results$dense$pbn$gpagen)
  results$dense$hb$gdf   <- geomorph.data.frame(results$dense$hb$gpagen)
  results$dense$hn$gdf   <- geomorph.data.frame(results$dense$hn$gpagen)
  results$dense$pb$gdf   <- geomorph.data.frame(results$dense$pb$gpagen)
  results$dense$pn$gdf   <- geomorph.data.frame(results$dense$pn$gpagen)
  
  # Use  semilandmarks
  results$sparse$hpbn$gdf <- geomorph.data.frame(results$sparse$hpbn$gpagen)
  results$sparse$hpb$gdf  <- geomorph.data.frame(results$sparse$hpb$gpagen)
  results$sparse$hpn$gdf  <- geomorph.data.frame(results$sparse$hpn$gpagen)
  results$sparse$hbn$gdf  <- geomorph.data.frame(results$sparse$hbn$gpagen)
  results$sparse$pbn$gdf  <- geomorph.data.frame(results$sparse$pbn$gpagen)
  results$sparse$hb$gdf   <- geomorph.data.frame(results$sparse$hb$gpagen)
  results$sparse$hn$gdf   <- geomorph.data.frame(results$sparse$hn$gpagen)
  results$sparse$pb$gdf   <- geomorph.data.frame(results$sparse$pb$gpagen)
  results$sparse$pn$gdf   <- geomorph.data.frame(results$sparse$pn$gpagen)

  # Use extremes semilandmarks
  results$extremes$hpbn$gdf <- geomorph.data.frame(results$extremes$hpbn$gpagen)
  results$extremes$hpb$gdf  <- geomorph.data.frame(results$extremes$hpb$gpagen)
  results$extremes$hpn$gdf  <- geomorph.data.frame(results$extremes$hpn$gpagen)
  results$extremes$hbn$gdf  <- geomorph.data.frame(results$extremes$hbn$gpagen)
  results$extremes$pbn$gdf  <- geomorph.data.frame(results$extremes$pbn$gpagen)
  results$extremes$hb$gdf   <- geomorph.data.frame(results$extremes$hb$gpagen)
  results$extremes$hn$gdf   <- geomorph.data.frame(results$extremes$hn$gpagen)
  results$extremes$pb$gdf   <- geomorph.data.frame(results$extremes$pb$gpagen)
  results$extremes$pn$gdf   <- geomorph.data.frame(results$extremes$pn$gpagen)
  
  # _______________________________________________________________________________________________
  # Perform regression of shape on centroid size

  # Use dense semilandmarks
  results$dense$hpbn$regression <- procD.lm(coords ~ Csize, data = results$dense$hpbn$gdf,iter = 9999, RRPP = FALSE)
  results$dense$hpb$regression  <- procD.lm(coords ~ Csize, data = results$dense$hpb$gdf, iter = 9999, RRPP = FALSE)
  results$dense$hpn$regression  <- procD.lm(coords ~ Csize, data = results$dense$hpn$gdf, iter = 9999, RRPP = FALSE)
  results$dense$hbn$regression  <- procD.lm(coords ~ Csize, data = results$dense$hbn$gdf, iter = 9999, RRPP = FALSE)
  results$dense$pbn$regression  <- procD.lm(coords ~ Csize, data = results$dense$pbn$gdf, iter = 9999, RRPP = FALSE)
  results$dense$hb$regression   <- procD.lm(coords ~ Csize, data = results$dense$hb$gdf,  iter = 9999, RRPP = FALSE)
  results$dense$hn$regression   <- procD.lm(coords ~ Csize, data = results$dense$hn$gdf,  iter = 9999, RRPP = FALSE)
  results$dense$pb$regression   <- procD.lm(coords ~ Csize, data = results$dense$pb$gdf,  iter = 9999, RRPP = FALSE)
  results$dense$pn$regression   <- procD.lm(coords ~ Csize, data = results$dense$pn$gdf,  iter = 9999, RRPP = FALSE)

  # Use  semilandmarks
  results$sparse$hpbn$regression <- procD.lm(coords ~ Csize, data = results$sparse$hpbn$gdf,iter = 9999, RRPP = FALSE)
  results$sparse$hpb$regression  <- procD.lm(coords ~ Csize, data = results$sparse$hpb$gdf, iter = 9999, RRPP = FALSE)
  results$sparse$hpn$regression  <- procD.lm(coords ~ Csize, data = results$sparse$hpn$gdf, iter = 9999, RRPP = FALSE)
  results$sparse$hbn$regression  <- procD.lm(coords ~ Csize, data = results$sparse$hbn$gdf, iter = 9999, RRPP = FALSE)
  results$sparse$pbn$regression  <- procD.lm(coords ~ Csize, data = results$sparse$pbn$gdf, iter = 9999, RRPP = FALSE)
  results$sparse$hb$regression   <- procD.lm(coords ~ Csize, data = results$sparse$hb$gdf,  iter = 9999, RRPP = FALSE)
  results$sparse$hn$regression   <- procD.lm(coords ~ Csize, data = results$sparse$hn$gdf,  iter = 9999, RRPP = FALSE)
  results$sparse$pb$regression   <- procD.lm(coords ~ Csize, data = results$sparse$pb$gdf,  iter = 9999, RRPP = FALSE)
  results$sparse$pn$regression   <- procD.lm(coords ~ Csize, data = results$sparse$pn$gdf,  iter = 9999, RRPP = FALSE)

  # Use extremes semilandmarks
  results$extremes$hpbn$regression <- procD.lm(coords ~ Csize, data = results$extremes$hpbn$gdf,iter = 9999, RRPP = FALSE)
  results$extremes$hpb$regression  <- procD.lm(coords ~ Csize, data = results$extremes$hpb$gdf, iter = 9999, RRPP = FALSE)
  results$extremes$hpn$regression  <- procD.lm(coords ~ Csize, data = results$extremes$hpn$gdf, iter = 9999, RRPP = FALSE)
  results$extremes$hbn$regression  <- procD.lm(coords ~ Csize, data = results$extremes$hbn$gdf, iter = 9999, RRPP = FALSE)
  results$extremes$pbn$regression  <- procD.lm(coords ~ Csize, data = results$extremes$pbn$gdf, iter = 9999, RRPP = FALSE)
  results$extremes$hb$regression   <- procD.lm(coords ~ Csize, data = results$extremes$hb$gdf,  iter = 9999, RRPP = FALSE)
  results$extremes$hn$regression   <- procD.lm(coords ~ Csize, data = results$extremes$hn$gdf,  iter = 9999, RRPP = FALSE)
  results$extremes$pb$regression   <- procD.lm(coords ~ Csize, data = results$extremes$pb$gdf,  iter = 9999, RRPP = FALSE)
  results$extremes$pn$regression   <- procD.lm(coords ~ Csize, data = results$extremes$pn$gdf,  iter = 9999, RRPP = FALSE)
    
  return(results)
  
}