# wxOntogenyRegression
# This functions computes a regression of brain and neurocranial
# on age of the specimens
# José Luis Alatorre Warren, University of Zurich, 2018

# Important note:
# Currently, the dataset tableOfSpecimens assumes the age of the following specimens:
# K003, angie, bart, billie, cissie, cybil, evelyne, and jake. Currently, the age of 
# these specimens is not known and, therefore, changes should be made in this function
# so that the correlations are computed solely on the individuals with known age.

wxOntogenyRegression <- function(rawData,range,tableOfSpecimens) {
  
  # Create list with results
  results <- list()
  
  # Create subsets (species and modules) with raw data
  results$hpbn$rawData <- rawData
  results$hpb$rawData  <- rawData[ranges$dense$rl$all$B,,]
  results$hpn$rawData  <- rawData[ranges$dense$rl$all$N,,]
  results$hbn$rawData  <- rawData[                     ,,ranges$dense$rl$all$H]
  results$pbn$rawData  <- rawData[                     ,,ranges$dense$rl$all$P]
  results$hb$rawData   <- rawData[ranges$dense$rl$all$B,,ranges$dense$rl$all$H]
  results$hn$rawData   <- rawData[ranges$dense$rl$all$N,,ranges$dense$rl$all$H]
  results$pb$rawData   <- rawData[ranges$dense$rl$all$B,,ranges$dense$rl$all$P]
  results$pn$rawData   <- rawData[ranges$dense$rl$all$N,,ranges$dense$rl$all$P]

  # Perform procrustes analysis as implemented in geormorph's gpagen
  results$hpbn$gpagen <- gpagen(results$hpbn$rawData)
  results$hpb$gpagen  <- gpagen(results$hpb$rawData)
  results$hpn$gpagen  <- gpagen(results$hpn$rawData)
  results$hbn$gpagen  <- gpagen(results$hbn$rawData)
  results$pbn$gpagen  <- gpagen(results$pbn$rawData)
  results$hb$gpagen   <- gpagen(results$hb$rawData)
  results$hn$gpagen   <- gpagen(results$hn$rawData)
  results$pb$gpagen   <- gpagen(results$pb$rawData)
  results$pn$gpagen   <- gpagen(results$pn$rawData)
  
  # Create geomorph data frames
  results$hpbn$gdf <- geomorph.data.frame(results$hpbn$gpagen)
  results$hpb$gdf  <- geomorph.data.frame(results$hpb$gpagen)
  results$hpn$gdf  <- geomorph.data.frame(results$hpn$gpagen)
  results$hbn$gdf  <- geomorph.data.frame(results$hbn$gpagen)
  results$pbn$gdf  <- geomorph.data.frame(results$pbn$gpagen)
  results$hb$gdf   <- geomorph.data.frame(results$hb$gpagen)
  results$hn$gdf   <- geomorph.data.frame(results$hn$gpagen)
  results$pb$gdf   <- geomorph.data.frame(results$pb$gpagen)
  results$pn$gdf   <- geomorph.data.frame(results$pn$gpagen)
  
  results$hpbn$gdf$age <- tableOfSpecimens$Age
  results$hpb$gdf$age  <- tableOfSpecimens$Age
  results$hpn$gdf$age  <- tableOfSpecimens$Age
  results$hbn$gdf$age  <- tableOfSpecimens$Age[ranges$dense$rl$all$H]
  results$pbn$gdf$age  <- tableOfSpecimens$Age[ranges$dense$rl$all$P]
  results$hb$gdf$age   <- tableOfSpecimens$Age[ranges$dense$rl$all$H]
  results$hn$gdf$age   <- tableOfSpecimens$Age[ranges$dense$rl$all$H]
  results$pb$gdf$age   <- tableOfSpecimens$Age[ranges$dense$rl$all$P]
  results$pn$gdf$age   <- tableOfSpecimens$Age[ranges$dense$rl$all$P]

  # Perform regression of shape on centroid size
  results$hpbn$regression <- procD.lm(coords ~ age, data = results$hpbn$gdf,iter = 999, RRPP = FALSE)
  results$hpb$regression  <- procD.lm(coords ~ age, data = results$hpb$gdf, iter = 999, RRPP = FALSE)
  results$hpn$regression  <- procD.lm(coords ~ age, data = results$hpn$gdf, iter = 999, RRPP = FALSE)
  results$hbn$regression  <- procD.lm(coords ~ age, data = results$hbn$gdf, iter = 999, RRPP = FALSE)
  results$pbn$regression  <- procD.lm(coords ~ age, data = results$pbn$gdf, iter = 999, RRPP = FALSE)
  results$hb$regression   <- procD.lm(coords ~ age, data = results$hb$gdf,  iter = 999, RRPP = FALSE)
  results$hn$regression   <- procD.lm(coords ~ age, data = results$hn$gdf,  iter = 999, RRPP = FALSE)
  results$pb$regression   <- procD.lm(coords ~ age, data = results$pb$gdf,  iter = 999, RRPP = FALSE)
  results$pn$regression   <- procD.lm(coords ~ age, data = results$pn$gdf,  iter = 999, RRPP = FALSE)
  
  return(results)
  
}