# wxDisplayMultilevelDistanceSets
# This functions displays the distance between 2 sets of landmarks using
# multiple colors to denote corresponding distance (level) sets
# José Luis Alatorre Warren, University of Zurich, 2018

wxDisplayMultilevelDistanceSets <- function(setA,
                                            setB,
                                            phantomA,
                                            phantomB,
                                            levelOption,
                                            colorOption,
                                            addOption) {

  # Add libraries
  library(viridis)
    
  # Assign colors from the viridis color map
  if (colorOption == "viridis"){
    colorMap <- viridis(n=16,alpha=1,begin=0.40,end=0.80,option="D")
  }
  
  if (levelOption == "uniform16"){

    # Create main sets
    
    set00 <- setA
    set16 <- setB
    
    # Create intermediate sets
    
    set08 <- array(dim = c(dim(set00),2))
    set08[,,1] <- set00
    set08[,,2] <- set16
    set08 <- arrMean3(set08)
    
    set04 <- array(dim = c(dim(set00),2))
    set04[,,1] <- set00
    set04[,,2] <- set08
    set04 <- arrMean3(set04)
    
    set12 <- array(dim = c(dim(set00),2))
    set12[,,1] <- set08
    set12[,,2] <- set16
    set12 <- arrMean3(set12)
    
    set02 <- array(dim = c(dim(set00),2))
    set02[,,1] <- set00
    set02[,,2] <- set04
    set02 <- arrMean3(set02)
    
    set06 <- array(dim = c(dim(set00),2))
    set06[,,1] <- set04
    set06[,,2] <- set08
    set06 <- arrMean3(set06)
    
    set10 <- array(dim = c(dim(set00),2))
    set10[,,1] <- set08
    set10[,,2] <- set12
    set10 <- arrMean3(set10)
    
    set14 <- array(dim = c(dim(set00),2))
    set14[,,1] <- set12
    set14[,,2] <- set16
    set14 <- arrMean3(set14)
    
    set01 <- array(dim = c(dim(set00),2))
    set01[,,1] <- set00
    set01[,,2] <- set02
    set01 <- arrMean3(set01)
    
    set03 <- array(dim = c(dim(set00),2))
    set03[,,1] <- set02
    set03[,,2] <- set04
    set03 <- arrMean3(set03)
    
    set05 <- array(dim = c(dim(set00),2))
    set05[,,1] <- set04
    set05[,,2] <- set06
    set05 <- arrMean3(set05)
    
    set07 <- array(dim = c(dim(set00),2))
    set07[,,1] <- set06
    set07[,,2] <- set08
    set07 <- arrMean3(set07)
    
    set09 <- array(dim = c(dim(set00),2))
    set09[,,1] <- set08
    set09[,,2] <- set10
    set09 <- arrMean3(set09)
    
    set11 <- array(dim = c(dim(set00),2))
    set11[,,1] <- set10
    set11[,,2] <- set12
    set11 <- arrMean3(set11)
    
    set13 <- array(dim = c(dim(set00),2))
    set13[,,1] <- set12
    set13[,,2] <- set14
    set13 <- arrMean3(set13)
    
    set15 <- array(dim = c(dim(set00),2))
    set15[,,1] <- set14
    set15[,,2] <- set16
    set15 <- arrMean3(set15)    
        
  }
  
  # Display multilevel distance sets

  wxDeformGrid3d(set00,
                 set01,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[1], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set01,
                 set02,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[2], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set02,
                 set03,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[3], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set03,
                 set04,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[4], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set04,
                 set05,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[5], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set05,
                 set06,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[6], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set06,
                 set07,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[7], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set07,
                 set08,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[8], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set08,
                 set09,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[9], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set09,
                 set10,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[10], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set10,
                 set11,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[11], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set11,
                 set12,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[12], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set12,
                 set13,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[13], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set13,
                 set14,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[14], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set14,
                 set15,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[15], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set15,
                 set16,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorMap[16], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(set00,
                 set16,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = FALSE, lcol = colorMap[16], add = TRUE, col1 = colorMap[1],
                 col2 = colorMap[16], type = c("s", "p"), size = 0.40, pcaxis = FALSE)
  
  # Add phantom
  # Here, a phantom is a set that will not be visualized because it will be described by
  # extremely small spheres. This is done to get the right 3D perspective.
  
  wxDeformGrid3d(phantomA,
                 phantomB,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = FALSE, lcol = "white", add = TRUE, col1 = colorMap[1],
                 col2 = "white", type = c("s", "p"), size = 0.0001, pcaxis = FALSE)
  
  # Return 0
  zero <- 0
  return(zero)
  
}