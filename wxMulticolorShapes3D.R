# wxMulticolorShapes3D
# This function displays distances between two 3D shapes using multicolor schemes

# José Luis Alatorre Warren, University of Zurich, 2018

wxMulticolorShapes3D <- function(shape00,
                                 shape16,
                                 colorType,
                                 perspectiveType,
                                 phantomNeurocranium01,
                                 phantomNeurocranium02) {

  # Get color scheme
  
  if (colorType == "RedBlue") {
    colorScheme <- c("#a50026","#d73027","#f46d43","#fdae61",
                     "#ffffbf",
                     "#abd9e9","#74add1","#4575b4","#313695")
  } else if (colorType == "PurpleAcqua") {
    colorScheme <- c("#40004b","#762a83","#9970ab","#c2a5cf",
                     "#ffffbf",
                     "#80cdc1","#35978f","#01665e","#003c30")
  } else if (colorType == "PurpleGreen") {
    colorScheme <- c("#40004b","#762a83","#9970ab","#c2a5cf",
                     "#f7f7f7",
                     "#a6dba0","#5aae61","#1b7837","#00441b")
  } else {  
    colorScheme <- c("#a50026","#d73027","#f46d43","#fdae61",
                     "#ffffbf",
                     "#abd9e9","#74add1","#4575b4","#313695")
  }
  
  # Load RGL 3D perspectives
  pathPerspectiveRGL <- "E:/DoctorWorks/RScripts/EvoBioPhD/"
  perspectiveRGL01 <- readRDS(paste(pathPerspectiveRGL,"perspectiveRGL01.rds",sep=""))
  perspectiveRGL02 <- readRDS(paste(pathPerspectiveRGL,"perspectiveRGL02.rds",sep=""))
  
  # Get RGL 3D perspective
  if (perspectiveType == "perspectiveLeft") {
    perspectiveRGL <- perspectiveRGL01
  } else if (perspectiveType == "perspectiveSuperior") {
    perspectiveRGL <- perspectiveRGL02
  } else {  
    perspectiveRGL <- perspectiveRGL01
  }  
  
  # Create intermediate shapes

  shape08 <- array(dim = c(dim(shape00),2))
  shape08[,,1] <- shape00
  shape08[,,2] <- shape16
  shape08 <- arrMean3(shape08)
  
  shape04 <- array(dim = c(dim(shape00),2))
  shape04[,,1] <- shape00
  shape04[,,2] <- shape08
  shape04 <- arrMean3(shape04)
  
  shape12 <- array(dim = c(dim(shape00),2))
  shape12[,,1] <- shape08
  shape12[,,2] <- shape16
  shape12 <- arrMean3(shape12)
  
  shape02 <- array(dim = c(dim(shape00),2))
  shape02[,,1] <- shape00
  shape02[,,2] <- shape04
  shape02 <- arrMean3(shape02)
  
  shape06 <- array(dim = c(dim(shape00),2))
  shape06[,,1] <- shape04
  shape06[,,2] <- shape08
  shape06 <- arrMean3(shape06)
  
  shape10 <- array(dim = c(dim(shape00),2))
  shape10[,,1] <- shape08
  shape10[,,2] <- shape12
  shape10 <- arrMean3(shape10)
  
  shape14 <- array(dim = c(dim(shape00),2))
  shape14[,,1] <- shape12
  shape14[,,2] <- shape16
  shape14 <- arrMean3(shape14)
  
  shape01 <- array(dim = c(dim(shape00),2))
  shape01[,,1] <- shape00
  shape01[,,2] <- shape02
  shape01 <- arrMean3(shape01)
  
  shape03 <- array(dim = c(dim(shape00),2))
  shape03[,,1] <- shape02
  shape03[,,2] <- shape04
  shape03 <- arrMean3(shape03)
  
  shape05 <- array(dim = c(dim(shape00),2))
  shape05[,,1] <- shape04
  shape05[,,2] <- shape06
  shape05 <- arrMean3(shape05)
  
  shape07 <- array(dim = c(dim(shape00),2))
  shape07[,,1] <- shape06
  shape07[,,2] <- shape08
  shape07 <- arrMean3(shape07)
  
  shape09 <- array(dim = c(dim(shape00),2))
  shape09[,,1] <- shape08
  shape09[,,2] <- shape10
  shape09 <- arrMean3(shape09)
  
  shape11 <- array(dim = c(dim(shape00),2))
  shape11[,,1] <- shape10
  shape11[,,2] <- shape12
  shape11 <- arrMean3(shape11)
  
  shape13 <- array(dim = c(dim(shape00),2))
  shape13[,,1] <- shape12
  shape13[,,2] <- shape14
  shape13 <- arrMean3(shape13)
  
  shape15 <- array(dim = c(dim(shape00),2))
  shape15[,,1] <- shape14
  shape15[,,2] <- shape16
  shape15 <- arrMean3(shape15)
  
  shape15 <- array(dim = c(dim(shape00),2))
  shape15[,,1] <- shape14
  shape15[,,2] <- shape16
  shape15 <- arrMean3(shape15)
  
  shape78 <- array(dim = c(dim(shape00),2))
  shape78[,,1] <- shape07
  shape78[,,2] <- shape08
  shape78 <- arrMean3(shape78)
  
  shape89 <- array(dim = c(dim(shape00),2))
  shape89[,,1] <- shape08
  shape89[,,2] <- shape09
  shape89 <- arrMean3(shape89)
  
  # Display the distance lines connecting the 3D shapes 
  wxDeformGrid3d(shape00,
                 shape02,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[1], add = FALSE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape02,
                 shape04,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[2], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape04,
                 shape06,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[3], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape06,
                 shape07,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[4], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape07,
                 shape78,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[4], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape78,
                 shape08,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[5], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape08,
                 shape89,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[5], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape89,
                 shape09,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[6], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape09,
                 shape10,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[6], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape10,
                 shape12,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[7], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape12,
                 shape14,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[8], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  wxDeformGrid3d(shape14,
                 shape16,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = TRUE, lcol = colorScheme[9], add = TRUE, col1 = "yellow",
                 col2 = "red", type = c("s", "p"), size = 0.00, pcaxis = FALSE)
  
  # Display the spheres comprising the 3D shapes   
  wxDeformGrid3d(shape00,
                 shape16,
                 ngrid = 0, lwd = 7.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = FALSE, lcol = "orange", add = TRUE, col1 = colorScheme[1],
                 col2 = colorScheme[9], type = c("s", "p"), size = 0.40, pcaxis = FALSE)
  
  # Display neurocranium phantom
  # This is needed to get the correct RGL 3D perspective
  wxDeformGrid3d(phantomNeurocranium01,
                 phantomNeurocranium02,
                 ngrid = 0, lwd = 14.0, showaxis = c(1, 2),
                 show = c(1, 2), lines = FALSE, lcol = "orange", add = TRUE, col2 = "#a50026",
                 col1 = "#313695", type = c("s", "p"), size = 0.0001, pcaxis = FALSE)
  
  # Apply RGL 3D perspective: lateral view (perspectiveRGL)
  rgl.viewpoint(userMatrix = perspectiveRGL$userMatrix, fov = perspectiveRGL$FOV, zoom = perspectiveRGL$zoom)
  par3d(windowRect = c(0, 23, 3840, 2140))

  return(colorScheme)
  
}