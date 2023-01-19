# wxPlotTangentSpace
# This function is based on plotTangentSpace but with some minor modifications

wxPlotTangentSpace <- function (A, axis1 = 1, axis2 = 2, warpgrids = TRUE, mesh = NULL, 
          label = NULL, groups = NULL, verbose = FALSE) 
{
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').")
  }
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if (!is.null(groups)) {
    groups <- factor(groups)
  }
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  pc.res <- prcomp(x)
  pcdata <- pc.res$x
  if (warpgrids == FALSE) {
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21, 
         bg = "black", cex = 2, xlab = paste("PC ", axis1), 
         ylab = paste("PC ", axis2))
    if (!is.null(groups)) {
      points(pcdata[, axis1], pcdata[, axis2], pch = 21, 
             bg = groups, cex = 2)
    }
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 
             0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), 
             lty = 2, lwd = 1)
    if (length(label != 0)) {
      text(pcdata[, axis1], pcdata[, axis2], seq(1, n), 
           adj = c(-0.7, -0.7))
    }
  }
  pcaxis.min.1 <- min(pcdata[, axis1])
  pcaxis.max.1 <- max(pcdata[, axis1])
  pc.min.1 <- pc.max.1 <- rep(0, dim(pcdata)[2])
  pc.min.1[axis1] <- pcaxis.min.1
  pc.max.1[axis1] <- pcaxis.max.1
  shape.min.1 <- arrayspecs(as.matrix(pc.min.1 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  shape.max.1 <- arrayspecs(as.matrix(pc.max.1 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  pcaxis.min.2 <- min(pcdata[, axis2])
  pcaxis.max.2 <- max(pcdata[, axis2])
  pc.min.2 <- pc.max.2 <- rep(0, dim(pcdata)[2])
  pc.min.2[axis2] <- pcaxis.min.2
  pc.max.2[axis2] <- pcaxis.max.2
  shape.min.2 <- arrayspecs(as.matrix(pc.min.2 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  shape.max.2 <- arrayspecs(as.matrix(pc.max.2 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  shapes <- list(shape.min.1, shape.max.1, shape.min.2, shape.max.2)
  names(shapes) <- c(paste("PC", axis1, "min", sep = ""), 
                     paste("PC", axis1, "max", sep = ""), paste("PC", axis2, 
                                                                "min", sep = ""), paste("PC", axis2, "max", sep = ""))
  if (warpgrids == TRUE) {
    if (k == 2) {
      layout(t(matrix(c(2, 1, 1, 1, 1, 1, 1, 1, 3), 3, 
                      3)))
    }
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21, 
         bg = "black", cex = 2, xlab = paste("PC ", axis1), 
         ylab = paste("PC ", axis2))
    if (!is.null(groups)) {
      points(pcdata[, axis1], pcdata[, axis2], pch = 21, 
             bg = groups, cex = 2)
    }
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 
             0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), 
             lty = 2, lwd = 1)
    if (length(label != 0)) {
      text(pcdata[, axis1], pcdata[, axis2], seq(1, n), 
           adj = c(-0.7, -0.7))
    }
    if (k == 2) {
      arrows(min(pcdata[, axis1]), (0.7 * max(pcdata[, 
                                                     axis2])), min(pcdata[, axis1]), 0, length = 0.1, 
             lwd = 2)
      arrows(max(pcdata[, axis1]), (0.7 * min(pcdata[, 
                                                     axis2])), max(pcdata[, axis1]), 0, length = 0.1, 
             lwd = 2)
      tps(ref, shape.min.1, 20)
      tps(ref, shape.max.1, 20)
    }
    if (k == 3) {
      if (is.null(mesh) == TRUE) {
        open3d()
        plot3d(shape.min.1, type = "s", col = "orange", 
               main = paste("PC ", axis1, " negative"), size = 0.25, 
               aspect = FALSE)
        open3d()
        plot3d(shape.max.1, type = "s", col = "orange", 
               main = paste("PC ", axis1, " positive"), size = 0.25, 
               aspect = FALSE)
      }
      if (is.null(mesh) == FALSE) {
        print("Warping mesh to axis 1 minima and maxima...")
        plotRefToTarget(ref, shape.min.1, mesh, method = "surface")
        title3d(main = paste("PC ", axis1, " negative"))
        plotRefToTarget(ref, shape.max.1, mesh, method = "surface")
        title3d(main = paste("PC ", axis1, " positive"))
      }
    }
  }
  layout(1)
  if (verbose == TRUE) {
    return(list(pc.summary = summary(pc.res), pc.scores = pcdata, 
                pc.shapes = shapes))
  }
  if (verbose == FALSE) {
    return(pc.summary = summary(pc.res))
  }
}