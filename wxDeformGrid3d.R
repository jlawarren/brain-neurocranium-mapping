# wxDifferencePlot
# This function is based on deformGrid3d.R but with some tiny modifications
# Differences compared to deformGrid3d.R:
# + wxDifferencePlot allows to set the width of the lines
# + wxDifferencePlot asks for fast rendering when there are more than 1000000 points (instead of 1000)

wxDeformGrid3d <- function (matrix, tarmatrix, ngrid = 0, lwd = 1, showaxis = c(1, 
    2), show = c(1, 2), lines = TRUE, lcol = 1, add = FALSE, 
    col1 = 2, col2 = 3, type = c("s", "p"), size = NULL, pcaxis = FALSE) 
{
    if (inherits(matrix, "mesh3d")) 
        matrix <- vert2points(matrix)
    if (inherits(tarmatrix, "mesh3d")) 
        tarmatrix <- vert2points(tarmatrix)
    type <- type[1]
    if (dim(matrix)[1] > 1000000 && type == "s" && (is.null(size) || 
        size > 0)) {
        answer <- readline("You have a lot of landmarks\n Render them as points (faster)? (yes/NO)\n")
        if (!substr(answer, 1L, 1L) %in% c("n", "N")) 
            type <- "p"
    }
    out3d <- spheres3d
    if (type == "p") 
        out3d <- points3d
    if (!add) 
        open3d()
    k <- dim(matrix)[1]
    if (is.null(size)) {
        if (type != "p") 
            sz <- (cSize(matrix)/sqrt(k) * (1/80))
        else sz <- 10
    }
    else {
        sz <- size
    }
    if (is.numeric(sz) && sz > 0) 
        if (1 %in% show) 
            out3d(matrix, col = col1, radius = sz, size = sz)
    if (2 %in% show) 
        if (is.numeric(sz) && sz > 0) 
            out3d(tarmatrix, col = col2, radius = sz, size = sz)
    if (lines) {
        linemesh <- list()
        linemesh$vb <- t(cbind(rbind(matrix, tarmatrix), 1))
        linemesh$it <- t(cbind(1:k, 1:k, (1:k) + k))
        class(linemesh) <- "mesh3d"
        wire3d(linemesh, lwd = lwd, col = lcol, lit = FALSE)
    }
    x2 <- x1 <- x3 <- c(0:(ngrid - 1))/ngrid
    x0 <- as.matrix(expand.grid(x1, x2, x3))
    cent.mat <- scale(matrix, scale = FALSE)
    mean.mat <- colMeans(matrix)
    if (ngrid > 1) {
        xrange <- diff(range(matrix[, 1]))
        yrange <- diff(range(matrix[, 2]))
        zrange <- diff(range(matrix[, 3]))
        xrange1 <- diff(range(tarmatrix[, 1]))
        yrange1 <- diff(range(tarmatrix[, 2]))
        zrange1 <- diff(range(tarmatrix[, 3]))
        maxi <- max(c(xrange, yrange, zrange, xrange1, yrange1, 
            zrange1))
        maxi <- 1.2 * maxi
        x0 <- maxi * x0
        x0 <- scale(x0, scale = FALSE)
        if (pcaxis) 
            space <- eigen(crossprod(cent.mat))$vectors
        else space <- diag(3)
        x0 <- t(t(x0 %*% space) + mean.mat)
        x0 <- tps3d(x0, matrix, tarmatrix, lambda = 1e-08)
        outmesh <- list(vb = rbind(t(x0), 1))
        class(outmesh) <- "mesh3d"
        yinit <- xinit <- zinit <- NULL
        if (2 %in% showaxis) {
            xinit0 <- xinit <- (c(1, 2, 2 + ngrid, 1 + ngrid))
            for (i in 1:(ngrid - 2)) xinit <- cbind(xinit, (xinit0 + 
                i))
            xinit0 <- xinit
            for (i in 1:(ngrid - 2)) xinit <- cbind(xinit, xinit0 + 
                (i * ngrid))
            xinit0 <- xinit
            for (i in 1:(ngrid - 1)) xinit <- cbind(xinit, xinit0 + 
                (i * ngrid^2))
        }
        if (1 %in% showaxis) {
            yinit0 <- yinit <- c(ngrid, ngrid + ngrid^2, 2 * 
                ngrid + ngrid^2, 2 * ngrid)
            for (i in 1:(ngrid - 2)) yinit <- cbind(yinit, yinit0 + 
                i * ngrid)
            yinit0 <- yinit
            for (i in 1:(ngrid - 2)) yinit <- cbind(yinit, yinit0 + 
                i * ngrid^2)
            yinit0 <- yinit
            for (i in 1:(ngrid - 1)) yinit <- cbind(yinit, yinit0 - 
                i)
        }
        if (3 %in% showaxis) {
            zinit0 <- zinit <- (c(2, 1, 1 + ngrid^2, 2 + ngrid^2))
            for (i in 1:(ngrid - 2)) zinit <- cbind(zinit, (zinit0 + 
                i))
            zinit0 <- zinit
            for (i in 1:(ngrid - 2)) zinit <- cbind(zinit, zinit0 + 
                (i * ngrid^2))
            zinit0 <- zinit
            for (i in 1:(ngrid - 1)) zinit <- cbind(zinit, zinit0 + 
                (i * ngrid))
        }
        outmesh$ib <- cbind(xinit, yinit, zinit)
        wire3d(outmesh, lit = F)
    }
}