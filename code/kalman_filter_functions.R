#######
# This file contains functions for smoothing GPS data


library('sp')
library('rgdal')
library('OpenStreetMap')


shift_vec <- function(vec, shift) {
    # Shift vector
    #
    # Args:
    #   vec: vector to shift.
    #   shift: integer.
    #
    # Returns:
    #   shifted vector.
    if (length(vec) <= abs(shift)) {
        shift.vec <- rep(NA, length(vec))
    } else{
        if (shift > 0) {
            shift.vec <- c(rep(vec[1], shift), vec[1:(length(vec)-shift)])
            for (i in 1:shift) {
                shift.vec[i] <- NA
            }
        } else {
            shift.vec <- c(vec[(abs(shift)+1):length(vec)], rep(NA, abs(shift)))
        }
    }
    return(shift.vec)
}


longlatToUTM <- function(long, lat, zone=14) {
    # Convert longitude and latitude to UTM
    #
    # Args:
    #   long: number.
    #   lat: number.
    #   zone: UTM zone. The default is 14.
    #
    # Returns:
    #   UTM position.
    longlat <- data.frame(long=long, lat=lat)
    coordinates(longlat) <- c('long', 'lat')
    proj4string(longlat) <- CRS('+proj=longlat +datum=WGS84')
    utm <- spTransform(longlat, CRS(paste('+proj=utm +zone=', zone, '
                                          ellps=WGS84', sep='')))
    return(as.data.frame(utm))
}

UTMTolonglat <- function(utm, zone=14) {
    # Convert UTM to longitude and latitude.
    #
    # Args:
    #   utm: dataframe including UTM position.
    #   zone: UTM zone.
    #
    # Returns:
    #   data frame including longitude and latitude.
    utm.sp <- SpatialPoints(utm, proj4string=CRS(paste('+proj=utm +zone=',
                                                       zone, ' ellps=WGS84',
                                                       sep='')))
    longlat <- spTransform(utm.sp, CRS("+proj=longlat +datum=WGS84"))
    return(as.data.frame(longlat))
}


permutation_matrix <- function(n, N) {
    # Calculate permutation matrix.
    # Permutation matrix can convert the direct sum to a vector.
    #
    # Args:
    #   n: number. The number of rows of a matrix to be tranformed.
    #   N: number. The number of elements in direct sum.
    #
    # Returns:
    #   Permutation matrix.
    permutation.matrix <- matrix(0, nrow=(n*N)^2, ncol=n^2)
    diag <- diag(1, n)
    for (i in 1:N) {
        for (j in 1:n) {
            sub <- matrix(0, nrow=n*N, ncol=n^2)
            sub[((i-1)*n+1):(i*n), ((j-1)*n+1):(j*n)] <- diag
            permutation.matrix[((i-1)*n*n*N + (j-1)*n*N+1):((i-1)*n*n*N+j*n*N),
                               ] <- sub
        }
    }
    return (permutation.matrix)
}


direct_sum <- function(i.matrix, n) {
    # Calculate direct sum.
    #
    # Args:
    #   i.matrix: matrix in direct sum.
    #   n: number.
    #
    # Returns:
    #   direct sum.
    sum.matrix <- kronecker(diag(1, n), i.matrix)
    return (sum.matrix)
}


states_estimation <- function(data, f.matrix, k.guess, h.matrix, positions) {
    # Estimate states in Autocovariance least-squares method.
    #
    # Args:
    #   data: matrix. The 1st column must be delta.time; the 2nd column must be
    #         longitude and the 3erd column must be latitude.
    #   f.matrix: state transition matrix in Kalman filter.
    #   k.guess: guessed Kalman filter gain matrix.
    #   h.matrix: transformation matrix in Kalman filter.
    #   positions: a list including positions of dleta.time in matrix F
    #
    # Returns:
    # states estimation.
    n <- nrow(F.matrix)

    xhat.posterior <- matrix(0, nrow=n, ncol=ncol(data))
    xhat.prior <- matrix(0, nrow=n, ncol=ncol(data))
    xhat.prior[, 1] <- c(data[-1, 1], rep(0, 2))
    for (i in 1:(ncol(data)-1)) {
        # update
        xhat.posterior[, i] <- xhat.prior[, i] + K.guess %*% (data[-1, i] -
                                                              H.matrix %*%
                                                              xhat.prior[, i])
        # update matrix F
        delta.time <- data[1, i+1]
        for (j in length(positions)) {
            F.matrix[j[1], j[2]] <- delta.time
        }
        # prediction
        xhat.prior[, i+1] <- F.matrix %*% xhat.posterior[, i]
    }
    return (xhat.prior)
}


ALS <- function(data, H.matrix, F.matrix, K.guess, N, ignore=100, positions,
                tolerance=6.88e-22) {
    # Calculate covariance matrices of process and measurement noises in Kalman
    # filter.
    #
    # Args:
    #   data : matrix including delta.time, long, and lat.
    #   H.matrix: transformation matrix in Kalman filter.
    #   F.matrix: states transition matrix in Kalman filter.
    #   K.guess: guessed Kalman filter gain matrix.
    #   N: the number of lags.
    #   ignore: data to be ignored in the beginning until initial condition is
    #           negligiable. The default is 100.
    #   positions: a list including positions of dleta.time in matrix F
    #   tolerance: parameter in R function `solve`.
    #
    # Returns:
    #   matrix Q and matrix R.
    data.t <- t(data)
    # dimensions
    n <- nrow(F.matrix)
    p <- nrow(H.matrix)

    F.line <- F.matrix - F.matrix %*% K.guess %*% H.matrix
    # calculate matrix O and matrix Gamma
    O.matrix <- matrix(0, nrow=p*N, ncol=n)
    for (i in 1:N) {
        temp.HF <- H.matrix %*% F.line^(N-1)
        O.matrix[(p*(i-1)+1):(p*i), ] <- temp.HF
    }
    Gamma.matrix <- matrix(0, nrow=p*N, ncol=n*N)
    for (i in 1:(N-1)) {
        Gamma.matrix[(p*i+1):(p*N), (n*(i-1)+1):(n*i)] <- O.matrix[1:(p*(N-i)), ]
    }
    # calculate matrix Phi
    Phi.matrix <- Gamma.matrix %*% direct_sum(-F.matrix%*%K.guess, N)
    # calculate matrix D
    D.matrix <- kronecker(O.matrix, O.matrix) %*% solve(diag(i, n*n) - kronecker(F.line, F.line)) + kronecker(Gamma.matrix, Gamma.matrix) %*% permutation_matrix(n, N)
    # calculate matrix A
    A.sub1 <- D.matrix %*% kronecker(diag(1, n), diag(1, n))
    A.sub2 <- D.matrix %*% kronecker(F.matrix%*%K.guess, F.matrix%*%K.guess) +
        (kronecker(Phi.matrix, Phi.matrix) + diag(1, p^2*N^2)) %*%
        permutation_matrix(p, N)
    A.matrix <- cbind(A.sub1, A.sub2)
    A.plus <- solve(t(A.matrix) %*% A.matrix, tol=tolerance) %*% t(A.matrix)

    # calculate bhat
    # first estimate states
    bhat <- matrix(0, nrow=p*N, ncol=p*N)
    xhat.prior <- states_estimation(data.t, F.matrix, K.guess, H.matrix, positions)
    xhat <- xhat.prior
    # calculate K-innovations
    K.innovations <- data.t[-1, (ignore+1):ncol(data.t)] - H.matrix %*% xhat[,
                                                                             (ignore+1):ncol(xhat)]
    # calculate autocorrelations
    Nd <- ncol(data.t) - ignore
    for (i in 0:(N-1)) {
        temp <- K.innovations[, 1:(ncol(K.innovations)-i)] %*%
            t(K.innovations[, (i+1):ncol(K.innovations)])
        temp <- temp / (Nd-i)
        I.matrix <- matrix(0, nrow=N, ncol=N)
        if (i == 0) {
            diag(I.matrix) <- 1
            bhat <- bhat + kronecker(I.matrix, temp)
        } else if (i == (N-1)){
            I.matrix[1, N] <- 1
            bhat <- bhat + kronecker(I.matrix, temp)
            I.matrix <- matrix(0, nrow=N, ncol=N)
            I.matrix[N, 1] <- 1
            bhat <- bhat + kronecker(I.matrix, t(temp))
        } else {
            diag(I.matrix[, -(1:i)]) <- 1
            bhat <- bhat + kronecker(I.matrix, temp)
            I.matrix <- matrix(0, nrow=N, ncol=N)
            diag(I.matrix[-(1:i), ]) <- 1
            bhat <- bhat + kronecker(I.matrix, t(temp))
        }
    }
    bhat.s <- matrix(bhat, ncol=1)
    QR.s <- A.plus %*% bhat.s
    Q.matrix <- matrix(QR.s[1:(n*n), 1], nrow=n)
    R.matrix <- matrix(QR.s[(n*n+1):nrow(QR.s), 1], nrow=p)
    return (list('Q.matrix'=Q.matrix, 'R.matrix'=R.matrix))
}


kalman_update <- function(data, Xhat.matrix, P.matrix, H.matrix, w.matrix,
                          R.matrix) {
    # Run kalman filter.
    #
    # This function conjugates with w.matrix and R.matrix estimated from
    # empirical method.
    #
    # Args:
    #   data: vector including delta.time, long, and lat.
    #   Xhat.matrix: posterior values of longitude and latitude from time t-1.
    #   P.matrix: posterior values of  covariance matrix from time t-1.
    #   H.matrix: transformation matrix
    #   w.matrix: covariance of velocities on longitude and latitude.
    #   R.matrix: covariance of measurement noise
    #
    # Returns:
    #   Posterior values of longitude and latitude from time t.
    #   Posterior values of covariance matrix from time t.

    # initalization
    Xhat.posterior <- Xhat.matrix
    P.posterior <- P.matrix
    R.matrix <- R.matrix
    H.matrix <- H.matrix
    # update matrix F and G
    delta.time <- data[1, 1]
    F.matrix <- diag(1, 4)
    F.matrix[1, 3] <- delta.time
    F.matrix[2, 4] <- delta.time
    G.matrix <- matrix(0, nrow=4, ncol=2)
    G.matrix[1, 1] <- 1/2 * delta.time^2
    G.matrix[2, 2] <- 1/2 * delta.time^2
    G.matrix[3, 1] <- delta.time
    G.matrix[4, 2] <- delta.time
    # calculate matrix Q
    Q.matrix <- G.matrix %*% w.matrix %*% t(G.matrix)

    # prediction
    Xhat.prior <- F.matrix %*% Xhat.posterior
    P.prior <- F.matrix %*% P.posterior %*% t(F.matrix) + Q.matrix
    # update
    K = P.prior %*% t(H.matrix) %*% solve(H.matrix %*% P.prior %*% t(H.matrix)
                                          + R.matrix, tol=2.79069e-18)
    Xhat.posterior <- Xhat.prior + K %*% (data[-1, ]- H.matrix %*% Xhat.prior)
    KH <- K %*% H.matrix
    I = diag(1, dim(KH)[1])
    P.posterior <- (I - KH) %*% P.prior

    return (list('Xhat.posterior'=Xhat.posterior,
                 'P.posterior'=P.posterior))
}


kalman_update2 <- function(data, Xhat.matrix, P.matrix, F.matrix, H.matrix,
                           Q.matrix, R.matrix, positions) {
    # Run kalman filter
    #
    # This function conjugates with w.matrix and R.matrix estimated from
    # ALS method.
    #
    # Args:
    #   data: vector including delta.time, long, and lat.
    #   Xhat.matrix: posterior values of longitude and latitude from time t-1.
    #   P.matrix: posterior values of  covariance matrix from time t-1.
    #   F.matrix: state transition matrix.
    #   H.matrix: transformation matrix
    #   Q.matrix: covariance of velocities on longitude and latitude.
    #   R.matrix: covariance of measurement noise
    #   positions: a list including positions of dleta.time in matrix F
    #
    # Returns:
    #   Posterior values of longitude and latitude from time t.
    #   Posterior values of covariance matrix from time t.

    # initalization
    Xhat.posterior <- Xhat.matrix
    P.posterior <- P.matrix
    Q.matrix <- Q.matrix
    R.matrix <- R.matrix
    H.matrix <- H.matrix
    # update matrix F and G
    delta.time <- data[1, 1]
    F.matrix <- diag(1, 4)
    F.matrix[1, 3] <- delta.time
    F.matrix[2, 4] <- delta.time

    # prediction
    Xhat.prior <- F.matrix %*% Xhat.posterior
    P.prior <- F.matrix %*% P.posterior %*% t(F.matrix) + Q.matrix
    # update
    K = P.prior %*% t(H.matrix) %*% solve(H.matrix %*% P.prior %*% t(H.matrix)
                                          + R.matrix, tol=2.79069e-18)
    Xhat.posterior <- Xhat.prior + K %*% (data[-1, ]- H.matrix %*% Xhat.prior)
    KH <- K %*% H.matrix
    I = diag(1, dim(KH)[1])
    P.posterior <- (I - KH) %*% P.prior

    return (list('Xhat.posterior'=Xhat.posterior,
                 'P.posterior'=P.posterior))
}


plotGPS_png <- function(data, long=2, lat=3, plot.name, map.type='skobbler',
                        width=1000, height=800, res=100, line.type='l',
                        color=scales::alpha('blue', 0.5), line.width=4) {
    # Plot longitude and latitude on actual map.
    #
    # Args:
    #   data: dataframe or matrix. There must be one column of longitude and
    #   one column of latitude.
    #   long: the index of the column containing longitude information.
    #   lat: the index of column containing latitude information.
    #   plot.name: the name of the plot to be saved.
    #   map.type: the type of map to be used.
    #   width: width of the plot.
    #   height: height of the plot.
    #   res: resolution of the plot.
    #   line.type: line type.
    #   color: line color on the plot.
    #   line.width: line width.
    #
    # Returns:
    #   saved a png.
    upperLeft <- c(max(data[, lat]), min(data[, long]))
    lowerRight <- c(min(data[, lat]), max(data[, long]))
    map <- openmap(upperLeft, lowerRight, type=map.type)

    transmap <- openproj(map, projection='+proj=longlat')

    png(plot.name, width=width, height=height, res=res)
    plot(transmap, raster=T)
    lines(data[, long], data[, lat], type=line.type, col=color, lwd=line.width)
    dev.off()
}
